#!/usr/bin/env python
# Copyright 2019 Irwin Jungreis
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""
Determine the PhyloCSF Regions most likely to be unannotated coding regions, by excluding
  ones that overlap known coding regions in the same frame or antisense frame and known
  pseudogenes in any frame, regions less than 9 codons long, and regions deemed more
  likely to be novel antisense regions than novel coding ones.
Sort the resulting Novel PhyloCSF Regions starting with those most likely to be coding.
The PhyloCSF Regions must have been previously created using the HMM.
This script requires rpy2 for the SVM step. Before invoking that the e1071 package must be
    installed, which can be done by calling ClassSVM.install_e1071 from an interactive
    session; this only needs to be done once for each machine on which this will be run.
    
Steps 1-5 need to be executed sequentially, each step using the output from the previous
and sometimes other inputs. The resulting and intermediate files go in the same
directory "homeDir". The final results are:
    Regions.04.txt: A Tab-delimited spreadsheet with one line for each Novel PhyloCSF
                    Region, in order, including the various scores associated with it.
    PhyloCSFNovel.bed: A bed-format file with the Novel PhyloCSF Regions, in order.
"""
from __future__ import division
from __future__ import print_function
import sys, shutil, bisect, random, math
from CommonUtils import (pjoin, assure_dir, err_msg, ls, myopen, bed_line_to_intervals,
                         DictClass, anti_interval, Tab)
from DelimitedFileUtil import DFR, DFW
from ClassTranscript import (bed_line_to_tr, intervals_to_bed_line, subtract_regions,
                             construct_transcript, intersect_tr_intervals)
from ClassTrOverlapChecker import TrOverlapChecker as OverlapChecker
import ClassSVM

assert sys.version_info[:2] == (2, 7), 'This is intended to run in Python 2.7.'

MinRelBranchLength = 0.1
MinCodons = 9 # Exclude regions shorter than this

RegTypes = [
    'CodingOverlap',   # Overlaps annotated coding region in same frame, but no pseudogene
    'AntisenseOverlap',# Overlaps annotated coding region in antisense frame, but no pseudogene
    'PseudoOverlap',   # Overlaps annotated pseudogene
    'NoOverlap',       # In original region list and not in previous 3 categories
    'Extension',       # Sub interval of *Overlap not overlapping annotated coding or pseudo
    ]
Fields = [
    'Rank',    # 1-Based rank after exclusions, sorted by SvmScore
    'Name',    # e.g., 'chr1:200-500-'
    'RegType', # One of RegTypes
    'Chrom',
    'Start',   # Nominal minimum
    'End',
    'Strand',
    'ScorePerCodon',
    'NumCodons',
    'Bls',
    'AntiScorePerCodon', # Of region with same 3rd codon positions
    'ScoreDiff',         # ScorePerCodon - AntiScorePerCodon
    'Parent',            # Name of parent for Extension regions, otherwise 'NA'
    'SvmScore',          # SVM prediction using NumCodons, ScorePerCodon, Bls, ScoreDiff
    'SvmAntisenseScore',
    ]

for field in Fields + RegTypes :
    globals()[field] = field # Allows referring to field names without quotes
intFields = [Rank, NumCodons, Start, End]
floatFields = [ScorePerCodon, Bls, AntiScorePerCodon, ScoreDiff,
               SvmScore, SvmAntisenseScore]

def get_reader(fileName) :
    reader = DFR(fileName, intFields = intFields, floatFields = floatFields)
    return reader

def get_input_fileName(homeDir, index) :
    return pjoin(homeDir, 'Regions.%02d.txt' % index)

def get_output_fileName(homeDir, index) :
    return pjoin(homeDir, 'Regions.%02d.txt' % index)

def main() :
    """
    First run this script with -Step1.
    Then, run PhyloCSF on Regions.pcsf.in to produce Regions.pcsf.out.
    Then run this script with -OtherSteps.
    """
    def print_usage() :
        print('Usage: (-Step1 HOMEDIR REGIONS_DIR CODING_BED_FILE PSEUDO_BED_FILE | '\
                      '-OtherSteps HOMEDIR [NUM_TRAINING])', file = sys.stderr)
        return 1
    if len(sys.argv) < 2 :
        return(print_usage())
    if sys.argv[1] == '-Step1' :
        if len(sys.argv) != 6 :
            return(print_usage())
        homeDir, regionsDir, codingBedFileName, pseudoBedFileName = sys.argv[2 : 6]
        classify_regions(homeDir, regionsDir, codingBedFileName, pseudoBedFileName)
    elif sys.argv[1] == '-OtherSteps' :
        if not 3 <= len(sys.argv) <= 4 :
            return(print_usage())
        homeDir = sys.argv[2]
        numTraining = int(sys.argv[3]) if len(sys.argv) > 3 else None
        do_other_steps(homeDir, numTraining)
    else :
        return(print_usage())

def do_other_steps(homeDir, numTraining = None) :
    add_phyloCSF_fields(homeDir)
    add_svm_scores(homeDir, numTraining)
    exclude_and_sort(homeDir)
    make_bed(homeDir)

# *********** Step 1 ***************************
def classify_regions(homeDir, regionsDir, codingBedFileName, pseudoBedFileName) :
    """
    Inputs: - PhyloCSF Regions bed files in regionsDir, produced by the PhyloCSF HMM.
                  The file names should have the format:
                      {CHROMOSOME}.Strand{STRAND}.Frame{FRAME}.coding.bed.gz
                  where STRAND is + or - and FRAME is 0, 1, or 2.
            - Coding and pseudogene annotated transcripts.
    Output: Regions.01.txt or Regions.01.txt
    Find overlapping transcripts, classify regions, and create extension regions.
    Also create input file for running PhyloCSF using the strategy=mle option.
    """
    assure_dir(homeDir)
    outFileName = get_output_fileName(homeDir, 1)
    err_msg('Reading transcripts')
    codingTrs = [bed_line_to_tr(line) for line in myopen(codingBedFileName)]
    pseudoTrs = [bed_line_to_tr(line) for line in myopen(pseudoBedFileName)]
    
    err_msg('Creating overlap checkers')
    codingOverlapChecker = OverlapChecker(codingTrs, onlyCDS = True)
    pseudoOverlapChecker = OverlapChecker(pseudoTrs, onlyCDS = False)

    outRecs = []
    lineCounter = plusCounter = minusCounter = 0

    def get_rec_name(rec) :
        return '%s:%d-%d%s' % (rec.Chrom, rec.Start, rec.End, rec.Strand)

    err_msg('Processing input PhyloCSF Regions.')
    for bedFileName in ls(regionsDir) :
        if not bedFileName.endswith('.coding.bed.gz') :
            continue
        for line in myopen(pjoin(regionsDir, bedFileName)) :
            lineCounter += 1
            if '.Strand+' in bedFileName :
                plusCounter += 1
            else :
                minusCounter += 1
            dummyName, chrom, bedIntervals, strand = bed_line_to_intervals(line)
            assert len(bedIntervals) == 1
            interval = bedIntervals[0]
            rec = DictClass()
            rec.Chrom = chrom
            rec.Start = interval[0]
            rec.End   = interval[1]
            rec.NumCodons = (rec.End - rec.Start + 1) // 3
            rec.Strand = strand
            rec.Name = get_rec_name(rec)
            rec.Parent = 'NA'
            assert (rec.End - rec.Start + 1) % 3 == 0, rec.Name

            # Find overlaps and set up RegType
            codingOverlapTrs = codingOverlapChecker.overlapping_trs(rec.Chrom, '+-', interval)
            pseudoOverlapTrs = pseudoOverlapChecker.overlapping_trs(rec.Chrom, '+-', interval)
            if len(pseudoOverlapTrs) > 0 :
                rec.RegType = PseudoOverlap
            elif any_same_frame(interval, rec.Strand, codingOverlapTrs) :
                rec.RegType = CodingOverlap
            elif any_anti_frame(interval, rec.Strand, codingOverlapTrs) :
                rec.RegType = AntisenseOverlap
            else :
                rec.RegType = NoOverlap
            outRecs.append(rec)

            # Subtract from intervals overlaps with pseudogenes in any frame or
            #     coding in same or antisense frame.
            # Make a new record for each resulting segment
            if rec.RegType != NoOverlap :
                intervals = subtract_trs(interval, rec.Strand, pseudoOverlapTrs, codingOverlapTrs)
                for interval in intervals :
                    subRec = DictClass()
                    subRec.RegType = 'Extension'
                    subRec.Chrom = rec.Chrom
                    subRec.Start = interval[0]
                    subRec.End   = interval[1]
                    subRec.NumCodons = (subRec.End - subRec.Start + 1) // 3
                    subRec.Strand = rec.Strand
                    subRec.Name = get_rec_name(subRec)
                    subRec.Parent = rec.Name
                    assert (subRec.End - subRec.Start + 1) % 3 == 0, subRec.Name
                    outRecs.append(subRec)

    outRecs.sort(key = lambda rec : (RegTypes.index(rec.RegType), rec.Chrom, rec.Strand,
                                     rec.Start, rec.End))

    fields = list(Fields)
    fields.remove(Rank) # We'll insert it at the beginning later
    outDFW = DFW(outFileName, fields)
    for recInd, rec in enumerate(outRecs) :
        outDFW.write_line(rec)
    outDFW.close()

    err_msg('Writing PhyloCSF input file.')
    _write_phyloCSF_in(homeDir)

def _write_phyloCSF_in(homeDir) :
    """
    Input:  Regions.01.txt
    Output: Regions.pcsf.in containing each region and its antisense region
    """
    inFileName = get_input_fileName(homeDir, 1)
    outFileName = pjoin(homeDir, 'Regions.pcsf.in')
    with myopen(outFileName, 'w') as outf :
        for recInd, rec in enumerate(get_reader(inFileName)) :
            print(chromInt2Str(rec.Chrom, [(rec.Start, rec.End)]),
                  rec.Strand,
                  sep = Tab,
                  file = outf)
            antiInterval, antiStrand = anti_interval([rec.Start, rec.End], rec.Strand)
            print(chromInt2Str(rec.Chrom, [antiInterval]),
                  antiStrand,
                  sep = Tab,
                  file = outf)

# *********** Step 2 ***************************
"""
- Run PhyloCSF with --bls (using the default --strategy=mle) on the multispecies
  alignments of the intervals specified in Regions.pcsf.in.
  Output file name must be Regions.pcsf.out.
  Each line of Regions.pcsf.out output file must be of the form:
        CHROM:CODON_START-CODON_END   STRAND  score(decibans) PHYLOCSF_SCORE  BLS
      where PHYLOCSF_SCORE and BLS are the output of PhyloCSF using the --strategy=mle
      and --bls options. The separator between fields must be a Tab character. For
      codons for which there is no alignment or for which PhyloCSF fails, put
      'No_Alignment' or 'failure', respectively, instead of the last 3 fields.
- Then call add_phyloCSF_fields
"""
def add_phyloCSF_fields(homeDir) :
    """
    Input:  Regions.01.txt
            Regions.pcsf.out
    Output: Regions.02.txt
    Fill fields ScorePerCodon, NumCodons, Bls,
        AntiScorePerCodon, and ScoreDiff
    """
    inFileName = get_input_fileName(homeDir, 1)
    err_msg('Reading %s.' % inFileName)
    inDFR = get_reader(inFileName)
    inRecs = list(inDFR)

    err_msg('Reading PhyloCSF scores.')
    pcsfDict = _get_pcsf_dict(homeDir)
    
    err_msg('Setting PhyloCSF-related fields.')
    _fill_mle_pcsf_fields(inRecs, pcsfDict)

    outFileName = get_output_fileName(homeDir, 2)
    err_msg('Writing %s.' % outFileName)
    outDFW = DFW(outFileName, inDFR.get_fieldNames())
    for rec in inRecs :
        outDFW.write_line(rec)
    outDFW.close()

def _get_pcsf_dict(homeDir) :
    fieldNames = ['IntervalsStr', 'Strand', 'Status', 'RawScore', 'Bls']
    fileName = pjoin(homeDir, 'Regions.pcsf.out')
    return {rec.IntervalsStr + rec.Strand : (rec.Status, rec.RawScore, rec.Bls)
            for rec in DFR(fileName, fieldNames = fieldNames,
                           floatFields = ['RawScore', 'Bls'])}

def _fill_mle_pcsf_fields(inRecs, pcsfDict) :
    """
    Fill in the PhyloCSF related fields in inRecs.
    """
    for recInd, rec in enumerate(inRecs) :
        intervalStr = chromInt2Str(rec.Chrom, [(rec.Start, rec.End)])
        status, rawScore, bls = pcsfDict[intervalStr + rec.Strand]
        rec[ScorePerCodon] = rawScore / rec.NumCodons
        rec[Bls] = bls
        antiInterval, antiStrand = anti_interval([rec.Start, rec.End], rec.Strand)
        antiIntervalStr = chromInt2Str(rec.Chrom, [antiInterval])
        status, antiScore, antiBls = pcsfDict[antiIntervalStr + antiStrand]
        if status not in ['No_Alignment', 'failure'] :
            rec[AntiScorePerCodon] = antiScore / rec.NumCodons
        scorePerCodon     = rec[ScorePerCodon]
        antiScorePerCodon = rec[AntiScorePerCodon]
        
        if scorePerCodon != '' :
            if antiScorePerCodon != '' :
                rec[ScoreDiff] = scorePerCodon-antiScorePerCodon
            else :
                """
                In rare cases original interval has alignment but 2-base extension of
                anti-interval doesn't. In those cases, keep AntiScorePerCodon='' but treat
                it as 0 for calculating ScoreDiff so that SVM has something to work with.
                """
                rec[ScoreDiff] = scorePerCodon

# *********** Step 3 ***************************
def add_svm_scores(homeDir, numTraining = None) :
    """
    Input:  Regions.02.txt
    Output: Regions.03.txt
    Add fields SvmAntisenseScore, SvmScore
    """
    if numTraining == None :
        numTraining = 10000
    inFileName = get_input_fileName(homeDir, 2)
    err_msg('Reading %s.' % inFileName)
    inDFR = get_reader(inFileName)
    inRecs = list(inDFR)
    
    err_msg('Getting training recs.')
    codingTrainingRecs, antiTrainingRecs, otherTrainingRecs = \
        _get_svm_training_recs(inRecs, numTraining)

    for features, svmField, trainingRecs in [
        ([NumCodons, ScorePerCodon, ScoreDiff],
         SvmAntisenseScore, [codingTrainingRecs, antiTrainingRecs]),
         
        ([NumCodons, ScorePerCodon, Bls, ScoreDiff],
         SvmScore, [codingTrainingRecs, otherTrainingRecs]),
         ] :
        err_msg('Making training and test vectors for %s.' % svmField)
        trainingClasses = [1] * len(trainingRecs[0]) + [0] * len(trainingRecs[1])
        trainingVecs = [[rec[feature] for feature in features]
                        for rec in trainingRecs[0] + trainingRecs[1]]
        testVecs = [[rec[feature] for feature in features]
                     for rec in inRecs]
                     
        err_msg('Training and using SVM for %s.' % svmField)
        svmProbs = ClassSVM.ClassSVM(trainingVecs, trainingClasses)(testVecs)

        err_msg('Setting %s fields.' % svmField)
        for rec, prob in zip(inRecs, svmProbs) :
            rec[svmField] = prob

    outFileName = get_output_fileName(homeDir, 3)
    err_msg('Writing %s.' % outFileName)
    outDFW = DFW(outFileName, inDFR.get_fieldNames())
    for rec in inRecs :
        outDFW.write_line(rec)
    outDFW.close()

def _get_svm_training_recs(inRecs, numTraining) :
    """
    Return numTraining recs of each of CodingOverlap, AntisenseOverlap, and NoOverlap types
        chosen randomly from inRecs, in a reproducible way.
    """
    codingRecs  = [rec for rec in inRecs if rec.RegType == CodingOverlap]
    antiRecs    = [rec for rec in inRecs if rec.RegType == AntisenseOverlap]
    otherRecs   = [rec for rec in inRecs if rec.RegType == NoOverlap]
    assert numTraining <= len(codingRecs), (numTraining, len(codingRecs))
    assert numTraining <= len(antiRecs),   (numTraining, len(antiRecs))
    assert numTraining <= len(otherRecs),  (numTraining, len(otherRecs))
    random.seed(0) # Make it reproducible
    codingTrainingRecs = random.sample(codingRecs, numTraining)
    antiTrainingRecs   = random.sample(antiRecs,   numTraining)
    otherTrainingRecs  = random.sample(otherRecs,  numTraining)
    return codingTrainingRecs, antiTrainingRecs, otherTrainingRecs

# *********** Step 4 ***************************
def exclude_and_sort(homeDir) :
    """
    Input:  Regions.03.txt
    Output: Regions.04.txt
    Exclude recs with NumCodons <= 8, because SvmScore was weird for these,
        presumably because of insufficient training vectors of that size.
    Exclude recs with SvmAntisenseScore < 0.3. That would keep around 99% of CodingOverlap
        regions, and exclude around 94% of AntisenseOverlap regions (in hg38).
    Keep only NoOverlap and Extension recs, since we only used the others for training.
    Sort by decreasing SvmScore.
    Add 1-based Rank field.
    """
    inFileName = get_input_fileName(homeDir, 3)
    err_msg('Reading %s.' % inFileName)
    inDFR = get_reader(inFileName)
    inRecs = list(inDFR)
    outFileName = get_output_fileName(homeDir, 4)

    outRecs = [rec for rec in inRecs if rec.RegType in [NoOverlap, Extension] and
                                        rec.NumCodons >= MinCodons and
                                        rec.SvmAntisenseScore >= 0.3]
    outRecs.sort(key = lambda rec : -rec.SvmScore)

    outDFW = DFW(outFileName, [Rank] + inDFR.get_fieldNames())
    for index, rec in enumerate(outRecs) :
        rec.Rank = index + 1
        outDFW.write_line(rec)
    outDFW.close()

# *********** Step 5 ***************************
def make_bed(homeDir) :
    """
    Input: Regions.04.txt
    Output: PhyloCSFNovel.bed.
    Output novel PhyloCSF regions in bed format for browser tracks.
    
    """
    err_msg('Writing .bed file.')
    inFileName = get_input_fileName(homeDir, 4)
    novelFileName = pjoin(homeDir, 'PhyloCSFNovel.bed')
    recs = list(get_reader(inFileName))
    
    """
    Color regions on +/- strands green/red to match PhyloCSF tracks, and dim ones with
       higher ranks. Ranks matter more at the start so use a logarithmicish scale. Put the
       middle of the range at the somewhat arbitrary rank 5000, which is sort of where
       they aren't as useful.
    UCSC says to limit to 8 colors to keep browser working well.
    """
    numBins = 8
    midInd = 5000
    
    colorStrs = {'+' : '0,175,0', '-' : '200,0,0'}
    
    numRecs = len(recs) # Slightly more than the largest index
    def scale_rank(recInd) :
        # 0 -> 0, midInd -> 0.5, numRecs - 1 -> 1 - epsilon
        a = (numRecs - 2 * midInd) / midInd ** 2
        return math.log(1 + a * recInd) / math.log(1 + a * numRecs)
    def color_str(recInd, strand) :
        binInd = (int)(numBins * scale_rank(recInd))
        # bin 0 -> colorStrs.
        # bin numBins -> white = (255,255,255) (never happens cause bin < numBins)
        fullRGB = map(int, colorStrs[strand].split(','))
        whiteRGB = (255, 255, 255)
        return ','.join('%d' % (fullRGB[ii] * (1 - binInd / numBins) +
                                whiteRGB[ii] * binInd / numBins)
                        for ii in range(3))

    with myopen(novelFileName, 'w') as novelFile :
        for recInd, rec in enumerate(recs) :
            chrom = rec.Chrom
            bedLine = intervals_to_bed_line(chrom, [(rec.Start, rec.End)], rec.Strand,
                                            recInd + 1, rec.Start, rec.End,
                                            color = color_str(recInd, rec.Strand))
            print(bedLine, file = novelFile)

# *********** Utilities ***************************

def chromInt2Str(chrom, intervals) :
    return '+'.join(chrom + ':%d-%d' % tuple(interval) for interval in intervals)

def same_frame_CDSs(interval, strandOfInterval, tr, includeStop) :
    """
    Return a list of intervals on the chromosome representing the overlap between interval
        and the CDSs of tr if they are on the same strand and in the same frame.
    Treat the stop codon as part of the CDS if specified.
    It is assumed that interval is a multiple of 3, is on a codon boundary, and is on the 
        same chromosome as tr.
    """
    intervalLen = interval[1] - interval[0] + 1
    fakeTr = construct_transcript(tr.get_chromosome(), [interval], strandOfInterval,
                                  codingStart = 1, codingEnd = intervalLen)
    return intersect_tr_intervals(tr, fakeTr,
                   [tr.codingStart, tr.codingEnd - (0 if includeStop else 3)],
                   [1, intervalLen],
                   sameStrandOnly = True,
                   sameOrAntiFrameOnly = True)

def any_same_frame(interval, strandOfInterval, trs) :
    """Does any CDS (without stop) of any of the trs overlap the interval in same frame
       on strand? Assumes all trs are on the same chromosome as interval."""
    return any(len(same_frame_CDSs(interval, strandOfInterval, tr, includeStop = False)) > 0
               for tr in trs)

def any_anti_frame(interval, strandOfInterval, trs) :
    """Does any CDS of any of the trs overlap the interval in the antisense frame of
       strandOfInterval? Assumes all trs are on the same chromosome as interval."""
    return any_same_frame(*anti_interval(interval, strandOfInterval), trs = trs)

def subtract_trs(interval, strandOfInterval, pseudoOverlapTrs, codingOverlapTrs) :
    """
    Return new intervals in which overlaps with trs have been subtracted.
    Subtract overlaps with pseudoOverlapTrs exons in any frame on either strand.
    Subtract overlaps with codingOverlapTrs CDSs (with stop) only if they are in the same or 
        antisense frame.
    Trim to codon boundaries of the input intervals.
    It is assumed that interval is a multiple of 3.
    Assumes all trs and all intervals are on the same chromosome.
    """
    assert (interval[1] - interval[0] + 1) % 3 == 0, interval
    toSubtract = set([exon for tr in pseudoOverlapTrs for exon in tr.get_exon_intervals()])
    antiInterval, antiStrand = anti_interval(interval, strandOfInterval)
    for tr in codingOverlapTrs :
        if tr.get_strand() == strandOfInterval :
            toSubtract.update(same_frame_CDSs(interval, strandOfInterval, tr, includeStop = True))
        else :
            toSubtract.update(same_frame_CDSs(antiInterval, antiStrand, tr, includeStop = True))
    # Extend to codon boundaries of interval
    toSubtract = list([start, end] for start, end in toSubtract)
    for ind, (start, end) in enumerate(toSubtract) :
        start -= (start - interval[0]) % 3
        end += (interval[1] - end) % 3
        toSubtract[ind] = [start, end]
                              
    newIntervals = [list(interval)]
    for intToSubtract in toSubtract :
        newIntervals = subtract_regions(newIntervals, [intToSubtract])

    assert all((s - interval[0]) % 3 == 0 and (interval[1] - e) % 3 == 0 for s, e in newIntervals)

    return newIntervals

if __name__ == '__main__' :
    main()
