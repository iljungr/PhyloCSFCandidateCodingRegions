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
from __future__ import division
import sys
from itertools import islice, imap
from math import log10, floor
from CommonUtils import pjoin, myopen
from ClassHMM import HMM

assert sys.version_info[:2] == (2, 7), 'This is intended to run in Python 2.7.'

MaxLogOdds = 15.0
def prob_to_log_odds(prob) :
    if prob < 10 ** -MaxLogOdds :
        return -MaxLogOdds
    if prob > 1 - 10 ** -MaxLogOdds :
        return MaxLogOdds
    return log10(prob / (1 - prob))

def create_PhyloCSF_Regions(
    hmmParams,         # Parameters for get_coding_hmm
    phyloCSFoutputDir, # Directory containing input file with phyloCSF scores
    phyloCSFregionDir, # Directory to contain output bed file
    chrom, strand, frame) :
    """
    Given a file with PhyloCSF scores (including bls) of consecutive codons for one frame
        of one strand of a chromosome, create corresponding .bed file of the PhyloCSF
        Regions in most likely path.
    The PhyloCSF scores file must be named {chrom}.Strand{strand}.Frame{frame}.fixed.out,
        e.g., chr1.Strand+.Frame2.fixed.out. Frame is 0, 1, or 2.
    Each line of the PhyloCSF scores file must be of the form:
        CHROM:CODON_START-CODON_END   STRAND  score(decibans) PHYLOCSF_SCORE  BLS
    where PHYLOCSF_SCORE and BLS are the output of PhyloCSF using the --strategy=fixed
    and --bls options. Codons should be consecutive, in increasing nominal order (i.e.,
    the numbers are increasing even for scores on the minus strand). The separator between
    fields must be a Tab character. For codons for which there is no alignment, put
    a single field containing "No_Alignment" instead of the last 3 fields.
    For the resulting bed file to be usable in the UCSC browser, chromosome names should
    be UCSC chromosome names, and positions must be 1-based.
    For example:
        chr1:10915-10917        -       No_Alignment
        chr1:10918-10920        -       score(decibans) 0.6017  0.0023
    The output bed file will be named {chrom}.Strand{strand}.Frame{frame}.coding.bed.
    """
    phyloCSFfileName = pjoin(phyloCSFoutputDir,
                             '%s.Strand%s.Frame%s.fixed.out'  % (chrom, strand, frame))
    outputBedFileName = pjoin(phyloCSFregionDir,
                             '%s.Strand%s.Frame%s.coding.bed' % (chrom, strand, frame))
                             
    print >>sys.stderr, 'Processing %s.' % phyloCSFfileName
    minRelBranchLength = 0.1
    bedFile = myopen(outputBedFileName, 'w')

    def processScores(scores, blockStartPos, chrom, strand) :
         # Process and write one block of scores.
        if len(scores) == 0 :
            return
        
        hmm = get_coding_hmm(*hmmParams)
        codingProbabilities = hmm.state_probabilities(scores)[:, 0]
        bestPath = [state == 0 for state in hmm.best_path(scores)]

        curCodonCount = 0
        for chunk in _chunkify(bestPath, lambda elt1, elt2 : elt1 == elt2) :
            # Iterate through chunks of positions in the same state
            if chunk[0] : # Only write coding chunks
                maxCodingProb = max(islice(codingProbabilities, curCodonCount,
                                           curCodonCount + len(chunk)))
                maxLogOdds = prob_to_log_odds(maxCodingProb)
                # 8 possible gray scales, 0, 30, ..., 210 (more than 210 is too light)
                # Divide into 8 bins based on maxLogOdds from 0-8 (but handle extremes)
                # (In sample region, 90th percentile of positive scores was 7.)
                grayScale = 210 if maxLogOdds < 1 else \
                            0 if maxLogOdds > 7 else   \
                            210 - 30 * int(floor(maxLogOdds))
                chunkStartPos = blockStartPos + 3 * curCodonCount
                chunkEndPos = chunkStartPos + 3 * len(chunk) - 1
                print >>bedFile, '\t'.join(map(str,[
                    chrom,
                    chunkStartPos - 1, # Bed counts from 0 instead of 1.
                    chunkEndPos, # Count from 0, but chromEnd is first position _after_ end
                    '%s:%d-%d' % (chrom, chunkStartPos, chunkEndPos), # Name
                    0,
                    strand,
                    chunkStartPos - 1, # thickStart
                    chunkEndPos,       # thickEnd
                    '%d,%d,%d' % (grayScale, grayScale, grayScale), # itemRgb
                    ]))
            curCodonCount += len(chunk)
    
    scores = []
    strand = None
    frame = None
    for lineCount, line in enumerate(open(phyloCSFfileName)) :
        words = line.split()
        if words[2] == 'No_Alignment' :
            continue
        region = words[0]
        score = float(words[3])
        relBranchLength = float(words[4])
        chrom, interval = region.split(':')
        pos = int(interval.split('-')[0])
        if strand == None :
            strand = words[1]
        else :
            assert words[1] == strand, 'Strand mismatch.'
        if frame == None :
            frame = pos % 3
        else :
            assert pos % 3 == frame, 'Position %d does not match frame %d.' % (pos, frame)
        
        if len(scores) > 0 and (
                relBranchLength <= minRelBranchLength or
                chrom != prevChrom or
                pos != prevPos + 3) :
            processScores(scores, blockStartPos, prevChrom, strand)
            del scores[:] # Reset to 0 but allow Python to reuse list (don't know if it does)
            
        if relBranchLength > minRelBranchLength :
            if len(scores) == 0 :
                blockStartPos = pos
            scores.append(score)
            prevChrom = chrom
            prevPos = pos
    else :
        if len(scores) > 0 :
            processScores(scores, blockStartPos, chrom, strand)
    bedFile.close()

def get_coding_hmm(codingPrior, codingNumCodons, nonCodingWeights, nonCodingNumCodonsList) :
    """
    Return an hmm to compute the probability that a particular codon is coding given the PhyloCSF scores
        of a set of consecutive codons. Ignores sequence properties like stop codons, splice sites, etc
        and just treats the scores as independent log likelihoods for individual codons and relies
        on the fact that coding codons come in consecutive clusters.
    There is one coding state but there can be several noncoding states. The noncoding states all have the
        same emissions probabilities, but they can have different probabilities of transitioning to coding
        to represent noncoding regions of different lengths (introns, intergenic, etc); they never
        transition to other non-coding regions.
        
    Input:
        codingPrior: probability that a random codon is coding
        codingNumCodons: typical length of a coding region, in codons; the reciprocal of
            the probability of going from a coding state to a noncoding state.
        nonCodingWeights: a list representing the probability that a coding to noncoding
            transition will go into a particular noncoding state (sum is 1).
            len(nonCodingWeights) determines the number of noncoding states.
        nonCodingNumCodonsList: typical region length for each of the noncoding states, in
            codons; the reciprocal of the probability of going from that state to the
            coding state. Must have same length as nonCodingWeights.
            Priors for each noncoding state will be computed from the weights and lengths.
    Output:
        HMM with state 0 representing coding
    """
    numNonCoding = len(nonCodingWeights)
    assert numNonCoding == len(nonCodingNumCodonsList)
    assert abs(1 - sum(nonCodingWeights)) < 1e-9
    unnormalizedNoncodingPriors = [weight * length
                                   for weight, length in zip(nonCodingWeights,
                                                             nonCodingNumCodonsList)]
    initProbs = [codingPrior] + [(1 - codingPrior) *
                                 unnormalizedNoncodingPrior / sum(unnormalizedNoncodingPriors)
                                 for unnormalizedNoncodingPrior in unnormalizedNoncodingPriors]
    cToNCprobs = [nonCodingWeight / codingNumCodons for nonCodingWeight in nonCodingWeights]
    ncToCprobs = [1 / length for length in nonCodingNumCodonsList]
    transProbs = [[1 - sum(cToNCprobs)] + cToNCprobs] + \
                 [[ncToCprob] + [1 - ncToCprob if ii == jj else 0
                                 for jj in range(numNonCoding)]
                  for ii, ncToCprob in enumerate(ncToCprobs)]
    
    return HMM(numStates = 1 + numNonCoding,
               initProbs = initProbs,
               transProbs = transProbs,
               emitProb = (lambda state, score : 10 ** (score / 10) if state == 0 else 1)
                    # Recall that 10**(score/10) = P(alignment|coding) / P(alignment|noncoding).
                    # If we think of the alignment as the emission, then the emission
                    # probabilities should be:
                    # P(alignment | noncoding) and P(alignment | coding).
                    # But scaling all of them doesn't affect best path or state
                    # probabilities, so scale by the noncoding probability.
              )

def _chunkify(iterable, continueChunk) :
    """
    Break up the iterable into chunks, given a function that species whether to continue
        a chunk from one element to the next.
    continueChunk = lambda prevElt, thisElt : True or False
    Yield lists of consecutive elements.
    """
    curChunk = []
    for elt in iterable :
        if curChunk == [] or continueChunk(curChunk[-1], elt) :
            curChunk.append(elt)
        else :
            yield curChunk
            curChunk = [elt]
    if curChunk != [] :
        yield curChunk

