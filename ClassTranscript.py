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
This file contains data structures and utilities for working with transcripts.
"""
from __future__ import division
import sys, os
from CommonUtils import stop, bed_line_to_intervals, neighbors

assert sys.version_info[:2] == (2, 7), 'This is intended to run in Python 2.7.'

##############################################################################
## class Transcript:
##############################################################################

class Transcript:
    '''   
    trName
    chromosome
    strand
    start
    end
    exonMap # sorted list of triples:
            #     (exonStart, exonEnd, num bases in exons up to and including this)
    length  # cache: total number of bases in all exons combined

    # Transcript position (TrPos) is a 1-up index into the union of exons in
    #   the strand direction.
    # TrPos < 1 or > length is allowed and is computed as if
    #    the transcript were extended beyond the 1st or last exon
    #    assuming no further splicing.
    codingStart (TrPos of 1st base of 1st CDS)
    codingEnd (TrPos of last base of last CDS)
    '''
    # Methods to extract transcript data
    def get_name(self) :
        "Return transcript ID."
        return self.trName
    def get_chromosome(self) :
        return self.chromosome
    def get_strand(self) :
        return self.strand
    def get_start(self) :
        return self.start
    def get_end(self) :
        return self.end
    def get_length(self) :
        return self.length
    def get_coding_start(self) :
        return self.codingStart
    def get_coding_end(self) :
        return self.codingEnd
    def __str__(self) :
        return 'Transcript(%s)' % self.get_name()

    # Methods to convert between transcript position and chromosomal position
    def trPos_to_pos(self, trPos) :
        if trPos == None :
            return None # Allows pos == trPos_to_pos(pos_to_trPos(pos)) to detect pos in intron
        if self.strand == '-' :
            trPos = self.length - trPos + 1
        for exonStart, exonEnd, cumNumBases in self.exonMap :
            if trPos <= cumNumBases :
                return exonEnd + trPos - cumNumBases
        exonStart, exonEnd, cumNumBases = self.exonMap[-1]
        return exonEnd + trPos - cumNumBases
    def pos_to_trPos(self, pos) :
        """
        Given a chromosome position, return the corresponding transcript position.
        Treat the transcript as extending infinitely without splicing past the start and end.
        If pos is within an intron return None.
        """
        trPos = 0
        for exonStart, exonEnd, cumNumBases in self.exonMap :
            if exonEnd < pos :
                trPos += exonEnd - exonStart + 1
            else :
                if pos < exonStart and trPos > 0 :
                    # Before the start of a non-first exon -> in an intron
                    return None
                trPos += pos - exonStart + 1 # Handles case of pos < transcript start
                break
        else : # Past end of last exon
            trPos += pos - self.exonMap[-1][1] # Add distance from last exon end
        if self.strand == '-' :
            trPos = self.length - trPos + 1
        return trPos
    def trInterval_to_intervals(self, start, end) :
        # Input: start and end trPos's with start <= end + 1 (result is empty if start == end + 1)
        # Output: (chromosome, intervals, strand)
        #   intervals is a tuple of pairs, (start, end) on the chromosome
        #   with start <= end, in increasing nomimal order
        # Implementation Note: variable beginning with TrU represent undirected trPos
        assert start <= end + 1, ('trInterval_to_intervals: start > end %d %d %s' % (start, end,
                                                                                     self.trName))
        if start > end :
            return self.chromosome, (), self.strand
        strand = self.strand
        intervals = [] # Will convert to tuple at the end
        if strand == '-' :
            trUStart, trUEnd = self.length - end + 1, self.length - start + 1
        else :
            trUStart, trUEnd = start, end
        prevExonTrUEnd = 0
        for index, (exonStart, exonEnd, cumNumBases) in enumerate(self.exonMap) :
            first = index == 0 
            last = index == len(self.exonMap) - 1
            exonTrUStart = prevExonTrUEnd + 1
            prevExonTrUEnd = exonTrUEnd = cumNumBases
            if exonTrUEnd < trUStart and not last :
                continue
            if trUEnd < exonTrUStart and not first :
                break
            if first : # Treat 1st exon as extending infinitely to the left
                segTrUStart = trUStart
            else :
                segTrUStart = max(trUStart, exonTrUStart)
            if last : # Treat last exon as extending infinitely to the right
                segTrUEnd = trUEnd
            else :
                segTrUEnd   = min(trUEnd,   exonTrUEnd  )
            offset = exonStart - exonTrUStart
            intervals.append((offset + segTrUStart, offset + segTrUEnd))
        return self.chromosome, tuple(intervals), strand
    
    ## Methods to get feature intervals (positions on chromosome, not trPos):
    def get_exon_intervals(self) :
        # Return [(start, end), ...] for exons in the transcript, sorted along positive strand.
        return [(start, end) for start, end, cumNumBases in self.exonMap]
    def get_intron_intervals(self) :
        # Return [(start, end), ...] for introns in the transcript, sorted along positive strand.
        return [(pEnd + 1, nStart - 1) for ((pStart, pEnd, pDummy), (nStart, nEnd, nDummy)) in 
                                            zip(self.exonMap[: -1], self.exonMap[1:])]
    def get_CDS_intervals(self, includeStop = True, trimToCodons = False) :
        """
        Return [(start, end), ...] for CDSs in the transcript, sorted along + strand.
        If includeStop = False, exclude the last 3 bases.
        If codingStart is None, return None.
        If trimToCodons, trim to codon boundaries.
        """
        start = self.codingStart
        if start == None :
            return None
        end = self.codingEnd if includeStop else self.codingEnd - 3
        if trimToCodons :
            start, end = self.trim_trInt_to_codons(start, end)
            if start == None : # 0 length after trimming
                return ()
        if start > end :
            assert trimToCodons or not includeStop, (start, end, self)
            return []
        return self.trInterval_to_intervals(start, end)[1]

    # Miscellaneous methods
    def frame_of_trPos(self, trPos) :
        """
        Return 0, 1, or 2 depending on whether trPos is the 1st, 2nd, or 3rd
          base of a codon in the reading frame of self, respectively
          (assume translation is not restricted by the start and stop codons).
        """
        return (trPos - self.codingStart) % 3
    def frame_of_base(self, pos) :
        """
        Return frame of a position on self.chromosome, as defined in frame_of_trPos.
        Return None if pos is in an intron or codingStart is None.
        If pos is before the beginning or or after the end of the transcript,
          consider the transcript to be extended indefinitely beyond its ends without
          introns, as in trPos.
        """
        trPos = self.pos_to_trPos(pos)
        if trPos == None : # In an intron
            return None
        return self.frame_of_trPos(trPos)
    def num_codons(self, includeStop = True) :
        """
        Return number of codons in the CDS, or None if codingStart is undefined.
        If not includeStop, exclude final 3 bases.
        """
        numCodons = (self.codingEnd - self.codingStart + 1) // 3
        if includeStop :
            return numCodons
        else :
            return numCodons - 1
    def trim_trInt_to_codons(self, startTrPos, endTrPos, offset = 0) :
        return trim_trInt_to_codons(self, startTrPos, endTrPos, offset)

##############################################################################
####### Utilities for making transcripts ########
##############################################################################

def construct_transcript(chromosome, exonIntervals, strand, trName = '',
                         codingStart = None, codingEnd = None,
                         codingStartChromPos = None, codingEndChromPos = None) :
    """
    exonIntervals should be increasing along the + strand.
    """
    tr = Transcript()
    tr.trName = trName
    tr.chromosome = chromosome
    tr.strand = strand
    tr.start = min(s for s, e in exonIntervals)
    tr.end   = max(e for s, e in exonIntervals)
    tr.exonMap = []
    totalBasesSoFar = 0
    for s, e in exonIntervals :
        totalBasesSoFar += e - s + 1
        tr.exonMap.append((s, e, totalBasesSoFar))
    tr.length = sum(e - s + 1 for s, e in exonIntervals)
    if codingStartChromPos != None :
        assert codingStart == None, (codingStart, codingStartChromPos, tr)
        codingStart = tr.pos_to_trPos(codingStartChromPos)
    if codingEndChromPos   != None :
        assert codingEnd   == None, (codingEnd,   codingEndChromPos,   tr)
        codingEnd   = tr.pos_to_trPos(codingEndChromPos)
    tr.codingStart = codingStart
    tr.codingEnd = codingEnd
    return tr

def bed_line_to_tr(line) :
    """
    Create a transcript from a line in bed format.
    If thickStart and thickEnd are included and unequal set codingStart and codingEnd.
    """
    name, chrom, intervals, strand = bed_line_to_intervals(line)

    # Set codingStart and codingEnd from thickStart and thickEnd
    codingStartChromPos = codingEndChromPos = None
    words = line.split('\t')
    if len(words) > 6 :
        thickStart = int(words[6])
        thickEnd   = int(words[7])
        if thickStart != thickEnd :
            thickStart += 1 # Convert from half-open 0-based to closed 1-based positions
            codingStartChromPos = thickEnd   if strand == '-' else thickStart
            codingEndChromPos   = thickStart if strand == '-' else thickEnd
            if strand == '-' :
                assert codingStartChromPos >= codingEndChromPos, line
            else :
                assert codingStartChromPos <= codingEndChromPos, line

    tr = construct_transcript(chrom, intervals, strand, name,
                              codingStartChromPos = codingStartChromPos,
                              codingEndChromPos = codingEndChromPos)
    return tr

##############################################################################
####### Utilities for set operations on sets of intervals ########
##############################################################################

def _keep_region_types(result, *inputs) :
    # Cast result region to have the sequence types of the inputs (first non-empty)
    # For example if inputs are lists of tuples, make result be a list of tuples.
    # Note that both inner and outer level objects get copied.
    outerType = type(inputs[0])
    innerType = type(())
    for input in inputs :
        if len(input) > 0 :
            innerType = type(input[0])
            break
    return outerType(map(innerType, result))

def intersect_regions(*listOfRegions) :
    """
    Return the region of intersection of the two or more input regions.
    Regions are represented as a sorted sequence of intervals
        [(start1, end1), (start2, end2), ...] with start1 <= end1 < start2, etc.
    When result contains same interval as one of inputs it is a copy.
    """
    def _intersect_regions() :
        if len(listOfRegions) == 0 :
            return []
        elif len(listOfRegions) == 1 :
            return listOfRegions[0]
        elif len(listOfRegions) > 2 :
            return intersect_regions(listOfRegions[0],
                                     intersect_regions(*listOfRegions[1:]))
        # Otherwise exactly 2 regions
        intervals1, intervals2 = listOfRegions
        if len(intervals1) == 0 or len(intervals2) == 0 :
            return []
        if intervals2[0][0] > intervals1[-1][1] or intervals1[0][0] > intervals2[-1][1] :
            return []
        result = []
        intervals1iter = iter(intervals1)
        intervals2iter = iter(intervals2)
        try :
            start1, end1 = intervals1iter.next()
            start2, end2 = intervals2iter.next()
            while True :
                intersectionStart = max(start1, start2)
                intersectionEnd = min(end1, end2)
                if intersectionStart <= intersectionEnd :
                    result.append((intersectionStart,intersectionEnd))
                if end1 < end2 :
                    start1, end1 = intervals1iter.next()
                else :
                    start2, end2 = intervals2iter.next()
        except StopIteration :
            return result
    return _keep_region_types(_intersect_regions(), *listOfRegions)

def complement_region(intervals) :
    """
    Return [-sys.maxint, sys.maxint] - intervals. Regions are represented as
      a sorted sequence of integer intervals [(start1, end1), (start2, end2), ...] with
      start1 <= end1 < start2, etc.
    """
    return _keep_region_types(subtract_regions([(-sys.maxint, sys.maxint)], intervals),
                              intervals)

def union_regions(*listOfRegions) :
    """
    Return the region of the union of the two or more input regions. Regions are
    represented as a sorted sequence of integer intervals
        [(start1, end1), (start2, end2), ...] with start1 <= end1 < start2, etc.
    Assumes that intervals are within [-sys.maxint, sys.maxint].
    When result contains same interval as one of inputs it is a copy.
    """
    if len(listOfRegions) == 0 :
        return []
    result = complement_region(intersect_regions(*map(complement_region, listOfRegions)))
    return _keep_region_types(result, *listOfRegions)

def subtract_regions(intervals1, intervals2) :
    """
    Return intervals1 minus intervals2. Regions are represented as a sorted 
      sequence of integer intervals [(start1, end1), (start2, end2), ...] with 
      start1 <= end1 < start2, etc.
    Inputs can be any sequence type of sequence types; output will match.
    When result contains same interval as one of inputs it is a copy.
    """
    if len(intervals1) == 0 or len(intervals2) == 0 :
        return _keep_region_types(list(intervals1), intervals1) # Make a copy of same type
    intervals2Comp = []
    if intervals1[0][0] < intervals2[0][0] :
        intervals2Comp += [(intervals1[0][0], intervals2[0][0] - 1)]
    intervals2Comp += [(b + 1, c - 1) for ((a, b), (c, d)) in neighbors(intervals2)]
    if intervals2[-1][1] < intervals1[-1][1] :
        intervals2Comp += [(intervals2[-1][1] + 1, intervals1[-1][1])]
    result = intersect_regions(intervals1, intervals2Comp)
    return _keep_region_types(result, intervals1)

##############################################################################
####### Miscellaneous utilities ########
##############################################################################

def trim_trInt_to_codons(tr, startTrPos, endTrPos, offset = 0) :
    """
    Trim the interval (startTrPos, endTrPos) to codon boundaries and 
        return (trimmedStartTrPos, trimmedEndTrPos).
    If offset is non-zero, it must be 1 or 2. In that case, trim to boundaries of "codons"
        starting 1 or 2 bases downstream of the annotated codons, respectively.
    Return None, None if the codingStart is None, or if the resulting interval
        has length <= 0.
    """
    assert offset in (0, 1, 2), offset
    if tr.codingStart == None :
        return None, None
    startTrPos += (3 - tr.frame_of_trPos(startTrPos) + offset) % 3
    endTrPos -= (tr.frame_of_trPos(endTrPos) - 2 - offset) % 3
    return (startTrPos, endTrPos) if startTrPos < endTrPos else (None, None)

def intervals_to_bed_line(chrom, intervals, strand, name = '',
                          thickStartChromPos = None, thickEndChromPos = None,
                          score = 0, color = 0) :
    """
    Return a .bed formatted line with specified data.
    If thickStartChromPos and thickEndChromPos are not none, they should be the minimum and
        maximum 1-up chromosome coordinates of the bases to be shown as thick.
    """
    assert (thickStartChromPos == None) == (thickEndChromPos == None)
    chromStart = intervals[0][0] - 1 # Bed counts from 0 instead of 1.
    chromEnd = intervals[-1][1] # Count from 0, but chromEnd is position _after_ end
    if thickStartChromPos == None :
        thickStart = thickEnd = chromStart
    else :
        thickStart = thickStartChromPos - 1
        thickEnd = thickEndChromPos
    return '\t'.join(map(str, [
        chrom, chromStart, chromEnd, name, score, strand,
        thickStart, thickEnd, color, len(intervals),
        ','.join(map(str, [end - start + 1 for start, end in intervals])),
        ','.join(map(str, [start - chromStart - 1 for start, end in intervals])),
        ]))

def is_same_or_anti_frame(tr1, tr2, chromPos) :
    """
    If tr1 and tr2 are on the same strand, return True if they will read the chromosomal
        position in the same frame. Otherwise, return True if they will read it
        in the antisense frame, i.e., the one in which 3rd codon positions align.
    It is assumed that the transcripts are on the same chromosome.
    Return None if codingStart is None.
    """
    assert tr1.get_chromosome() == tr2.get_chromosome(), (tr1, tr2, chromPos)
    frame1 = tr1.frame_of_base(chromPos)
    frame2 = tr2.frame_of_base(chromPos)
    if frame1 == None or frame2 == None : # codingStart == None
        return None
    if tr1.strand == tr2.strand :
        return frame1 == frame2
    else :
        return frame1 == (1 - frame2) % 3

def intersect_tr_intervals(tr1, tr2, trInterval1, trInterval2,
                           sameStrandOnly = True, sameOrAntiFrameOnly = False) :
    """
    Return the region of overlap between the two transcript intervals in the form:
        [(start1, end1), (start2, end2), ...] with start1 <= end1 < start2, etc.
    trInterval is an interval of transcript positions on the respective transcript.
    If sameStrandOnly, ignore overlaps on opposite strand.
    If sameOrAntiFrameOnly, ignore overlaps that are not in the same frame (if trs are
        on the same strand) or antisense frame (if trs are on opposite strands), and
        return None if codingStart is None for either transcript.
    """
    if tr1.get_chromosome() != tr2.get_chromosome() :
        return []
    if sameStrandOnly and tr1.get_strand() != tr2.get_strand() :
        return []
    intervals1 = tr1.trInterval_to_intervals(*trInterval1)[1]
    intervals2 = tr2.trInterval_to_intervals(*trInterval2)[1]
    if len(intervals1) == 0 or len(intervals2) == 0 :
        return []
    if intervals1[0][0] > intervals2[-1][1] or intervals2[0][0] > intervals1[-1][1] :
        return []                 # Quick reject that shouldn't affect results
    overlapIntervals = intersect_regions(intervals1, intervals2)
    if sameOrAntiFrameOnly :
        if tr1.codingStart == None or tr2.codingStart == None :
            return None
        overlapIntervals = [interval
                            for interval in overlapIntervals
                            if is_same_or_anti_frame(tr1, tr2, interval[0])]
    return overlapIntervals
