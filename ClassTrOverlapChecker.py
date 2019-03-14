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
Cache to rapidly check if a set of intervals overlaps any CDS or exon in a bunch of transcripts.
"""
from __future__ import division
from __future__ import print_function
import sys

assert sys.version_info[:2] == (2, 7), 'This is intended to run in Python 2.7.'

from ClassPointsInBoxCache import PointsInBoxCache

class TrOverlapChecker :
    """
    Cache to rapidly check if a set of intervals overlaps any CDS or exon in a bunch of transcripts.
    """
    def __init__(self, transcripts, onlyCDS = True, includeStop = True) :
        """
        If onlyCDS, only looks for overlap with CDS, including stop if specified, 
        otherwise looks for overlap with exons.
        """
        self.transcripts = transcripts
        self.caches = {}  # (chrom, strand) : PointsInBoxCache for that chrom and strand
        self.onlyCDS = onlyCDS
        self.includeStop = includeStop # Ignored unless onlyCDS
        
    def overlapping_trs(self, chrom, strands, interval) :
        """Return the list of all transcripts whose exons overlap the interval on the specified 
               strand or strands. (Strand can be '+', '-', or '+-')
           If self.onlyCDS, only look for overlaps in the CDS portion of the exons.
        """
        resultSet = set()
        for strand in strands :
            cache = self._get_pointsInBoxCache(chrom, strand)
            for begin, end, tr in cache.points_in_box(0, interval[1], interval[0], sys.maxint) :
                resultSet.add(tr)
        return sorted(list(resultSet),
                      key = lambda tr : tr.get_name()) # Sort for reproducibility

    def overlapping_trs_intervals(self, chrom, strands, intervals) :
        """Return the list of all transcripts whose exons overlap at least one of the intervals
               on the specified strand or strands. (Strand can be '+', '-', or '+-')
           If self.onlyCDS, only look for overlaps in the CDS portion of the exons.
        """
        return list(set(sum([self.overlapping_trs(chrom, strands, interval)
                             for interval in intervals], [])))

    def _get_pointsInBoxCache(self, chrom, strand) :
        key = (chrom, strand)
        if key not in self.caches :
            relevantTrs = (tr for tr in self.transcripts
                              if tr.get_chromosome() == chrom and tr.get_strand() == strand)
            points = ((begin, end, tr)
                      for tr in relevantTrs
                      for (begin, end) in (tr.get_CDS_intervals(self.includeStop)
                                           if self.onlyCDS else tr.get_exon_intervals()))
            self.caches[key] = PointsInBoxCache(points)
        return self.caches[key]
