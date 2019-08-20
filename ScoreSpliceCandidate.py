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
ScoreSpliceCandidate.py

Given the context of a GT or AG in a string of DNA bases, and precomputed coefficients
files, return a score that indicates how likely it is to be a splice donor or acceptor,
respectively, using the maximum entropy method [1].

The donor coefficient file consists of 4^7 numbers stored in little-endian double format.
The index into this file for a particular splice candidate is computed from the string
of 3 bases 5' of the GT concatenated to the 4 bases 3' of the GT. The index is a base-4
number created from the indices of the bases of this base string in ACGT with the first
being most significant. For example, the index for a potential splice site context of:
ACG-GT-TACG is 0*4^6 + 1*4^5 + 2*4^4 + 3*4^3 + 0*4^2 + 1*4^1 + 2*4^0
Once the coefficient is extracted, the score is
	log(16.302010666666664 * coeff, 2)

The acceptor coefficient file consists of 9 concatenated sequences containing
4^7, 4^7, 4^7, 4^7, 4^7, 4^3, 4^4, 4^3, and 4^4 numbers. The relevant context
is 18 bases 5' of the AG concatenated to 3 bases 3'. The 9 coefficients are extracted
from these arrays using indices computed from the substrings starting at the following
offsets with the following lengths:
0 7, 7 7, 14 7 (including the 3 bases 3' of the AG), 4 7, 11 7, 4 3, 7 4, 11 3, 14 4
Once these 9 coefficients are found, the score is:
log(16.3482025 * (prod 1st 5 coeffs) / (prod other 4 coeffs), 2)

[1] Yeo, G., & Burge, C. B. (2004). Maximum entropy modeling of short sequence motifs with
applications to RNA splicing signals. Journal of Computational Biology : a Journal of
Computational Molecular Cell Biology, 11(2-3), 377-394. doi:10.1089/1066527041410418
"""
from __future__ import division
from __future__ import print_function
import os, struct, math

class DonorPredictor(object) :
    def __init__(self, donorCoefFileName) :
        self.file = open(os.path.abspath(os.path.expanduser(donorCoefFileName)))
        self.file.seek(0, os.SEEK_END)
        assert self.file.tell() / RecordSize == 16384, \
            '%s is not a valid donor coefficients file.' % donorCoefFileName

    def __call__(self, prev3bases, next4bases) :
        # Score potential GT splice donor site given 3 prev bases and 4 next ones.
        index = _bases_to_number(prev3bases + next4bases)
        self.file.seek(index * RecordSize)
        coeff = struct.unpack(RecordFormat, self.file.read(RecordSize))[0]
        return math.log(16.302010666666664 * coeff, 2)

class AcceptorPredictor(object) :
    def __init__(self, acceptorCoefFileName) :
        self.file = open(os.path.abspath(os.path.expanduser(acceptorCoefFileName)))
        self.file.seek(0, os.SEEK_END)
        assert self.file.tell() / RecordSize == 82560, \
            '%s is not a valid acceptor coefficients file.' % acceptorCoefFileName

    def __call__(self, prev18bases, next3bases) :
        # Score potential AG splice acceptor site given 18 prev bases and 3 next ones.
        bases = prev18bases + next3bases
        coeffsCombination = 1
        for ii, (start, end) in enumerate(AcceptorStartEnds) :
            index = _bases_to_number(bases[start : end + 1])
            self.file.seek((AcceptorArrayLengthSums[ii] + index) * RecordSize)
            coeff = struct.unpack(RecordFormat, self.file.read(RecordSize))[0]
            if ii < 5 :
                coeffsCombination *= coeff
            else :
                coeffsCombination /= coeff
        return math.log(16.3482025 * coeffsCombination, 2)

def _bases_to_number(bases) :
    "Convert a string of DNA bases to a base-4 number."
    BaseMap = {'A' : 0, 'C' : 1, 'G' : 2, 'T' : 3}
    result = 0
    for b in bases :
        result = 4 * result + BaseMap[b]
    return result

RecordFormat = '<d' # little-endian double
RecordSize = struct.calcsize(RecordFormat) # Typically 8

# Relevant intervals for acceptor site prediction
AcceptorStartEnds = [(0, 6), (7, 13), (14, 20), (4, 10), (11, 17),
                     (4, 6), (7, 10), (11, 13), (14, 17)]
# Lengths of the coefficient arrays for acceptor sites.
AcceptorArrayLengths = [4 ** (end - start + 1) for start, end in AcceptorStartEnds]
AcceptorArrayLengthSums = [sum(AcceptorArrayLengths[:ii])
                           for ii in range(len(AcceptorArrayLengths))]
