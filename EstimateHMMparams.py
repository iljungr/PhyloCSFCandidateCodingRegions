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
EstimateHMMparams.py
Estimate the best parameters for the HMM used to calculate PhyloCSF Regions.
"""
from __future__ import division
from __future__ import print_function
import sys, random
from math import log10, log
from CommonUtils import myopen, neighbors

assert sys.version_info[:2] == (2, 7), 'This is intended to run in Python 2.7.'

def estimate_hmm_params_for_genome(codingExonsFileName, genomeLength) :
    """
    codingExonsFileName is the name of a file containing information about every annotated
        coding exon in the genome. Each line contains five Tab-separated fields, namely
            chromosome, strand, frame, start, end (start <= end)
        for the annotated coding portion of every exon in the genome.
        Frame is the remainder mod 3 of the chromosomal coordinate of the first base of a
        codon, first being counted along the + strand even if the exon is on the - strand.
    genomeLength is the sum of the number of nucleotides of all chromosomes (and
        scaffolds) in the genome assembly.
        
    Return a 4-tuple containing the parameters needed by get_coding_hmm, namely:
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
    """
    exonList = []
    for line in myopen(codingExonsFileName) :
        line = line.strip()
        words = line.split('\t')
        exonList.append((words[0], words[1], int(words[2]), int(words[3]), int(words[4])))

    exonsByChrStrFr = {} # {(chrom, strand, frame) : [(start1, end1), (start2, end2),...]
    exonList.sort()
    for (chrom, strand, frame, start, end) in exonList :
        exonsByChrStrFr.setdefault((chrom, strand, frame), []).append((start, end))
    gapsNT = [] # Gaps between consecutive non-overlapping coding exons in same frame
                # Note: ideally this should be adjusting e1 and s2 to codon boundaries,
                # since we are interested in the gap between in-frame codons (even if we
                # are reporting it in NT).
    numExons = 0
    totalCodingLengthNT = 0
    for pairsList in exonsByChrStrFr.values() : # Deleting from sublist is much faster
        index = 0
        while index < len(pairsList) - 1 :
            start1, end1 = pairsList[index]
            start2, end2 = pairsList[index + 1]
            if start2 <= end1 :
                if end1 - start1 >= end2 - start2 :
                    del pairsList[index + 1]
                else :
                    del pairsList[index]
            else :
                index += 1
        gapsNT.extend(start2 - end1 - 1
                      for (start1, end1), (start2, end2) in neighbors(pairsList)
                      if start2 > end1 + 1) # Exclude 0 length
        numExons += len(pairsList)
        totalCodingLengthNT += sum(end - start + 1 for start, end in pairsList)

    # Estimate distribution of gaps between coding regions as a mixture of exponential distributions
    nonCodingLengthsNT, nonCodingWeights = estimate_gap_mixture_model(gapsNT, 3, numSteps = 20)
    nonCodingLengthsInCodons = [x / 3 for x in nonCodingLengthsNT]

    codingPrior = totalCodingLengthNT / genomeLength / 6 # Prior for being coding in a particular frame
    codingLengthInCodons = totalCodingLengthNT / numExons / 3 # Mean length in codons
    
    return codingPrior, codingLengthInCodons, nonCodingWeights, nonCodingLengthsInCodons

def estimate_gap_mixture_model(gapsNT, numInMix = 1, numSteps = 10, relxtol = 0.001) :
    # Return the mean lengths (in NT) and priors for a mixture model of the distribution of gaps
    #    between same-frame exons.
    assert(1 <= numInMix <= 3)
    
    from MixtureModelUtil import infer_mixture, ParametrizedLogDistribution
    
    MAX_NUM_GAPS = 20000  # 20000 is more than enough; no need to waste time with more
    if len(gapsNT) > MAX_NUM_GAPS :
        random.seed(0)
        gapsNT = [gapsNT[random.randint(0, len(gapsNT) - 1)]
                  for ii in range(MAX_NUM_GAPS)]
    
    # Starting guesses found by eyeballing mouse gap distribution.
    guessLengths = [3e3, 8e4, 1e2][:numInMix]
    guessPriors = [30, 10, 1][:numInMix]
    guessPriors = [x / sum(guessPriors) for x in guessPriors]
    
    class ExpDist(ParametrizedLogDistribution) :
        def log_density_and_grad(self, point, params) :
            # Exponential distribution with logarithmic parameter.
            # Density(x) = 10 ** params[0]
            tau = 10 ** params[0]
            return - point / tau - log(tau)

    expDist = ExpDist()
    paramGuesses = [[log10(x)] for x in guessLengths]
    
    paramSeqs, classPriors, classProbs = infer_mixture(gapsNT, [expDist] * numInMix, paramGuesses,
                                                                    guessPriors,
                                                                    numSteps = numSteps, relxtol = relxtol,
                                                                    dump = False)
    gapMeanLengthsNT = [10 ** paramSeq[0] for paramSeq in paramSeqs]
    return gapMeanLengthsNT, classPriors
