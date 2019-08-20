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
Examples.py
Run examples that use the key algorithms of the PhyloCSFCandidateCodingRegions project.
"""
from __future__ import division
from __future__ import print_function
import os, sys
from itertools import izip_longest
from CommonUtils import pushd, popd, equal_files, err_msg, cp, rm, pjoin
from HMMforPhyloCSFRegions import create_PhyloCSF_Regions
from HMMparamsForMudge2019 import HMMparams
from EstimateHMMparams import estimate_hmm_params_for_genome
from SVMsForPhyloCSFCandidateCodingRegions import classify_regions, do_other_steps, get_reader
from ScoreSpliceCandidate import DonorPredictor, AcceptorPredictor

assert sys.version_info[:2] == (2, 7), 'This is intended to run in Python 2.7.'

def main() :
    if len(sys.argv) != 2 or sys.argv[1] not in ['svm', 'hmm', 'splice'] :
        print('Run with argument "svm", "hmm", "splice".', file = sys.stderr)
    elif sys.argv[1] == 'hmm' :
        hmm_example()
    elif sys.argv[1] == 'svm' :
        svm_example()
    else :
        splice_example()

def hmm_example() :
    "Calculate HMM parameters and then run HMM for one reading frame of one scaffold"
    
    pushd(os.path.dirname(__file__))
    
    # This file has information about every annotated coding exon in human GENCODE v29
    codingExonsFileName = 'ExampleFiles/EstimateHMMparamsExample/HumanCodingExonsV29.txt'
    
    humanGenomeLength = 3252208893 # Sum of chromosome/scaffold lengths in hg38 assembly
    err_msg('Estimating HMM parameters.')
    humanHMMparams = estimate_hmm_params_for_genome(codingExonsFileName, humanGenomeLength)
    
    print('Here are the HMM parameters computed for the Human Genome from GENCODE v29.')
    print(humanHMMparams)
    print()
    """
    (0.0018352187276509814, 56.07319200374541,
        [0.6533594654727364, 0.20618155683835557, 0.140458977688905],
        [3299.5676063682404, 118314.79205686758, 314.61620563310396])
    """
    
    print('For comparison, here are the ones used for the paper (based on v16).')
    print(HMMparams['Human'])
    print()
    """
    (0.0018577729491349902, 56.83017481932195,
        [0.6417673350475721, 0.2150725662530122, 0.1431600986994177],
        [3376.1207312280044, 122153.58921094926, 328.20142626163494])
    """

    # The result of running PhyloCSF on every codon in frame 1 on the minus
    # strand of chr4_GL000008v2_random is in
    # ExampleFiles/HMMexample/chr4_GL000008v2_random.Strand-.Frame1.fixed.out
    
    phyloCSFoutputDir = phyloCSFregionDir = 'ExampleFiles/HMMexample'
    
    err_msg('Computing PhyloCSF Regions.')
    create_PhyloCSF_Regions(humanHMMparams, phyloCSFoutputDir, phyloCSFregionDir,
                            chrom = 'chr4_GL000008v2_random', strand = '-', frame = 1)
    
    phyloCSFregionFileName = 'ExampleFiles/HMMexample/chr4_GL000008v2_random.Strand-.Frame1.coding.bed'
    assert equal_files(phyloCSFregionFileName,
                       'ExampleFiles/Results/chr4_GL000008v2_random.Strand-.Frame1.coding.bed')
    err_msg('\nYay! Result file matches precomputed file.')
    rm(phyloCSFregionFileName)
    popd()
          
def svm_example() :
    pushd(os.path.dirname(__file__))
    homeDir = 'ExampleFiles/SVMexample'
    regionsDir = homeDir # Where previously computed PhyloCSF Regions bed files are
    
    # The following bed files have all GENCODE v29 coding and pseudogene transcripts on chr17
    codingBedFileName = 'ExampleFiles/SVMexample/CodingTranscripts.GENCODEv29.chr17.bed'
    pseudoBedFileName = 'ExampleFiles/SVMexample/PseudogeneTranscripts.GENCODEv29.chr17.bed'
    
    # Classify regions based on overlap with annotations and write input file for PhyloCSF
    classify_regions(homeDir, regionsDir, codingBedFileName, pseudoBedFileName)

    # Simulate running PhyloCSF with strategy=mle by copying the resulting file
    err_msg('Simulating running PhyloCSF with strategy=mle by copying the resulting file')
    cp('ExampleFiles/SVMexample/SimulatedPhyloCSFOutput/Regions.pcsf.out',
       'ExampleFiles/SVMexample')

    # Fill in PhyloCSF fields, train and run SVMs, prune and sort regions, make bed file
    do_other_steps(homeDir, 100) # Real SVM used 10000 training vectors, but this is faster

    for fileName in ['Regions.01.txt', 'Regions.02.txt', 'Regions.04.txt',
                     'Regions.pcsf.in', 'Regions.pcsf.out', 'PhyloCSFNovel.bed'] :
        assert equal_files(pjoin('ExampleFiles/SVMexample', fileName),
                           pjoin('ExampleFiles/Results',    fileName)), 'Files differ: %s' % fileName
        rm(pjoin('ExampleFiles/SVMexample', fileName))
    # Some values in Regions.03.txt differ slightly on linux and mac, so can't use equal_files.
    for lineNew, lineOld in izip_longest(get_reader('ExampleFiles/SVMexample/Regions.03.txt'),
                                         get_reader('ExampleFiles/Results/Regions.03.txt')) :
        for key in lineNew.keys() :
            if isinstance(lineNew[key], float) :
                assert abs(lineNew[key] - lineOld[key]) < 1e-9
            else :
                assert lineNew[key] == lineOld[key]
    rm('ExampleFiles/SVMexample/Regions.03.txt')
                     
    err_msg('\nYay! Result files match precomputed files.')
    popd()

def splice_example() :
    donorPredictor    = DonorPredictor(   'ExampleFiles/SpliceExample/Hsap.donor.mecoef')
    acceptorPredictor = AcceptorPredictor('ExampleFiles/SpliceExample/Hsap.acceptor.mecoef')

    # Predict donor score for TCA-GT-AAGG
    bases = 'TCAGTAAGG'
    score = donorPredictor(bases[0 : 3], bases[5 : 9])
    assert abs(score - 5.685921922927156) < 1e-9, score

    # Predict acceptor score for CAATGGTTAGTTTCAGTA-AG-GAA
    bases = 'CAATGGTTAGTTTCAGTAAGGAA'
    score = acceptorPredictor(bases[0 : 18], bases[20 : 23])
    assert abs(score - -3.033034001453567) < 1e-9, score

    err_msg('\nYay! Results match precomputed values.')

if __name__ == '__main__' :
    main()
