The scripts in this git repository contain the main algorithms for computing PhyloCSF Candidate Coding Regions, as defined in "Discovery of protein-coding genes and exons by whole-genome PhyloCSF helps elucidate 118 GWAS loci", by Jonathan M. Mudge, Irwin Jungreis, et al. There is also a script that computes the scores used to create the splice site prediction browser track.

The scripts should be run using python 2.7. The HMM-related scripts require prior installation of numpy from scipy.org. The SVM-related scripts require prior installation of rpy2 from https://bitbucket.org/rpy2/rpy2. Note that the most recent versions of rpy2 do not work with python 2.7, so an older version must be used. Before running the SVM, the e1071 R package must be installed, which can be done by calling ClassSVM.install_e1071 from an interactive python window; this only needs to be done once on each machine.

Summary of the workflow for creating PhyloCSF Candidate Coding Regions:
1. Run PhyloCSF on every codon in genome.
2. Determine parameters for the HMM.
3. Run the HMM to compute the PhyloCSF Regions.
4. Intersect PhyloCSF Regions with known annotations.
5. Run PhyloCSF with the mle option on the PhyloCSF Regions and their antisense regions.
6. Prune and sort the list of PhyloCSF Regions to find PhyloCSF Candidate Coding Regions.

Details:

Step 1 must be completed by the user. Run PhyloCSF with the "strategy=fixed" and "bls" options on the multispecies alignment of every codon in every frame of each strand of every chromosome or scaffold of the genome. The resulting scores (and bls) need to be put in one file per chromosome/strand/frame, all in the same directory. The required file names and format is described in the header comment for HMMforPhyloCSFRegions.create_PhyloCSF_Regions. Instructions for installing PhyloCSF may be found at https://github.com/mlin/PhyloCSF/wiki. The user will need to extract multispecies alignments for every codon in every frame of the genome, which will typically involve a bespoke procedure for each genome for which PhyloCSF Candidate Coding Regions are being computed. Step 1 is computationally intensive. For example, it required more than 10,000 CPU-hours for the human genome.

Step 2. Compute initial state, transition, and emission probabilities for the HMM by calling EstimateHMMparams.estimate_hmm_params_for_genome with user-supplied information about the genome, namely the total length as well as information about every coding exon (described in the comment at the head of the function). The parameters that were actually used in the paper are included in HMMparamsForMudge2019.py.

Step 3. For each chromosome or scaffold, each frame, and both strands, call HMMforPhyloCSFRegions.create_PhyloCSF_Regions with the HMM parameters computed in step 2 and the PhyloCSF output from step 1 to compute the intervals, "PhyloCSF Regions", in the coding state in the most likely path through the HMM.

Step 4. Call SVMsForPhyloCSFCandidateCodingRegions.classify_regions, supplying the path to the PhyloCSF Region files produced in step 3, and user-supplied bed-format files containing the coding and pseudogene annotations for the genome. This will intersect the PhyloCSF Regions with the annotations, produce an initial spreadsheet, and produce the input file for running PhyloCSF in step 5.

Step 5 must be completed by the user. Run PhyloCSF with "strategy=mle" and "bls" on the alignments of the intervals in the file Regions.pcsf.in produced by step 4, which are the PhyloCSF Regions and their antisense regions. Output the results in Regions.pcsf.out.

Step 6. Call SVMsForPhyloCSFCandidateCodingRegions.do_other_steps. This trains SVMs to distinguish regions most likely to be coding; eliminates regions overlapping annotations, short regions, and regions more likely to be antisense; and sorts the remaining regions, "PhyloCSF Candidate Coding Regions", according to how likely they are to be true coding regions. The result is a bed-format file, PhyloCSFNovel.bed, containing the sorted list of PhyloCSF Candidate Coding Regions, and a spreadsheet, Regions.04.txt, containing the associated scores.

The script Examples.py has an example of running through steps 2, 3, 4, and 6. 
Run:
	python Examples.py hmm
	python Examples.py svm

Examples.py also has an example of splice site scoring.
Run:
	python Examples.py splice
