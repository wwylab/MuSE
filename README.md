# MuSE
The detection of somatic point mutations is a key component of cancer genomic research, which has been rapidly developing since next-generation sequencing (NGS) technology revealed its potential for describing genetic alterations in cancer. We present MuSE, a novel approach to mutation calling based on the F81 Markov substitution model for molecular evolution, which models the evolution of the reference allele to the allelic composition of the matched tumor and normal tissue at each genomic locus. To improve overall accuracy, we further adopt a sample-specific error model to identify cutoffs, reflecting the variation in tumor heterogeneity among samples.

To compile MuSE, type the following command: make

More information can be found at http://bioinformatics.mdanderson.org/main/MuSE
