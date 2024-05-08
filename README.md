# MuSE

An accurate and ultra-fast somatic mutation calling tool for whole-genome sequencing (WGS) and whole-exome sequencing (WES) data from heterogeneous tumor samples. This tool is unique in accounting for tumor heterogeneity using a sample-specific error model that improves sensitivity and specificity in mutation calling from sequencing data. The latest version of this software is **v2.1**.

## News

- **We are thrilled to share that MuSE 2 is published on Genome Research.** Find the paper at [https://genome.cshlp.org/content/early/2024/05/03/gr.278456.123.long](https://genome.cshlp.org/content/early/2024/05/03/gr.278456.123.long).

- **We are excited to announce the launch of an automated pipeline designed for rapid consensus mutation calling, MuSE.Snakemake.** This pipeline starts with BAM or FASTQ files from tumor-normal pairs of a cancer patient cohort, followed by preprocessing stages for the sequencing reads and the intersection of calls from MuSE 2 and [Strelka2](https://github.com/Illumina/strelka). It also includes postprocessing stages for the read depth adjustment and functional annotation for the calls. This pipeline is optimized for High-Performance Computing environments, reducing manual task curation and complexity, thereby making genetic variant analysis accessible to users like clinicians without deep bioinformatics expertise. Please visit the [README of MuSE.Snakemake](https://github.com/wwylab/MuSE/blob/master/MuSE.Snakemake/README.md) for tutorial. 


## Introduction

Detection of somatic point mutations is a key component of cancer genomics research, which has been rapidly developing since next-generation sequencing (NGS) technology revealed its potential for describing genetic alterations in cancer. We previously launched MuSE 1<sup>1</sup>, a statistical approach for mutation calling based on a Markov substitution model for molecular evolution. It has been used as a major contributing caller in a consensus calling strategy by the TCGA PanCanAtlas project<sup>2</sup> and the ICGC Pan-Cancer Analysis of Whole Genomes (PCAWG) initiative<sup>3</sup>.

We have now released MuSE 2<sup>4</sup>, which is powered by a multi-threaded producer-consumer model and an efficient way of memory allocation. MuSE 2 speeds up 50 times than MuSE 1 and 8-80 times than the other callers adopted in the Genomic Data Commons DNA-seq analysis pipeline, i.e., [MuTect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2), [SomaticSniper](https://gmt.genome.wustl.edu/packages/somatic-sniper/) and [VarScan2](https://varscan.sourceforge.net/). MuSE 2 can reduce the computing time cost of a somatic mutation calling project from ∼40 hours to < 1 hour for WGS data, and from 2-4 hours to ~5 minutes for WES data, from each pair of tumor-normal samples. We also performed a benchmarking study, which suggests combining MuSE 2 and the recently accelerated [Strelka2](https://github.com/Illumina/strelka) can almost fully recover PCAWG consensus mutation calls (based on 4 popular callers), as well as recover a majority of the TCGA consensus mutation calls (based on 5 popular callers). 

## Platform
1.	MuSE 1 supports both Linux system and MacOS.
2.	MuSE 2 only supports Linux system with `gcc=7.0` and `git=2.0` or above.

## Installation
```
git clone https://github.com/wwylab/MuSE.git
cd MuSE
./install_muse.sh
```
The executable file `MuSE` will be generated in the same directory.

A Docker file is also provided in the repository for building and running MuSE 2 in a Docker container. 

## Pre-processing
Before running MuSE, raw WES/WGS data need to be processed with the following software, as outlined in the following flowchart. Please refer to GDC best practice guidelines (https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/) for a detailed description of the pre-processing pipeline.

<img src="etc/preprocessing_flowchart.png" alt="preprocessing_flowchart" width="600"/>

1. `Biobambam` (https://github.com/gt1/biobambam) to convert to FASTQ if the input is a binary sequence alignment map (BAM) file
2. `cutadapt` (https://cutadapt.readthedocs.io/en/stable/) for adapter and low-quality bases trimming
3. `bwa` (http://bio-bwa.sourceforge.net/) for alignment against the reference genome
4. `samtools` (http://www.htslib.org/) to sort and index BAM files
5. `picard` (https://broadinstitute.github.io/picard/) to mark duplicates
6. `GATK 3.7` (https://github.com/broadgsa/gatk) for base quality score recalibration

## Input data
Same as MuSE 1, MuSE 2 requires the following files as input:

1. Indexed reference genome FASTA file.
2. BAM formatted sequence data from the pair of tumor and normal DNA samples that have gone through the pre-processing step.
3. dbSNP variant call format (VCF) file that should be bgzip compressed, tabix indexed, and based on the same reference genome.

## Running MuSE
MuSE runs in two steps. The first step, `MuSE call`, takes Files (1) and (2) as the input. This step carries out pre-filtering and calculating position-specific summary statistics using the Markov substitution model.

```
Usage:   MuSE call [options] tumor.bam matched_normal.bam
Options:
         -f FILE    faidx indexed reference genome file
         -O STR     output file name (suffix '.MuSE.txt' is
                    automatically added)
         -n INT     number of cores specified (default=1)

Example: 
MuSE call -f Reference.Genome -O Output.Prefix -n 20 Tumor.bam Matched.Normal.bam
```

The second step, `MuSE sump`, takes the output file from `MuSE call` and File (3) as the input. This step computes tier-based cutoffs from a sample-specific error model. We provide two options for building the model, one for WES data (option `-E`), and the other for WGS data (option `-G`).

```
Usage:   MuSE sump [options]
Options:
         -I FILE    single input file generated by 'MuSE call'
         -G         input generated from whole genome sequencing data
         -E         input generated from whole exome sequencing data
         -O STR     output file name (VCF format)
         -n int     number of cores specified (default=1)
         -D FILE    dbSNP vcf file that should be bgzip compressed,
                    tabix indexed and based on the same reference
                    genome used in 'MuSE call'

Example:
WGS
MuSE sump -I Output.Prefix.MuSE.txt -O Output.Prefix.vcf -G -n 10 -D dbsnp.vcf.gz

or WES
MuSE sump -I Output.Prefix.MuSE.txt -O Output.Prefix.vcf -E -n 10 -D dbsnp.vcf.gz
```

## Output of MuSE
The final output of MuSE is a VCF file (v4.1) that lists the identified somatic variants along with tiered rankings (in the `FILTER` field) for these mutations. The rankings range from `PASS` which is the highest confidence category, followed by `Tiers 1-5`, with `Tier 5` being the tier at the lowest confidence. The INFO field of the VCF file is always `SOMATIC`.


**Note:**
1. For WGS data, we recommend to include the calls of all categories from MuSE for downstream analysis.
2. For WES data, we recommended to include the calls of all categories except Tier 5 from MuSE for downstream analysis. We also recommend the user to filter out non-targeted calls from MuSE using the exome capture kit BED file, which is usually available in the experiment.

## Validation of installation

We provide downsampled WGS data from a cell line tumor-normal pair COLO829/COLO829BL and the corresponding mutation calls from MuSE 2 for the user to validate if the software is run properly. The BAM files can be downloaded from ([tumor.bam](https://drive.google.com/file/d/1K7RA81LhbI1YHX-Z7JT6joI6eZKVHpzc/view?usp=sharing), [tumor.bam.bai](https://drive.google.com/file/d/1LTuGKw5ZAxRdQhomDEjDodSXmry2-HaY/view?usp=sharing), [normal.bam](https://drive.google.com/file/d/1sOGbKPtlCMNXBG8MzeAv49FMghOFmw2s/view?usp=sharing), [normal.bam.bai](https://drive.google.com/file/d/1k3OWBHViGzp8Wv87X2Pf3B4L3IMpurxT/view?usp=sharing), 130MB in total). The mutation calls from MuSE 2 before any filtering can be found in the example folder: `COLO829_illumina_1_100000_1000000.vcf`. 


## Report an issue/bug
Please follow the [issue report template](https://github.com/wwylab/MuSE/blob/master/.github/ISSUE_TEMPLATE/issue-report.md) to report your issue/bug when use MuSE 2/1, which can help us fix it quickly. 

## Acknowledgement
We thank Mehrzad Samadi and his team from Nvidia Corporation, including Tong Zhu, Timothy Harkins and Ankit Sethia, for their contributions towards implementing accelerating techniques in the ` MuSE call` step in MuSE2.

## Reference

1.  <ins>Fan Y</ins>, Xi L, Hughes DST, Zhang J, Zhang J, Futreal PA, Wheeler DA and **Wang W**. MuSE: accounting for tumor heterogeneity using a sample-specific error model improves sensitivity and specificity in mutation calling from sequencing data. Genome Biology. 2016. 17:178. doi: 10.1186/s13059-016-1029-6.

2. Ellrott K, Bailey MH, Saksena G, Covington KR, Kandoth C, Stewart C, Hess J, Ma S, Chiotti KE, McLellan MD, et al. Scalable Open Science Approach for Mutation Calling of Tumor Exomes Using Multiple Genomic Pipelines. Cell Systems. 2018. 271-281. doi: 10.1016/j.cels.2018.03.002.

3. The ICGC/TCGA Pan-Cancer Analysis of Whole Genomes Consortium. 2020. Pan-cancer analysis of whole genomes. Nature 578: 82–93.  doi: 10.1038/s41586-020-1969-6.

4. <ins>Ji S</ins>, Zhu T, Sethia A, **Wang W**. Accelerated somatic mutation calling for whole-genome and whole-exome sequencing data from heterogenous tumor samples. Genome Res. 2024 May 3;. doi: 10.1101/gr.278456.123.

