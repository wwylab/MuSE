# A fully-automated pipeline for fast consensus mutation calling from sequencing BAM files

This pipeline efficiently generates consensus somatic calls for next-generation sequencing (NGS) data from a cancer patient cohort. It begins by processing BAM or FASTQ files from tumor-normal pairs through a sequence of finetuning steps. Following this, it performs consensus somatic SNV calls by intersecting the calls from two accelerated methods, MuSE 2<sup>1</sup> and Strelka2<sup>2</sup>. Additional post-processing is then carried out to finalize the results. The pipeline is designed to optimize CPU and memory use on a High-Performance Computing (HPC) environment by parallelizing independent tasks. This pipeline eliminates the manual curation of each task, hence lowering the complexity for users, such as clinicians, who seek quick access to variant calls from a cohort but may lack extensive expertise in bioinformatics or computational biology. The latest version of this pipeline is **v1.0**.

## Download the pipeline
The user can download the pipeline from the package of MuSE 2. 

``` 
git clone https://github.com/wwylab/MuSE.git
cd MuSE.Snakemake
```

## Software dependencies

1. **Snakemake (v7.0 or above)** install it from https://snakemake.readthedocs.io/en/stable/getting_started/installation.html.
2. **MuSE 2**: install it following the README: https://github.com/wwylab/MuSE.
3. **VEP (v101 or above)**: install it and set up the reference files following https://useast.ensembl.org/info/docs/tools/vep/script/index.html.
4. **vcf2maf (v1.6.18 or above)**: install it from https://github.com/mskcc/vcf2maf.

Other required software are automatically installed by the pipeline. 

## Reference files

This pipeline requires the following reference files to run:
1. Indexed reference genome FASTA file.
2. database of known SNPs (VCF format) generated based on the same reference genome as (1).
3. database of known indels (VCF format) generated based on the same reference genome as (1).

We suggest to download them from the Broad Institute Resource Bundle, and save them in the same folder: 

1. For genome build hg38, please go to https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0. Download these files:
   - Homo_sapiens_assembly38.fasta
   - Homo_sapiens_assembly38.dict
   - Homo_sapiens_assembly38.fasta.fai
   - Homo_sapiens_assembly38.dbsnp138.vcf
   - Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
   - Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi

2. For genome build hg19, please go to https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg19/v0. Download these files:
   - Homo_sapiens_assembly19.fasta
   - Homo_sapiens_assembly19.dict
   - Homo_sapiens_assembly19.fasta.fai
   - Homo_sapiens_assembly19.dbsnp.vcf
   - Mills_and_1000G_gold_standard.indels.b37.vcf.gz
   - Mills_and_1000G_gold_standard.indels.b37.vcf.gz.tbi

**Note:** For Homo_sapiens_assembly38.dbsnp138.vcf and Homo_sapiens_assembly19.dbsnp.vcf, please use the following commands to compress and index:

```
bgzip -c Homo_sapiens_assembly38.dbsnp138.vcf > Homo_sapiens_assembly38.dbsnp138.vcf.gz
tabix -p vcf Homo_sapiens_assembly38.dbsnp138.vcf.gz
```
or 
```
bgzip -c Homo_sapiens_assembly19.dbsnp.vcf > Homo_sapiens_assembly19.dbsnp.vcf.gz
tabix -p vcf Homo_sapiens_assembly19.dbsnp.vcf.gz
```

Additionaly, Strelka2 requires a bed file to specific the contigs to call mutations. One can download it from here: hg38 (download both [hg38.bed.gz](https://drive.google.com/file/d/1vrZuTrkRfP6e1agexpHJdST-JZpRmpjc/view?usp=sharing) and [hg38.bed.gz.tbi](https://drive.google.com/file/d/1PXq-AnqUmZHNfPpxfMwFed0D3TkU6pOS/view?usp=sharing)), hg19 (download both [hg19.bed.gz](https://drive.google.com/file/d/1kgpFMnw2h8duU7ts2DHFj3Ksewovv5cb/view?usp=sharing) and [hg19.bed.gz.tbi](https://drive.google.com/file/d/1yzb4K9J7ignDBCWzNBDJJmJpSpn886c5/view?usp=sharing)). Keep them in the same folder as the reference files.


## Input file
This pipeline requires the user to create a `project_info.tsv` file to provide the information about the data to run. It has the following format:

If the raw data input is BAM:

| PatientID | TumorName | TumorPath | NormalName | NormalPath | DataType |
| ---------- | --------- | ---------- | --------- | ---------- | --------- |
| P1  | tumor  | ABSOLUTE_PATH_tumor1.bam  |  normal | ABSOLUTE_PATH_normal1.bam  |  bam |
| P2  | tumor  | ABSOLUTE_PATH_tumor2.bam  |  normal | ABSOLUTE_PATH_normal2.bam  |  bam |

If the raw data input is FASTQ from pair-end sequencing:

| PatientID | TumorName | TumorPath | NormalName | NormalPath | DataType |
| ---------- | --------- | ---------- | --------- | ---------- | --------- |
| P1  | tumor  | ABSOLUTE_PATH_tumor1_R1.fq.gz,ABSOLUTE_PATH_tumor1_R2.fq.gz  |  normal | ABSOLUTE_PATH_normal1_R1.fq.gz,ABSOLUTE_PATH_normal1_R2.fq.gz  |  fastq |
| P2  | tumor  | ABSOLUTE_PATH_tumor2_R1.fq.gz,ABSOLUTE_PATH_tumor2_R2.fq.gz  |  normal | ABSOLUTE_PATH_normal2_R1.fq.gz,ABSOLUTE_PATH_normal2_R2.fq.gz  |  fastq |

## Environment configuration

Before running the pipeline, we need to set the paths of the reference files and `project_info.tsv` in the configuration file `config/config.yaml`. 

- `workdir`: the path of the root directory for running this pipeline. 

- `project_info`: the path of `project_info.tsv`

Since we run this pipeline on a HPC, we need to set up a profile to specify how Snakemake behaves during the job running: how many jobs can be run at the same time? If  Snakemake can use conda or singularity to create virtual environments? How many times to restart if a job fails? The user can follow this tutorial https://github.com/Snakemake-Profiles/lsf#snakemake-lsf-profile. 

An example setup:

```
LSF_UNIT_FOR_LIMITS: GB
UNKWN_behaviour: wait
ZOMBI_behaviour: ignore
latency_wait: 5
use_conda: True
using_singularity: True
restart_times: 2
print_shell_command: True
jobs: 500
default_mem_mb: 1024
default_cluster_logdir: default
default_queue: None
default_project: None
max_status_checks_per_second: 10
max_jobs_per_second: 100
max_status_checks: 1
wait_between_tries: 0.001
profile_name: lsf
```

These settings will be saved to the home directory: `~/.config/snakemake/lsf/`

## Run the pipeline

Go to the root directory of this pipeline. Run the following command to start the pipeline.

```nohup bash snakemake_jobscript.lsf > log.txt &```

You can find the job running log in the `log.txt`. The detailed log files of each step in the pipeline are saved in the `workdir/logs/cluster/`.

## Output

The intermediate and final files are stored in different folders under the `workdir` directory:
  - `Preprocessing`: intermediate files generated during the preprocessing step (e.g., bam2fastq, triming, alignment, marking
duplicates, coalignment and base quality recalibration).
  - `SNVCalling`: this folder saves the following files:
    - Calls from each individual caller - MuSE 2 and Strelka2.
    - Consensus calls by intersecting MuSE 2 and Strelka2 calls
    - Consensus calls annotated by VEP<sup>3</sup>
    - A final MAF that combines the consensus calls and the VEP based annotations for all the patient samples in the cohort.
  - `QC`: 
    - fastqc reports for the fastq files before and after triming.
    - Read alignment summaries for the BAM files after the alignment and after marking duplicates.

There is also a `Raw` folder, which includes all the raw input BAM/FASTQ data, with the ones from the same patient saved in the same subfolder. The pipeline only creates symbolic links for these files pointing to the original locations.

**Note:** This pipeline was developed based on our benchmarking results from Ji et al.<sup>1</sup>, which works for both WES and WGS data. If the user is processing WES data and wants to save computing resources and time, he/she can comment out the following lines in the code to output the mutation calling result of MuSE 2, which was shown to be the same accuracy level as the consensus calls. For WES data, we recommend filtering out SNVs at the lowest quality tier - `Tier5` from the vcf output of MuSE 2(check https://github.com/wwylab/MuSE).  

```
## line 218-221 in the Snakemake file
output_file_collection.append("SNVCalling/Strelka2/%s_%s_%s.vcf.gz" % (patient, "normal", "tumor"))
output_file_collection.append("SNVCalling/RawMerge/%s_%s_%s.vcf" % (patient, "normal", "tumor"))
output_file_collection.append("SNVCalling/RawMergeAdjustedDP/%s_%s_%s.vcf" % (patient, "normal", "tumor")) 
output_file_collection.append("SNVCalling/IndividialMAF/%s_%s_%s.maf" % (patient, "normal", "tumor"))

## line 227 in the Snakemake file
output_file_collection.append("SNVCalling/FinalMAF/final.maf")
```


## Reference

1. <ins>Ji S</ins>, Zhu T, Sethia A, **Wang W**. Accelerated somatic mutation calling for whole-genome and whole-exome sequencing data from heterogenous tumor samples. Genome Res. 2024 May 3;. doi: 10.1101/gr.278456.123.

2. Kim S, et al. Strelka2: fast and accurate calling of germline and somatic variants. Nature Methods. 2018 Aug 15. 591-594. doi: 10.1038/s41592-018-0051-x.

3. McLaren W, et al. The Ensembl Variant Effect Predictor. Genome Biology. 2016 Jun 6. 1-14. doi: 10.1186/S13059-016-0974-4/TABLES/8.

