#!/bin/bash
#BSUB -J PROJECT
#BSUB -W 240:00
#BSUB -eo logs/%J.stderr                                    
#BSUB -q long                                                                                                          
#BSUB -n 1                                                         
#BSUB -M 4                                                               
#BSUB -R rusage[mem=4]                                                   
#BSUB -u XX@XX

module load vcf2maf/1.6.18
module load MuSE2
module load snakemake/7.2.0

snakemake \
  --profile lsf \
  --keep-going \
  --restart-times 1 \
  --rerun-incomplete \
  --use-conda
