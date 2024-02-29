from snakemake.shell import shell
import os

ngs = snakemake.params.get('ngs', 'WES')
pon = snakemake.params.get('mutect2_pon')
input_maf = snakemake.params.get('input_maf')
