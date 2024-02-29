from snakemake.shell import shell
import pandas as pd
import sys
import shutil
import os

input_vcf = snakemake.input.get("merge_vcf")
allele_counter_loci = snakemake.input.get("allele_counter_loci")
gatk4_normal_bam = snakemake.input.get("gatk4_normal_bam")
gatk4_tumor_bam = snakemake.input.get("gatk4_tumor_bam")
reference_genome = snakemake.input.get("reference_genome", "")

normal_read_count = snakemake.params.get("normal_read_count", "")
tumor_read_count = snakemake.params.get("tumor_read_count", "")

dp_adjusted_vcf = snakemake.output.get("dp_adjusted_vcf", "")

shell(
    """
    alleleCounter -b {gatk4_normal_bam} -q 0 -l {allele_counter_loci} -o {normal_read_count} -r {reference_genome}
    alleleCounter -b {gatk4_tumor_bam} -q 0 -l {allele_counter_loci} -o {tumor_read_count} -r {reference_genome}
    """
    )

dp_adjusted_normal = dict()
dp_adjusted_tumor = dict()

acgt_order = ['A', 'C', 'G', 'T']

with open(normal_read_count, 'r') as handle:
    for line in handle:
        line = line.strip()
        
        if line.startswith('#'):
            line = line.split('\t')
            if len(line) == 7:
                acgt_order = line[2:6]
                acgt_order = [item.replace('Count_', '') for item in acgt_order]
            continue
        
        line = line.split('\t')
        if len(line) != 7:
            continue
        _id = line[0] + ':' + line[1]
        if _id not in dp_adjusted_normal:
            dp_adjusted_normal[_id] = dict()
        
        _dp = line[2:6]
        dp_adjusted_normal[_id] = {key:value for key, value in zip(acgt_order, _dp)}
        
        
acgt_order = ['A', 'C', 'G', 'T']

with open(tumor_read_count, 'r') as handle:
    for line in handle:
        line = line.strip()
        
        if line.startswith('#'):
            line = line.split('\t')
            if len(line) == 7:
                acgt_order = line[2:6]
                acgt_order = [item.replace('Count_', '') for item in acgt_order]
            continue
        
        line = line.split('\t')
        if len(line) != 7:
            continue
        _id = line[0] + ':' + line[1]
        if _id not in dp_adjusted_tumor:
            dp_adjusted_tumor[_id] = dict()
        
        _dp = line[2:6]
        dp_adjusted_tumor[_id] = {key:value for key, value in zip(acgt_order, _dp)}
        
output_handle = open(dp_adjusted_vcf, 'w')
 
with open(input_vcf, 'r') as handle:
    for line in handle:
        line = line.strip()
        if line.startswith('#'):
            output_handle.write(line + '\n')
            continue
        
        line = line.split('\t')
        _id = line[0] + ':' + line[1]
        _ref = line[3]
        _alt = line[4]
        
        normal_orig_dp = line[-2]
        tumor_orig_dp = line[-1]
        dp = line[7]
        try:
            normal_adjusted = dp_adjusted_normal[_id]
            tumor_adjusted = dp_adjusted_tumor[_id]
            
            normal_adjusted_ref = int(normal_adjusted[_ref])
            normal_adjusted_alt = int(normal_adjusted[_alt])
            tumor_adjusted_ref = int(tumor_adjusted[_ref])
            tumor_adjusted_alt = int(tumor_adjusted[_alt])
            
            dp = 'SOMATIC;DP=%s' % (normal_adjusted_ref + normal_adjusted_alt + tumor_adjusted_ref + tumor_adjusted_alt)
            
            normal_orig_dp = '0/0:%s,%s:%s:%s' % (normal_adjusted_ref, normal_adjusted_alt, round(float(normal_adjusted_alt)/float(normal_adjusted_ref + normal_adjusted_alt), 3), (normal_adjusted_ref + normal_adjusted_alt))
            tumor_orig_dp = '0/1:%s,%s:%s:%s' % (tumor_adjusted_ref, tumor_adjusted_alt, round(float(tumor_adjusted_alt)/float(tumor_adjusted_ref + tumor_adjusted_alt), 3), (tumor_adjusted_ref + tumor_adjusted_alt))
            
            line[7] = dp
            line[-2] = normal_orig_dp
            line[-1] = tumor_orig_dp
        except:
            pass
        
        output_handle.write('\t'.join(line) + '\n')
output_handle.close()
