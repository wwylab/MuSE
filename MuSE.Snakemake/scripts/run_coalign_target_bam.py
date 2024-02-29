from snakemake.shell import shell
import sys
import os

dedup_bam = snakemake.input.get("dedup_bam", "")

dedup_bam_trimmed = []

for bam in dedup_bam:
    bam = bam.split('/')[2:]
    bam = '/'.join(bam)
    dedup_bam_trimmed.append(bam)
    
dedup_bam = ' '.join(['-I ' + bam for bam in dedup_bam_trimmed])

reference_genome = snakemake.input.get("reference_genome", "")
reference_known_indel = snakemake.input.get("reference_known_indel", "")
realign_target_interval = snakemake.params.get("realign_target_interval", "")
java_path = snakemake.params.get("java_path", "")
java_options = snakemake.params.get("java_options", "")
bam_compression_level = snakemake.params.get("bam_compression_level", "")
output_dir = snakemake.params.get("output_dir", "")

os.chdir(output_dir)

shell(
"""
    java {java_options} -jar {java_path} -T IndelRealigner -R {reference_genome} {dedup_bam} -targetIntervals {realign_target_interval} -known {reference_known_indel} -nWayOut _coalign.bam -compress {bam_compression_level}
"""
)

os.chdir("../../")