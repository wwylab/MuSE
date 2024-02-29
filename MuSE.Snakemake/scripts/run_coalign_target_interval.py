from snakemake.shell import shell
import sys
import os

dedup_bam = snakemake.input.get("dedup_bam", "")
dedup_bam = ' '.join(['-I ' + bam for bam in dedup_bam])

reference_genome = snakemake.input.get("reference_genome", "")
reference_known_indel = snakemake.input.get("reference_known_indel", "")
realign_target_interval = snakemake.output.get("realign_target_interval", "")
java_path = snakemake.params.get("java_path", "")
java_options = snakemake.params.get("java_options", "")

with open('/rsrch3/home/bcb/sji/Document/Programming/snakemake/snakemake_mutation_calling_pipeline_v1_muse/temp.txt', 'w') as handle:
    handle.write(dedup_bam + '\n')
    # handle.write(reference_known_indel + '\n')

shell(
"""
    module load gatk/3.7
    java {java_options} -jar $gatk37 -T RealignerTargetCreator -R {reference_genome} {dedup_bam}  -known {reference_known_indel} -o {realign_target_interval} --num_threads {snakemake.threads}
    module unload gatk/3.7
"""
)
