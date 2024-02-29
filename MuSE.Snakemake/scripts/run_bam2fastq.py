from snakemake.shell import shell
import sys
import os


raw_bam = snakemake.params.get("raw_bam", "")
if_input_bam_symlink = os.path.islink(raw_bam)

if if_input_bam_symlink:
    raw_bam_source_path = os.readlink(raw_bam)
    # raw_bam_temp: destination of copying the source file to
    raw_bam_temp = raw_bam + '.bam'
    
    os.system("cp %s %s" % (raw_bam_source_path, raw_bam_temp))
    
    # shell(
    # """
    # bamtofastq collate=1 exclude=QCFAIL,SECONDARY,SUPPLEMENTARY filename={raw_bam_temp} gz=1 inputformat=bam level=5  \
    # outputperreadgroup=1 \
    # outputdir={snakemake.params.output_dir} \
    # outputperreadgroupsuffixF={snakemake.params.fastq_suffix}_R1.fq.gz \
    # outputperreadgroupsuffixF2={snakemake.params.fastq_suffix}_R2.fq.gz \
    # tryoq=1

    # mv {snakemake.params.output_dir}/*{snakemake.params.fastq_suffix}_R1.fq.gz {snakemake.output.fastq_r1}
    # mv {snakemake.params.output_dir}/*{snakemake.params.fastq_suffix}_R2.fq.gz {snakemake.output.fastq_r2}
    # # rm {raw_bam_temp}
    # """
    # )

else:
    shell(
    """
    bamtofastq collate=1 exclude=QCFAIL,SECONDARY,SUPPLEMENTARY filename={snakemake.params.raw_bam} gz=1 inputformat=bam level=5  \
    outputperreadgroup=1 \
    outputdir={snakemake.params.output_dir} \
    outputperreadgroupsuffixF={snakemake.params.fastq_suffix}_R1.fq.gz \
    outputperreadgroupsuffixF2={snakemake.params.fastq_suffix}_R2.fq.gz \
    tryoq=1

    mv {snakemake.params.output_dir}/*{snakemake.params.fastq_suffix}_R1.fq.gz {snakemake.output.fastq_r1}
    mv {snakemake.params.output_dir}/*{snakemake.params.fastq_suffix}_R2.fq.gz {snakemake.output.fastq_r2}
    """
    )
