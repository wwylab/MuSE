rule bam2fastq:
    input:
        raw_bam = expand("{raw_data}/{{patient}}/{{patient}}_{{sample}}.bam", raw_data=config["raw_data"])
    output:
        fastq_r1 = temp("Preprocessing/{patient}/{patient}_{sample}_R1.fq.gz"),
        fastq_r2 = temp("Preprocessing/{patient}/{patient}_{sample}_R2.fq.gz"),
	    fastq_s = temp("Preprocessing/{patient}/{patient}_{sample}.fq.gz")
    params:
        fastq_r1="{patient}_{sample}_R1.fq.gz",
        fastq_r2="{patient}_{sample}_R2.fq.gz",
        fastq_s="{patient}_{sample}.fq.gz",
        output_dir="Preprocessing/{patient}",
        fastq_suffix="Preprocessing/{patient}/{sample}"
    resources:
        mem_mb = resource_alloc['mem_mb'],
        mem_mbi = resource_alloc['mem_mbi'],
        mem_gb = resource_alloc['mem_gb'],
    threads: resource_alloc['small']
    conda: "../envs/env-biobambam.yaml"
    shell:
        """
 bamtofastq collate=1 exclude=QCFAIL,SECONDARY,SUPPLEMENTARY filename={input.raw_bam} gz=1 inputformat=bam level=5 outputperreadgroup=0 F={output.fastq_r1} F2={output.fastq_r2} O={params.fastq_suffix}_o1.fq.gz O2={params.fastq_suffix}_o2.fq.gz S={output.fastq_s}
        """

rule fastqc_before_trimming:
    input:
        fastq_r1 = "Preprocessing/{patient}/{patient}_{sample}_R1.fq.gz",
        fastq_r2 = "Preprocessing/{patient}/{patient}_{sample}_R2.fq.gz"
    output:
        fastqc_r1_html = "QC/PreQC/{patient}_{sample}_R1_fastqc.html",
        fastqc_r1_zip = "QC/PreQC/{patient}_{sample}_R1_fastqc.zip",
        fastqc_r2_html = "QC/PreQC/{patient}_{sample}_R2_fastqc.html",
        fastqc_r2_zip = "QC/PreQC/{patient}_{sample}_R2_fastqc.zip"      
    resources:
        mem_mb = resource_alloc['mem_mb'],
        mem_mbi = resource_alloc['mem_mbi'],
        mem_gb = resource_alloc['mem_gb'],
    threads: resource_alloc['small']
    params:
        output_dir = "QC/PreQC/"
    conda: "../envs/env-fastqc.yaml"
    shell:
        """
fastqc {input.fastq_r1} -t {threads} -o {params.output_dir}
fastqc {input.fastq_r2} -t {threads} -o {params.output_dir}
        """

rule trim_reads:
    input:
        fastq_r1 = "Preprocessing/{patient}/{patient}_{sample}_R1.fq.gz",
        fastq_r2 = "Preprocessing/{patient}/{patient}_{sample}_R2.fq.gz"
    output:
        fastq_val_r1 = temp("Preprocessing/{patient}/{patient}_{sample}_R1_val_1.fq.gz"),
        fastq_val_r2 = temp("Preprocessing/{patient}/{patient}_{sample}_R2_val_2.fq.gz")
    resources:
        mem_mb = resource_alloc['mem_mb'],
        mem_mbi = resource_alloc['mem_mbi'],
        mem_gb = resource_alloc['mem_gb'],
    threads: resource_alloc['medium']
    log:
        "logs/trim_reads/{patient}{sample}.log"
    params: 
        output_dir = "Preprocessing/{patient}"
    conda: "../envs/env-trim_glore.yaml"
    shell:
        """
        cd {params.output_dir}
        trim_galore -j {threads} --paired ../../{input.fastq_r1} ../../{input.fastq_r2}
        cd ../../
        """

rule fastqc_after_trimming:
    input:
        fastq_val_r1 = "Preprocessing/{patient}/{patient}_{sample}_R1_val_1.fq.gz",
        fastq_val_r2 = "Preprocessing/{patient}/{patient}_{sample}_R2_val_2.fq.gz"
    output:
        fastqc_r1_html = "QC/PostQC/{patient}_{sample}_R1_val_1_fastqc.html",
        fastqc_r1_zip = "QC/PostQC/{patient}_{sample}_R1_val_1_fastqc.zip",
        fastqc_r2_html = "QC/PostQC/{patient}_{sample}_R2_val_2_fastqc.html",
        fastqc_r2_zip = "QC/PostQC/{patient}_{sample}_R2_val_2_fastqc.zip"      
    resources:
        mem_mb = resource_alloc['mem_mb'],
        mem_mbi = resource_alloc['mem_mbi'],
        mem_gb = resource_alloc['mem_gb'],
    threads: resource_alloc['small']
    conda: "../envs/env-fastqc.yaml"
    params:
        output_dir = "QC/PostQC/"
    shell:
        """
fastqc {input.fastq_val_r1} -t {threads} -o {params.output_dir}
fastqc {input.fastq_val_r2} -t {threads} -o {params.output_dir}
        """

rule align_reads:
    input: 
        fastq_val_r1 = "Preprocessing/{patient}/{patient}_{sample}_R1_val_1.fq.gz",
        fastq_val_r2 = "Preprocessing/{patient}/{patient}_{sample}_R2_val_2.fq.gz",
        reference_genome = expand("{root}/{version}/{genome}", root=config["ref"]["root"], version=config["ref"]["version"], genome=config["ref"]["genome"]), 
        reference_genome_faidx = expand("{root}/{version}/{genome}.fai", root=config["ref"]["root"], version=config["ref"]["version"], genome=config["ref"]["genome"])      
    output:
        align_bam = temp("Preprocessing/{patient}/{patient}_{sample}_realigned.bam")
    resources:
        mem_mb = resource_alloc['mem_mb'],
        mem_mbi = resource_alloc['mem_mbi'],
        mem_gb = resource_alloc['mem_gb'],
    threads: resource_alloc['large']
    log:
        "logs/alin_reads/{patient}{sample}.log"
    params:
        rg = r"@RG\tID:{sample}\tSM:{sample}\tLB:lib\tPL:ILLUMINA"
    conda: "../envs/env-bwa.yaml"
    shell:
        """
        bwa mem -t {threads} -T 0 -R '{params.rg}' {input.reference_genome} {input.fastq_val_r1} {input.fastq_val_r2} | samtools view -@ 2 -Shb -o {output.align_bam}
        """

rule sort_bam:
    input: 
        align_bam = "Preprocessing/{patient}/{patient}_{sample}_realigned.bam"
    output:
        sorted_bam = temp("Preprocessing/{patient}/{patient}_{sample}_realign-sort.bam")
    log:
        "logs/alin_reads/{patient}{sample}.log"
    params:
        extra="-m 2G",
    resources:
        mem_mb = resource_alloc['mem_mb'],
        mem_mbi = resource_alloc['mem_mbi'],
        mem_gb = resource_alloc['mem_gb'],
    threads: resource_alloc['small']
    conda: "../envs/env-samtools.yaml"
    shell: 
        """
        samtools sort -@ {threads} {params.extra} -o {output.sorted_bam} {input.align_bam}
        """

rule collect_alignment_summary_after_alignment:
    input: 
        sorted_bam = "Preprocessing/{patient}/{patient}_{sample}_realign-sort.bam",
        reference_genome = expand("{root}/{version}/{genome}", root=config["ref"]["root"], version=config["ref"]["version"], genome=config["ref"]["genome"])
    output:
        alignment_summary = "QC/AddRG/{patient}_{sample}_summary_metrics.txt"    
    resources:
        mem_mb = resource_alloc['mem_mb'],
        mem_mbi = resource_alloc['mem_mbi'],
        mem_gb = resource_alloc['mem_gb'],
    threads: resource_alloc['small']
    conda: "../envs/env-picard.yaml"
    shell:
        "picard CollectAlignmentSummaryMetrics R={input.reference_genome} I={input.sorted_bam} O={output.alignment_summary}"


rule mark_duplicates:
    input: 
        sorted_bam = "Preprocessing/{patient}/{patient}_{sample}_realign-sort.bam"
    output:
        dedup_bam = temp("Preprocessing/{patient}/{patient}_{sample}_realign-sort-rg-dedup.bam"),
        dedup_metric = temp("Preprocessing/{patient}/{patient}_{sample}_marked_dup_metrics.txt")
    params:
        tmp_dir="Preprocessing/{patient}/tmp_{patient}_{sample}"
    resources:
        mem_mb = resource_alloc['mem_mb'],
        mem_mbi = resource_alloc['mem_mbi'],
        mem_gb = resource_alloc['mem_gb'],
    threads: resource_alloc['small']
    log:
        "logs/mark_duplicates_gatk3/{patient}{sample}.log"
    conda: "../envs/env-picard.yaml"
    shell:
        """
        rm -rf {params.tmp_dir}
        mkdir -p {params.tmp_dir}
        picard MarkDuplicates I={input.sorted_bam} O={output.dedup_bam} M={output.dedup_metric} TMP_DIR={params.tmp_dir} VALIDATION_STRINGENCY=STRICT CREATE_INDEX=true
        """  

rule collect_alignment_summary_after_mark_duplicates:
    input: 
        dedup_bam = "Preprocessing/{patient}/{patient}_{sample}_realign-sort-rg-dedup.bam",
        reference_genome = expand("{root}/{version}/{genome}", root=config["ref"]["root"], version=config["ref"]["version"], genome=config["ref"]["genome"])
    output:
        alignment_summary = "QC/MarkDup/{patient}_{sample}_summary_metrics.txt"    
    resources:
        mem_mb = resource_alloc['mem_mb'],
        mem_mbi = resource_alloc['mem_mbi'],
        mem_gb = resource_alloc['mem_gb'],
    threads: resource_alloc['small']
    conda: "../envs/env-picard.yaml"
    shell:
        "picard CollectAlignmentSummaryMetrics R={input.reference_genome} I={input.dedup_bam} O={output.alignment_summary}"

rule coalign_target_interval_gatk3:
    input:
        dedup_bam = lambda wildcards: expand("Preprocessing/{{patient}}/{{patient}}_{sample}_realign-sort-rg-dedup.bam", sample=get_samples_for_patient(wildcards)),
        reference_genome = expand("{root}/{version}/{genome}", root=config["ref"]["root"], version=config["ref"]["version"], genome=config["ref"]["genome"]),
        reference_known_indel = expand("{root}/{version}/{known_indel}", root=config["ref"]["root"], version=config["ref"]["version"], known_indel=config["ref"]["known_indel"])
    output:
        realign_target_interval = temp("Preprocessing/{patient}/{patient}_realign-targets.intervals")
    params:
        dedup_bam = lambda wildcards: ' '.join(["-I {{patient}}_{sample}_realign-sort-rg-dedup.bam".format(sample=sample) for sample in get_samples_for_patient(wildcards)]),
        java_options = "-Xmx%sg" % (resource_alloc['mem_gb'] - 5)
    resources:
        mem_mb = resource_alloc['mem_mb'],
        mem_mbi = resource_alloc['mem_mbi'],
        mem_gb = resource_alloc['mem_gb'],
    threads: resource_alloc['medium']
    log:
        "logs/coalign_target_interval_gatk3/{patient}.log"
    conda: "../envs/env-java.yaml"

    script:
        "../scripts/run_coalign_target_interval.py"


rule coalign_target_bam_gatk3:
    input:
        dedup_bam = "Preprocessing/{patient}/{patient}_{sample}_realign-sort-rg-dedup.bam",
        reference_genome = expand("{root}/{version}/{genome}", root=config["ref"]["root"], version=config["ref"]["version"], genome=config["ref"]["genome"]),
        reference_known_indel = expand("{root}/{version}/{known_indel}", root=config["ref"]["root"], version=config["ref"]["version"], known_indel=config["ref"]["known_indel"]),
        realign_target_interval = "Preprocessing/{patient}/{patient}_realign-targets.intervals"
    output:
        realign_bam = temp("Preprocessing/{patient}/{patient}_{sample}_realign-sort-rg-dedup_coalign.bam")
    params:
        # dedup_bam = lambda wildcards: ' '.join(["-I {{patient}}_{sample}_realign-sort-rg-dedup.bam".format(sample=sample) for sample in get_samples_for_patient(wildcards)]),
        realign_target_interval = "{patient}_realign-targets.intervals",
        output_dir = "Preprocessing/{patient}",
        # java_path = config["software"]["GATK3"],
        java_options = "-Xmx%sg" % (resource_alloc['mem_gb'] - 5),
        bam_compression_level = config["bam_compression_level"]
    log:
        "logs/coalign_target_bam_gatk3/{patient}{sample}.log"
    conda: "../envs/env-java.yaml"
    resources:
        mem_mb = resource_alloc['mem_mb'],
        mem_mbi = resource_alloc['mem_mbi'],
        mem_gb = resource_alloc['mem_gb'],
    threads: resource_alloc['small']
    shell:
        """
        module load gatk/3.7
        java {params.java_options} -jar $gatk37 -T IndelRealigner -R {input.reference_genome} -I {input.dedup_bam} -targetIntervals {input.realign_target_interval} -known {input.reference_known_indel} -o {output.realign_bam} -compress {params.bam_compression_level}
        module unload gatk/3.7
        """

rule base_recalibration_gatk3:
    input: 
        realign_bam = "Preprocessing/{patient}/{patient}_{sample}_realign-sort-rg-dedup_coalign.bam",
        reference_genome = expand("{root}/{version}/{genome}", root=config["ref"]["root"], version=config["ref"]["version"], genome=config["ref"]["genome"]), 
        reference_known_snp = expand("{root}/{version}/{known_snp}", root=config["ref"]["root"], version=config["ref"]["version"], known_snp=config["ref"]["known_snp"])
    output:
        gatk3_rec_table = temp("Preprocessing/{patient}/{patient}_{sample}_output-recal-table-gatk3.txt")
    params:
        java_options = "-Xmx%sg" % (resource_alloc['mem_gb'] - 5)
    log:
        "logs/base_recalibration_gatk3/{patient}{sample}.log"
    resources:
        mem_mb = resource_alloc['mem_mb'],
        mem_mbi = resource_alloc['mem_mbi'],
        mem_gb = resource_alloc['mem_gb'],
    threads: resource_alloc['medium']
    conda: "../envs/env-java.yaml"
    shell:
        """
        module load gatk/3.7
        java {params.java_options} -jar $gatk37 -T BaseRecalibrator -R {input.reference_genome} -I {input.realign_bam} -knownSites {input.reference_known_snp} -o {output.gatk3_rec_table} -nct {threads}
        module unload gatk/3.7
        """       

rule apply_bqsr_gatk3:
    input: 
        realign_bam = "Preprocessing/{patient}/{patient}_{sample}_realign-sort-rg-dedup_coalign.bam",
        reference_genome = expand("{root}/{version}/{genome}", root=config["ref"]["root"], version=config["ref"]["version"], genome=config["ref"]["genome"]), 
        gatk3_rec_table = "Preprocessing/{patient}/{patient}_{sample}_output-recal-table-gatk3.txt"
    output:
        gatk3_bam = "Preprocessing/{patient}/{patient}_{sample}_output-GATK.bam",
        gatk3_bai = "Preprocessing/{patient}/{patient}_{sample}_output-GATK.bai"
    params:
        java_options = "-Xmx%sg" % (resource_alloc['mem_gb'] - 5)
    log:
        "logs/apply_bqsr_gatk3/{patient}{sample}.log"
    resources:
        mem_mb = resource_alloc['mem_mb'],
        mem_mbi = resource_alloc['mem_mbi'],
        mem_gb = resource_alloc['mem_gb'],
    threads: resource_alloc['medium']
    conda: "../envs/env-java.yaml"
    shell:
        """
        module load gatk/3.7
        java {params.java_options} -jar $gatk37 -T PrintReads -R {input.reference_genome} -I {input.realign_bam} -BQSR {input.gatk3_rec_table} -o {output.gatk3_bam} -nct {threads} -compress {config[bam_compression_level]}
        module unload gatk/3.7
        """ 