def get_all_individual_preqc_multiqc():
    qc_list = []
    if not os.path.exists('QC/PreQC'):
        return qc_list

    for f in listdir('QC/PreQC/'):
        if f.endswith('_fastqc.html'):
            qc_list.append(f.replace('_fastqc.html', ''))
    return qc_list

def get_all_individual_postqc_multiqc():
    qc_list = []
    if not os.path.exists('QC/PostQC'):
        return qc_list

    for f in listdir('QC/PostQC/'):
        if f.endswith('_fastqc.html'):
            qc_list.append(f.replace('_fastqc.html', ''))
    return qc_list

def get_all_individual_addrg_multiqc():
    qc_list = []
    if not os.path.exists('QC/AddRG'):
        return qc_list

    for f in listdir('QC/AddRG/'):
        if f.endswith('summary_metrics.txt'):
            qc_list.append(f)
    return qc_list

def get_all_individual_dedup_multiqc():
    qc_list = []
    if not os.path.exists('QC/MarkDup'):
        return qc_list

    for f in listdir('QC/MarkDup/'):
        if f.endswith('summary_metrics.txt'):
            qc_list.append(f)
    return qc_list

rule multiqc_preqc:
    input: expand("QC/PreQC/{all}_fastqc.html", all = get_all_individual_preqc_multiqc())
    output: "QC/MultiQC/pre_qc.html"
    params:
        input_dir = "QC/PreQC/",
        output_dir = "QC/MultiQC/",
        output_name = "pre_qc.html"
    conda: "../envs/env-fastqc.yaml" 
    resources:
        mem_mb = resource_alloc['mem_mb'],
        mem_mbi = resource_alloc['mem_mbi'],
        mem_gb = resource_alloc['mem_gb'],
    threads: resource_alloc['small']
    shell:
        """
        multiqc --force -o {params.output_dir} -n {params.output_name} {params.input_dir}
        """

rule multiqc_postqc:
    input: expand("QC/PostQC/{all}_fastqc.html", all = get_all_individual_postqc_multiqc())
    output: "QC/MultiQC/post_qc.html"
    params:
        input_dir = "QC/PostQC/",
        output_dir = "QC/MultiQC/",
        output_name = "post_qc.html"
    conda: "../envs/env-fastqc.yaml" 
    resources:
        mem_mb = resource_alloc['mem_mb'],
        mem_mbi = resource_alloc['mem_mbi'],
        mem_gb = resource_alloc['mem_gb'],
    threads: resource_alloc['small']
    shell:
        """
        multiqc --force -o {params.output_dir} -n {params.output_name} {params.input_dir}
        """

rule multiqc_aligned_bam:
    input: expand("QC/AddRG/{all}", all = get_all_individual_addrg_multiqc())
    output: "QC/MultiQC/aligned_bam_qc.html"
    params:
        input_dir = "QC/AddRG/",
        output_dir = "QC/MultiQC/",
        output_name = "aligned_bam_qc.html"
    conda: "../envs/env-fastqc.yaml" 
    resources:
        mem_mb = resource_alloc['mem_mb'],
        mem_mbi = resource_alloc['mem_mbi'],
        mem_gb = resource_alloc['mem_gb'],
    threads: resource_alloc['small']
    shell:
        """
        multiqc --force -o {params.output_dir} -n {params.output_name} {params.input_dir}
        """

rule multiqc_mark_dup:
    input: expand("QC/MarkDup/{all}", all = get_all_individual_dedup_multiqc())
    output: "QC/MultiQC/mark_duplicates_qc.html"
    params:
        input_dir = "QC/MarkDup/",
        output_dir = "QC/MultiQC/",
        output_name = "mark_duplicates_qc.html"
    conda: "../envs/env-fastqc.yaml" 
    resources:
        mem_mb = resource_alloc['mem_mb'],
        mem_mbi = resource_alloc['mem_mbi'],
        mem_gb = resource_alloc['mem_gb'],
    threads: resource_alloc['small']
    shell:
        """
        multiqc --force -o {params.output_dir} -n {params.output_name} {params.input_dir}
        """
