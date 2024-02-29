rule genome_faidx:
    input:
        expand("{root}/{version}/{genome}", root=config["ref"]["root"], version=config["ref"]["version"], genome=config["ref"]["genome"]),
    output:
        expand("{root}/{version}/{genome}.fai", root=config["ref"]["root"], version=config["ref"]["version"], genome=config["ref"]["genome"]),
    log:
        "logs/genome-faidx.log",
    conda: "../envs/env-samtools.yaml"
    resources:
        mem_mb = resource_alloc['mem_mb'],
        mem_mbi = resource_alloc['mem_mbi'],
        mem_gb = resource_alloc['mem_gb'],
    shell:
        "samtools faidx {input} > {output} 2> {log} "


rule genome_dict:
    input:
        expand("{root}/{version}/{genome}", root=config["ref"]["root"], version=config["ref"]["version"], genome=config["ref"]["genome"]),
    output:
        expand("{root}/{version}/{genome}.dict", root=config["ref"]["root"], version=config["ref"]["version"], genome=config["ref"]["genome"]),
    log:
        "logs/samtools/create_dict.log",
    conda: "../envs/env-samtools.yaml"
    resources:
        mem_mb = resource_alloc['mem_mb'],
        mem_mbi = resource_alloc['mem_mbi'],
        mem_gb = resource_alloc['mem_gb'],
    shell:
        "samtools dict {input} > {output} 2> {log} "