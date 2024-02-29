rule merge_snvs:
    input:
        muse2_vcf="SNVCalling/MuSE2/{patient}_{sample1}_{sample2}.vcf.gz",
        strelka2_vcf="SNVCalling/Strelka2/{patient}_{sample1}_{sample2}.vcf.gz"
    output:
        merge_vcf="SNVCalling/RawMerge/{patient}_{sample1}_{sample2}.vcf",
        allele_counter_loci="SNVCalling/RawMerge/{patient}_{sample1}_{sample2}.loci"
    params: 
        ngs=ngs_data_type
    resources:
        mem_mb = resource_alloc['mem_mb'],
        mem_mbi = resource_alloc['mem_mbi'],
        mem_gb = resource_alloc['mem_gb'],
    threads: resource_alloc['small']
    conda: "../envs/env-python.yaml"
    script:
        "../scripts/run_merge_snvs.py"