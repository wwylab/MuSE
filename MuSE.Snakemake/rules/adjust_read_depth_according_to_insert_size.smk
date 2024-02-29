rule adjust_read_depth:
    input:
        merge_vcf="SNVCalling/RawMerge/{patient}_{sample1}_{sample2}.vcf",
        allele_counter_loci="SNVCalling/RawMerge/{patient}_{sample1}_{sample2}.loci",
        gatk4_normal_bam = "Preprocessing/{patient}/{patient}_{sample1}_output-GATK.bam",
        gatk4_tumor_bam = "Preprocessing/{patient}/{patient}_{sample2}_output-GATK.bam",
        reference_genome = expand("{root}/{version}/{genome}", root=config["ref"]["root"], version=config["ref"]["version"], genome=config["ref"]["genome"])
    output:
        dp_adjusted_vcf="SNVCalling/RawMergeAdjustedDP/{patient}_{sample1}_{sample2}.vcf",
    params:
        normal_read_count="SNVCalling/RawMergeAdjustedDP/{patient}_{sample1}_{sample2}_normal_read_count.txt",
        tumor_read_count="SNVCalling/RawMergeAdjustedDP/{patient}_{sample1}_{sample2}_tumor_read_count.txt"
    resources:
        mem_mb = resource_alloc['mem_mb'],
        mem_mbi = resource_alloc['mem_mbi'],
        mem_gb = resource_alloc['mem_gb'],
    threads: resource_alloc['small']
    conda: "../envs/env-allelecount.yaml"
    script:
        "../scripts/run_adjust_dp.py"