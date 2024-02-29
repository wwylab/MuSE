rule vcf2maf:
    input: 
        merge_vcf="SNVCalling/RawMergeAdjustedDP/{patient}_{sample1}_{sample2}.vcf",
        reference_genome = expand("{root}/{version}/{genome}", root=config["ref"]["root"], version=config["ref"]["version"], genome=config["ref"]["genome"])
    output:
        maf="SNVCalling/IndividialMAF/{patient}_{sample1}_{sample2}.maf"
    params:
        patient = "{patient}",
        normal_sample = "{sample1}",
        tumor_sample = "{sample2}",
    resources:
        mem_mb = resource_alloc['mem_mb'],
        mem_mbi = resource_alloc['mem_mbi'],
        mem_gb = resource_alloc['mem_gb'],
    threads: resource_alloc['small']
    shell:
        """
vcf2maf.pl --maf-center MDA --input-vcf {input.merge_vcf} --output-maf {output.maf} --vcf-tumor-id TUMOR --tumor-id {params.tumor_sample} --vcf-normal-id NORMAL \
--normal-id {params.patient}_{params.normal_sample}  --ref-fasta {input.reference_genome} --vep-data /rsrch3/scratch/reflib/REFLIB_data/ensembl-vep-v101 \
--filter-vcf 0 --ncbi-build GRCh38 --vep-path /risapps/rhel7/ensembl-vep/v101/
        """
