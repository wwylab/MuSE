rule mutation_calling_muse2:
    input: 
        gatk_normal_bam = "Preprocessing/{patient}/{patient}_{sample1}_output-GATK.bam",
        gatk_tumor_bam = "Preprocessing/{patient}/{patient}_{sample2}_output-GATK.bam",
        reference_genome = expand("{root}/{version}/{genome}", root=config["ref"]["root"], version=config["ref"]["version"], genome=config["ref"]["genome"]),
        reference_known_snp = expand("{root}/{version}/{known_snp}", root=config["ref"]["root"], version=config["ref"]["version"], known_snp=config["ref"]["known_snp"])
    output:
        output_vcf = "SNVCalling/MuSE2/{patient}_{sample1}_{sample2}.vcf.gz"
    params:
        normal_sample = "{sample1}",
        tumor_sample = "{sample2}",
        output_intermediate_vcf = "SNVCalling/MuSE2/{patient}_{sample1}_{sample2}",
    resources:
        mem_mb = resource_alloc['mem_mb'],
        mem_mbi = resource_alloc['mem_mbi'],
        mem_gb = resource_alloc['mem_gb'],
    threads: resource_alloc['large']
    conda: "../envs/env-muse.yaml"
    shell:
        """
MuSE call -O {params.output_intermediate_vcf}_call_output -f {input.reference_genome} -n {threads} {input.gatk_tumor_bam} {input.gatk_normal_bam}
MuSE sump -I {params.output_intermediate_vcf}_call_output.MuSE.txt -O {params.output_intermediate_vcf}.vcf -E -D {input.reference_known_snp}
rm -f {params.output_intermediate_vcf}.vcf.gz
bgzip {params.output_intermediate_vcf}.vcf
        """