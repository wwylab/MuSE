rule mutation_calling_strelka2:
    input: 
        gatk_normal_bam = "Preprocessing/{patient}/{patient}_{sample1}_output-GATK.bam",
        gatk_tumor_bam = "Preprocessing/{patient}/{patient}_{sample2}_output-GATK.bam",
        reference_genome = expand("{root}/{version}/{genome}", root=config["ref"]["root"], version=config["ref"]["version"], genome=config["ref"]["genome"]),
        strelka2_bed = expand("{root}/{version}/{strelka2_bed}", root=config["ref"]["root"], version=config["ref"]["version"], strelka2_bed=config["ref"]["strelka_bed"])
    output:
        output_vcf = "SNVCalling/Strelka2/{patient}_{sample1}_{sample2}.vcf.gz"
    params:
        normal_sample = "{sample1}",
        tumor_sample = "{sample2}",
        raw_strelka2_dir = "SNVCalling/Strelka2/{patient}_{sample1}_{sample2}",
        output_intermediate_vcf = "SNVCalling/Strelka2/{patient}_{sample1}_{sample2}/results/variants/somatic.snvs.vcf.gz",
    resources:
        mem_mb = resource_alloc['mem_mb'],
        mem_mbi = resource_alloc['mem_mbi'],
        mem_gb = resource_alloc['mem_gb'],
    threads: resource_alloc['large']
    conda: "../envs/env-strelka.yaml"
    shell:
        """
rm -rf {params.raw_strelka2_dir}
mkdir -p {params.raw_strelka2_dir}
configureStrelkaSomaticWorkflow.py --normalBam {input.gatk_normal_bam} --tumorBam {input.gatk_tumor_bam} --referenceFasta {input.reference_genome} --runDir {params.raw_strelka2_dir} --callRegions {input.strelka2_bed}
{params.raw_strelka2_dir}/runWorkflow.py -m local -j {threads}
cp {params.output_intermediate_vcf} {output.output_vcf}
        """