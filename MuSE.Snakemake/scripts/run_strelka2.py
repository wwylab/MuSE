from snakemake.shell import shell
import pandas as pd
import sys
import shutil
import os

gatk4_bam = snakemake.input.get("gatk4_bam")
gatk4_bam = set(gatk4_bam)

reference_genome = snakemake.input.get("reference_genome", "")
strelka2_bed_file = snakemake.input.get("strelka2_bed", "")

output_vcf = snakemake.output.get("output_vcf", "")
project_info = snakemake.params.get("project_info", "")
threads = snakemake.params.get("threads", "")
patient_name = output_vcf.split('/')[-1].replace('.vcf.gz', '').strip()

threads = int(threads)

if patient_name == '':
    sys.exit(-1)
else:
    data = pd.read_table(project_info)
    data = data[data['PatientID'] == patient_name]
    num_of_pairs = data.shape[0]
    cores_of_each_pair = int(threads / num_of_pairs)
    if cores_of_each_pair < 1: cores_of_each_pair = 1
    
    raw_strelka2_dir_list = []
    
    for index, row in data.iterrows():
        tumor_sample = row['TumorName']
        normal_sample = row['NormalName']
        
        tumor_sample_path = "Preprocessing/%s/%s_%s_output-GATK4.bam" % (patient_name, patient_name, tumor_sample)
        normal_sample_path = "Preprocessing/%s/%s_%s_output-GATK4.bam" % (patient_name, patient_name, normal_sample)
        
        if tumor_sample_path not in gatk4_bam:
            continue
        if normal_sample_path not in gatk4_bam:
            continue
        
        raw_strelka2_dir = "SNVCalling/Strelka2/%s_%s-%s" %(patient_name, normal_sample, tumor_sample)
        
        shell(
            """
            rm -rf {raw_strelka2_dir}
            mkdir -p {raw_strelka2_dir}
            configureStrelkaSomaticWorkflow.py --normalBam {normal_sample_path} --tumorBam {tumor_sample_path} --referenceFasta {reference_genome} --runDir {raw_strelka2_dir} --callRegions {strelka2_bed_file}
            {raw_strelka2_dir}/runWorkflow.py -m local -j {cores_of_each_pair}
            """
        )
        if os.path.isfile(os.path.join(raw_strelka2_dir, "results/variants/somatic.snvs.vcf.gz")):
            if os.stat(os.path.join(raw_strelka2_dir, "results/variants/somatic.snvs.vcf.gz")).st_size > 0:
                raw_strelka2_dir_list.append(raw_strelka2_dir)
    
    if len(raw_strelka2_dir_list) == 1:
        shutil.copy(os.path.join(raw_strelka2_dir_list[0], "results/variants/somatic.snvs.vcf.gz"), output_vcf)
    
    else:
        raw_strelka2_vcfs = ["I=" + os.path.join(item, "results/variants/somatic.snvs.vcf.gz") for item in raw_strelka2_dir_list]
        raw_strelka2_vcfs = ' '.join(raw_strelka2_vcfs)
        shell (
            """
            java -jar picard.jar MergeVcfs {raw_strelka2_vcfs} O={output_vcf}
            """
        )