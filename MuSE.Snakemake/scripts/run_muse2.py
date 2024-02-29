from snakemake.shell import shell
import pandas as pd
import sys
import shutil
import os

gatk3_bam = snakemake.input.get("gatk3_bam")
gatk3_bam = set(gatk3_bam)

reference_genome = snakemake.input.get("reference_genome", "")
reference_known_snp = snakemake.input.get("reference_known_snp", "")

output_vcf = snakemake.output.get("output_vcf", "")
project_info = snakemake.params.get("project_info", "")
threads = snakemake.params.get("threads", 2)
muse2_path = snakemake.params.get("muse2_path", "")
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
    
    raw_muse2_file_list = []
    
    for index, row in data.iterrows():
        tumor_sample = row['TumorName']
        normal_sample = row['NormalName']
        
        tumor_sample_path = "Preprocessing/%s/%s_%s_output-GATK3.bam" % (patient_name, patient_name, tumor_sample)
        normal_sample_path = "Preprocessing/%s/%s_%s_output-GATK3.bam" % (patient_name, patient_name, normal_sample)
        
        if tumor_sample_path not in gatk3_bam:
            continue
        if normal_sample_path not in gatk3_bam:
            continue
        
        raw_muse2_file = "SNVCalling/MuSE2/%s_%s-%s" %(patient_name, normal_sample, tumor_sample)
        
        shell(
            """
            MuSE call -O {raw_muse2_file}_call_output -f {reference_genome} -n {cores_of_each_pair} {tumor_sample_path} {normal_sample_path}
            MuSE sump -I {raw_muse2_file}_call_output.MuSE.txt -O {raw_muse2_file}_sump_out_muse.vcf -E -D {reference_known_snp}
            rm -f {raw_muse2_file}_sump_out_muse.vcf.gz
            bgzip {raw_muse2_file}_sump_out_muse.vcf
            """
        )
        if os.path.isfile(raw_muse2_file + "_sump_out_muse.vcf.gz"):
            if os.stat(raw_muse2_file + "_sump_out_muse.vcf.gz").st_size > 0:
                raw_muse2_file_list.append(raw_muse2_file)
    
    if len(raw_muse2_file_list) == 1:
        shutil.copy(raw_muse2_file_list[0] + "_sump_out_muse.vcf.gz", output_vcf)
    
    else:
        raw_muse2_vcfs = ["I=" + os.path.join(item, "_sump_out_muse.vcf.gz") for item in raw_muse2_file_list]
        raw_muse2_vcfs = ' '.join(raw_muse2_vcfs)
        shell (
            """
            java -jar picard.jar MergeVcfs {raw_muse2_vcfs} O={output_vcf}
            """
        )