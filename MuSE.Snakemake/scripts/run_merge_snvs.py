import gzip

muse_f = snakemake.input.get("muse2_vcf", "")
strelka2_f = snakemake.input.get("strelka2_vcf", "")

merge_vcf = snakemake.output.get('merge_vcf', '')
loci = snakemake.output.get('allele_counter_loci', '')
ngs = snakemake.params.get('ngs', 'WES')

chroms_hg38 = ['chr' + str(i + 1) for i in range(22)]
chroms_hg38.append('chrX')
chroms_hg38.append('chrY')
chroms_hg38.append('chrM')

chroms_hg38 = {chrom: num for chrom, num in zip(chroms_hg38, range(25))}
chroms_hg19 = {chrom.replace('chr', ''): num for chrom, num in chroms_hg38.items()}
chroms_hg19[-1] = 'MT'

chrom_version_tag = 'hg38'

muse_set = set() # store muse calls
    
with gzip.open(muse_f, 'rt') as handle:
    for line in handle:
        line = line.strip()
        if line.startswith('#'):
            continue
        line = line.split('\t')
        if len(line) < 7:
            continue
            
        if ngs == "WES" and line[6] == 'Tier5':
            continue
        if chrom_version_tag == 'hg38':
            if line[0] not in chroms_hg38:
                continue
        else:
            if line[0] not in chroms_hg19:
                continue
        if len(line[3]) != 1 or len(line[4]) != 1:
            continue
            
        snv = line[0] + ':' + line[1] + ':' + line[3] + ':' + line[4]
        muse_set.add(snv)


output_vcf = merge_vcf.replace('.vcf.gz', '.vcf')
output_loci_handle = open(loci, 'w')
output_vcf_handle = open(output_vcf, 'w')
        
_entry_header = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'NORMAL', 'TUMOR']
_entry_header = '\t'.join(_entry_header)

vcf_header_predefined = []
vcf_header_predefined.append('##FILTER=<ID=PASS,Description="Site contains at least one allele that passes filters">')
vcf_header_predefined.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>")
vcf_header_predefined.append('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">')
vcf_header_predefined.append('##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions of alternate alleles in the tumor">')
vcf_header_predefined.append('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">')
vcf_header_predefined.append('##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">')
vcf_header_predefined.append('##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Indicates if record is a somatic mutation">')
vcf_header_predefined.append('##source=consensus calls from MuSE2 and Strelka2')                  
vcf_header_predefined.append('#' + _entry_header)

count = 0
with gzip.open(strelka2_f, 'rt') as handle:
    for line in handle:
        line = line.strip()
        if line.startswith('#'):
            if line.startswith('##fileformat'):
                output_vcf_handle.write(line + '\n')
            elif line.startswith('##contig'):
                output_vcf_handle.write(line + '\n')
            continue
        
        line = line.split('\t')
        if len(line) < 7:
            continue
            
        if line[6] != 'PASS':
            continue
            
        ref_allele = line[3]
        alt_allele = line[4]
        
        if len(ref_allele) != 1 or len(alt_allele) != 1:
            continue
            
        if line[0].startswith('chr'):
            if line[0] not in chroms_hg38:
                continue
        else:
             if line[0] not in chroms_hg19:
                continue  
            
        info_normal = line[-2].split(':')
        info_tumor = line[-1].split(':')
            
        info_normal = info_normal[4:]
        info_tumor = info_tumor[4:]
        
        bases = ['A', 'C', 'G', 'T']
        info_normal = [int(item.split(',')[0]) for item in info_normal]
        info_tumor = [int(item.split(',')[0]) for item in info_tumor]
        
        info_normal = {base: count for base, count in zip(bases, info_normal)}
        info_tumor = {base: count for base, count in zip(bases, info_tumor)}
        
        ref_normal = info_normal[ref_allele]
        alt_normal = info_normal[alt_allele]
        dp_normal = ref_normal + alt_normal
        af_normal = round(float(alt_normal)/float(dp_normal), 3)
        
        ref_tumor = info_tumor[ref_allele]
        alt_tumor = info_tumor[alt_allele]
        dp_tumor = ref_tumor + alt_tumor
        
        if dp_tumor < 3:
            continue
            
        af_tumor = round(float(alt_tumor)/float(dp_tumor), 3)
        
        total_dp = dp_normal + dp_tumor
        
        _output_normal = '0/0' + ':' + str(ref_normal) + ',' + str(alt_normal) + ':' + str(af_normal) + ':' + str(dp_normal)
        _output_tumor = '0/1' + ':' + str(ref_tumor) + ',' + str(alt_tumor) + ':' + str(af_tumor) + ':' + str(dp_tumor)
        
        _output = 'DP=%s\t' % (total_dp) + 'GT:AD:AF:DP\t' + _output_normal + '\t' + _output_tumor
        
        snv = line[0] + ':' + line[1] + ':' + line[3] + ':' + line[4]
        
        if snv in muse_set:
            if count == 0:
                for _header in vcf_header_predefined:
                    output_vcf_handle.write(_header + '\n')
                    
                if line[0].startswith('chr'):
                    chrom_version_tag = 'hg38'
                else:
                    chrom_version_tag = 'hg19'
                count += 1
                
            if chrom_version_tag == 'hg38':
                if line[0] not in chroms_hg38:
                    continue
            else:
                if line[0] not in chroms_hg19:
                    continue
                
            _output = [line[0], line[1], '.', line[3], line[4], '.', 'PASS', _output]
            _output = '\t'.join(_output)
            
            output_vcf_handle.write(_output + '\n')
            output_loci_handle.write(line[0] + '\t' + line[1] + '\n')
output_vcf_handle.close()       
output_loci_handle.close()