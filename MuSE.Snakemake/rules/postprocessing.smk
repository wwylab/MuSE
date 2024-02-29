def get_all_individual_maf():
    maf_list = []
    if not os.path.exists('SNVCalling/IndividialMAF'):
        return maf_list

    for f in listdir('SNVCalling/IndividialMAF/'):
        if f.endswith('.maf') and len(f.split('_')) == 3:
            maf_list.append(f)
    return maf_list

rule postprocessing:
    input: expand("SNVCalling/IndividialMAF/{all}", all = get_all_individual_maf())
    params:
        ngs=ngs_data_type
    output: "SNVCalling/FinalMAF/final.maf"
    run:
        output_handle = open(output[0], 'w')
        maf_header = ''
        for index, f in enumerate(input):
            with open(f, 'r') as handle:
                for line in handle:
                    line = line.strip()
                    if line.startswith('Hugo_Symbol') and index == 0:
                        output_handle.write(line + '\n')
                    elif line.startswith('Hugo_Symbol'):
                        continue
                    else:
                        output_handle.write(line + '\n')
        output_handle.close()