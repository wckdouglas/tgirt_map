import re

def atropos_trimming(config, input, output, params):
    ''' 
    atropos detected:
        read1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
        read2: GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT
        
    Include barcode and P5/P7:
        read1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG
        read2: GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT
    '''
    R2R = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG'
    R1R = 'GATCGTCGGACTGTAGAACTCTGAACGTGTAGA'
    R2 = 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
    if config['TTN']:
        option='-U 1'
        R2R = 'A' + R2R 
        R2 = re.sub('TCT$','TCTT',R2)
    else:
        option = ''

    #if R2R jumps to R2 RNA, template switch byproduct
    fwd_byproduct = R2[-14:]
    rvs_byproduct = R2R[:14]

    
    #if R2R jimps to R2 RNA or R2R DNA
    fwd_byproduct += ' -b GCACACGTCTGAACTCCAGTCAC -b {R2} '.format(R2 = R2)
    rvs_byproduct += ' -B GTGACTGGAGTTCAGACGTGTGC -b {R2R} '.format(R2R = R2R)

    if config['polyA']:
        smart_seq_CDS = 'AAGCAGTGGTATCAACGCAGAGTAC'
        switch_oligo = 'AGTGGTATCAACGCAGAGTACGGGG'

        fwd_byproduct += ' -a A{100} -a T{100} -g %s -g %s ' %( smart_seq_CDS, switch_oligo)
        rvs_byproduct += ' -A A{100} -A T{100} -G %s -G %s ' %( smart_seq_CDS, switch_oligo)


    single_end_adaptor = '--adapter={R2R} '.format(R2R=R2R)
    paired_end_adaptor = single_end_adaptor + \
            '-A {R1R} '.format(R1R=R1R)
    shared_options = '--minimum-length=15 --threads={threads} --no-cache-adapters '\
                    '--batch-size 100000 --process-timeout 300 '\
                    '--pair-filter=both  --report-file {prefix}.txt --op-order WCGQA '\
            .format(threads = config['threads'], prefix = output['FQ1'].split('.')[0])
    if not config['trim_aggressive']:
        shared_options += '--error-rate=0.1 --overlap 5 --quality-cutoff=20  --aligner insert '

    else:
        '''
            -B anywhere 
            -G front 
            -A adapter
        '''
        shared_options += '--overlap 6 --nextseq-trim=25 --times=2 --max-n=3 '\
                        '--error-rate=0.1 --front={front_adapter1} --anywhere={anywhere_adapter1} '\
                        '-G {front_adapter2} -B {anywhere_adapter2} --trim-n '\
                        .format(front_adapter1 = R2, anywhere_adapter1 = rvs_byproduct,
                                    front_adapter2 = R2R, anywhere_adapter2 = fwd_byproduct)  +\
                        ' -A T{100} -A A{100} -a A{100} -a T{100} '

    if config['umi'] == 0:
        command = 'atropos trim {option} {adaptors} {shared_options} '\
                    '-o {trimed1} -p {trimed2} -pe1 {file1} -pe2 {file2}'\
                    .format(option=option, adaptors=paired_end_adaptor, shared_options=shared_options,
                            trimed1 = output['FQ1'], trimed2 = output['FQ2'],
                            file1 = input['FQ1'], file2 = input['FQ2'])

    elif config['umi'] > 0:
        command = 'clip_fastq.py --fastq1={file1} --fastq2={file2} --idxBase={umi} '\
                    ' --barcodeCutOff=20 --out_file={TEMP} -r read1 --min_length 15 ' \
                '; atropos trim {option} {shared_options} {adaptors}  '\
                '--interleaved-input {TEMP} '\
                ' --quiet  --report-file /dev/stderr -f fastq '\
                '--interleaved-out /dev/stdout '\
                '| deinterleave_fastq.py -1 {trimed1} -2 {trimed2} --min_length 15 '\
                '; rm {TEMP}'\
                .format(file1= input['FQ1'], 
                        file2= input['FQ2'], 
                        umi = config['umi']*'X',
                        option = option,
                        adaptors = paired_end_adaptor, 
                        shared_options = shared_options,
                        trimed1 = output['FQ1'], 
                        trimed2 = output['FQ2'],
                        TEMP = params['TEMP_FQ'])
    return command


def fastp_trimming(config, input, output, params):

    option = ''
    if config['TTN']:
        option += ' --trim_front2 1 '

    threads = 16 if config['threads'] > 16 else config['threads']
    shared_option = '--overrepresentation_analysis --length_required  15 '\
                '--trim_poly_x --poly_x_min_len 8 '\
                '{option} 1 --low_complexity_filter --thread {threads} '\
                '--out1 {trimed1} --out2 {trimed2} '\
                '--complexity_threshold 30 '\
                '--html {prefix}.html --json {prefix}.json '\
                .format(option = option,
                        threads = threads,
                        trimed1 = output['FQ1'],
                        trimed2 = output['FQ2'],
                        prefix = output['FQ1'].split('.')[0])
                
        
    if config['umi'] > 0:
        command = 'clip_fastq.py --fastq1={file1} --fastq2={file2} --idxBase={umi} '\
                ' --barcodeCutOff=20 --out_file=- -r read1 --min_length 15 ' \
                '| fastp --stdin --interleaved_in '\
                '{shared_option} '\
                .format(file1= input['FQ1'], 
                        file2= input['FQ2'], 
                        umi = config['umi']*'X',
                        shared_option = shared_option)

    else:
        command = ' fastp --in1 {file2} --in2 {file2} {shared_option} '\
                .format(file1= input['FQ1'], 
                        file2= input['FQ2'], 
                        shared_option = shared_option)
    
    return command
