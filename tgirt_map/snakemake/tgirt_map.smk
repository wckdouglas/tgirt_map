from collections import defaultdict, Counter
from pandas import DataFrame
from tgirt_map.trim_function import atropos_trimming, fastp_trimming

# set up config
FASTQ1 = config["fastq1"]
FASTQ2 = config["fastq2"]
OUTPATH = config["outdir"]
HISAT_INDEX = config["hisat_index"]
BOWTIE2_INDEX = config["bowtie2_index"]
RRNA_MT_INDEX = config["rRNA_mt_index"]
RMSK_INDEX = config["repeats_index"]
UNIVEC_INDEX = config["univec_index"]
SMRNA_INDEX = config["smRNA_index"]
BEDPATH = config["bedpath"]
SPLICESITE = config["splicesite"]
THREADS = config["threads"]
RMSK = config["repeats"]
UMI = config["umi"]
TTN = config["TTN"]
TRIM_HARD = config["trim_aggressive"]
COUNT_ALL = config["count_all"]
NOVEL_SPLICE = config["novel_splice"]
POLYA = config["polyA"]
MULTI = config["multi"]
SAMPLENAME = config["samplename"]

# setup path
COUNT_PATH = OUTPATH + '/Counts'
RAW_COUNT_PATH = COUNT_PATH + '/RAW'
REPEAT_COUNT_PATH = COUNT_PATH + '/repeat_RAW'

SAMPLE_FOLDER = OUTPATH + '/' + SAMPLENAME
UNIVEC_FOLDER = SAMPLE_FOLDER + '/UniVec'
mt_rRNA_FOLDER = SAMPLE_FOLDER + '/rRNA_mt'
SMALL_RNA_FOLDER = SAMPLE_FOLDER + '/smallRNA'
HISAT_FOLDER = SAMPLE_FOLDER + '/hisat'
BOWTIE_FOLDER = SAMPLE_FOLDER + '/bowtie'
COMBINED_FOLDER = SAMPLE_FOLDER + '/Combined'
REPEAT_FOLDER = SAMPLE_FOLDER + '/repeats'

#setup file
TRIMMED_FQ1 = SAMPLE_FOLDER + '/Trimmed/trim.1.fq.gz'
TRIMMED_FQ2 = SAMPLE_FOLDER + '/Trimmed/trim.2.fq.gz'
UNIVEC_MAPPED_BAM = UNIVEC_FOLDER + '/aligned.bam'
UNIVEC_MAPPED_BED = UNIVEC_FOLDER + '/aligned.bed'
UNIVEC_UNMAPPED_FQ1 = UNIVEC_FOLDER + '/filtered.1.fq.gz'
UNIVEC_UNMAPPED_FQ2 = UNIVEC_FOLDER + '/filtered.2.fq.gz'
UNIVEC_COUNT = UNIVEC_FOLDER + '/aligned.count'


mt_rRNA_MAPPED_BAM = mt_rRNA_FOLDER + '/aligned.bam'
mt_rRNA_MAPPED_BED = mt_rRNA_FOLDER + '/aligned.bed'
mt_rRNA_UNMAPPED_FQ1 = mt_rRNA_FOLDER + '/filtered.1.fq.gz'
mt_rRNA_UNMAPPED_FQ2 = mt_rRNA_FOLDER + '/filtered.2.fq.gz'
mt_rRNA_COUNT = mt_rRNA_FOLDER + '/aligned.count'

SNC_RNA_MAPPED_BAM = SMALL_RNA_FOLDER + '/aligned.bam'
SNC_RNA_MAPPED_BED = SMALL_RNA_FOLDER + '/aligned.bed'
SNC_RNA_UNMAPPED_FQ1 = SMALL_RNA_FOLDER + '/filtered.1.fq.gz'
SNC_RNA_UNMAPPED_FQ2 = SMALL_RNA_FOLDER + '/filtered.2.fq.gz'
SNC_RNA_COUNT = SMALL_RNA_FOLDER + '/aligned.count'

HISAT_BAM = HISAT_FOLDER + '/aligned.bam'
HISAT_UNIQUE_BAM = HISAT_FOLDER + '/hisat.unique.bam'
HISAT_MULTI_BAM = HISAT_FOLDER + '/hisat.multi.bam'
HISAT_UNMAPPED_FQ1 = HISAT_FOLDER + '/unmapped.1.fq.gz'
HISAT_UNMAPPED_FQ2 = HISAT_FOLDER + '/unmapped.2.fq.gz'
BOWTIE_BAM = BOWTIE_FOLDER + '/aligned.bam'
BOWTIE_UNIQUE_BAM = BOWTIE_FOLDER + '/bowtie.unique.bam'
BOWTIE_MULTI_BAM = BOWTIE_FOLDER + '/bowtie.multi.bam'
MULTI_BAM = COMBINED_FOLDER + '/multi.bam'
FILTERED_MULTI_BAM = COMBINED_FOLDER + '/multi_filtered.bam'
PRIMARY_BAM = COMBINED_FOLDER + '/primary.bam'
PRIMARY_SNC_BAM = COMBINED_FOLDER + '/sncRNA.bam'
PRIMARY_rRNA_BAM = COMBINED_FOLDER + '/rRNA_primary.bam'
PRIMARY_tRNA_BAM = COMBINED_FOLDER + '/tRNA_primary.bam'
PRIMARY_NO_SNC_BAM = COMBINED_FOLDER + '/primary_no_sncRNA_tRNA_rRNA.bam'
PRIMARY_NO_SNC_REPEAT_BAM = COMBINED_FOLDER + '/primary_no_sncRNA_tRNA_rRNA_repeats.bam'
PRIMARY_REPEAT_BAM = COMBINED_FOLDER + '/repeats.bam'
PRIMARY_SNC_BED = COMBINED_FOLDER + '/sncRNA.bed'
PRIMARY_BED = COMBINED_FOLDER + '/primary_no_sRNAs.bed'

PRIMARY_BED = COMBINED_FOLDER + '/primary_no_sRNAs.bed'
RMSK_FQ1 = REPEAT_FOLDER + '/repeats.1.fq.gz'
RMSK_FQ2 = REPEAT_FOLDER + '/repeats.2.fq.gz'
RMSK_BED = REPEAT_FOLDER + '/repeats.bed'
RMSK_BAM = REPEAT_FOLDER + '/repeats_remap.bam'

SNC_COUNT = COMBINED_FOLDER + '/sncRNA.counts'
NO_SNC_COUNT = COMBINED_FOLDER + '/non_sRNAs.counts'
ALL_COUNT_FILE = RAW_COUNT_PATH + '/' + SAMPLENAME + '.counts'
REPEAT_COUNT_FILE = REPEAT_COUNT_PATH + '/' + SAMPLENAME + '.counts'


# mapping params
HISAT2 = 'hisat2 '\
        '--no-mixed --no-discordant '\
        ' --new-summary --dta --mp 4,2 '\
        '-p {threads} '.format(threads = config["threads"]) 

BOWTIE2 = ' bowtie2  ' \
        '--very-sensitive-local ' \
        '-L 8  --mp 4,2 -N 1 '\
        '--no-mixed --no-discordant --dovetail '\
        '-p {threads}'.format(threads = config['threads'])

#shared command
RNA_FILTER_ALIGN_COMMAND = BOWTIE2 + ' -k 1 -x {params.INDEX} '\
    '-1 {input.FQ1} -2 {input.FQ2} ' \
    '| samtools view -bS@ {params.THREADS} - '\
    '> {output.BAM} '

RNA_FILTER_FQ_COMMAND = ' samtools fastq -@ {params.THREADS} -nf4 '\
    '-1 {output.FQ1} -2 {output.FQ2} {input.BAM} '

RNA_BED_COMMAND = 'cat {input.BAM} '\
    '| samtools view -bF 2048 -F256 -F4 '\
    '| bam_to_bed.py -i - -o {output.BED} -m 5 -M 10000 -p '

BAM_EXTRACT_COMMAND = 'bedtools pairtobed -s -f 0.01 '\
        '-abam {input.BAM} '\
        '-type either -b {params.BED} '\
        '> {output.BAM}'

RNA_COUNT_COMMAND = 'bedtools coverage '\
        '-a {params.ANNOTATION} -b {input.BED} '\
        '-s -counts -F 0.1 > {output.COUNT_FILE}'

BAM_TO_BED_COMMAND = 'bam_to_bed.py -m 5 -M 1000000 -i {input.BAM} -o {output.BED}'


# pipeline starts
rule all:
    input:
        ALL_COUNT_FILE,
        REPEAT_COUNT_FILE,
        SNC_RNA_COUNT,
        UNIVEC_COUNT,
        mt_rRNA_COUNT,
        PRIMARY_NO_SNC_REPEAT_BAM,

rule generate_all_count:
    input:
        sRNA_COUNT = SNC_COUNT,
        non_sRNA_COUNT = NO_SNC_COUNT

    output:
        COUNT_FILE = ALL_COUNT_FILE 

    shell:
        'cat {input.sRNA_COUNT} {input.non_sRNA_COUNT} '\
        '> {output.COUNT_FILE} '
    
rule generate_lRNA_count:
    input:
        BED = PRIMARY_BED
    
    params:
        ANNOTATION = BEDPATH + '/genes_no_sncRNA_rRNA_tRNA.bed'
    
    output:
        COUNT_FILE = NO_SNC_COUNT
    
    shell:
        RNA_COUNT_COMMAND

rule generate_snc_count:
    input:
        BED = PRIMARY_SNC_BED

    params:
        ANNOTATION = BEDPATH + '/sncRNA_no_tRNA.bed'

    output:
        COUNT_FILE = SNC_COUNT
    
    shell:
        RNA_COUNT_COMMAND
        
rule make_sncRNA_bed:
    input:
        BAM = PRIMARY_SNC_BAM
    
    output:
        BED = PRIMARY_SNC_BED
    
    shell:
        BAM_TO_BED_COMMAND

rule make_primary_bed:
    input:
        BAM = PRIMARY_NO_SNC_BAM
    
    output:
        BED = PRIMARY_BED
    
    shell:
        BAM_TO_BED_COMMAND



rule count_repeats:
    input:
        RMSK_BED
    
    output:
        REPEAT_COUNT_FILE

    run:
        repeat_count = defaultdict(Counter)
        repeat_bed = input[0]
        repeat_count_file = output[0]

        print('Reading from %s' %repeat_bed, file=sys.stderr)
        with open(repeat_bed,'r') as bed:
            for line in bed:
                fields = line.strip().split('\t')
                repeat = fields[0]
                strand = fields[5]
                repeat_count[repeat][strand] += 1

        with open(repeat_count_file, 'w') as count_file:
            for repeat_name, strand_dict in repeat_count.items():
                for strand, value in strand_dict.items():
                    print('%s\t%s\t%i' %(repeat_name, strand, value), file=count_file)


rule filter_repeats:
    input:
        PRIMARY_NO_SNC_BAM

    params:
        BED = RMSK

    output:
        PRIMARY_NO_SNC_REPEAT_BAM
    
    shell:
        'bedtools pairtobed -type neither '\
        '-abam {input} '\
        '-b {params.BED}'
        '> {output}'

rule make_repeat_bed:
    input:
        FQ1 = RMSK_FQ1,
        FQ2 = RMSK_FQ2 
    
    params:
        INDEX = RMSK_INDEX,
        THREADS = THREADS

    output:
        BED = RMSK_BED,
        BAM = RMSK_BAM

    shell:
        BOWTIE2 + \
        ' -x {params.INDEX} -1 {input.FQ1} -2 {input.FQ2} '  \
        '| samtools view -F4 -b -@ {params.THREADS} '\
        '| tee {output.BAM} '\
        '| bam_to_bed.py -i - -o {output.BED}  -m 5 -M 10000'


rule make_repeat_fq:
    input:
        BAM = PRIMARY_REPEAT_BAM

    params:
        THREADS = THREADS

    output:
        FQ1 = RMSK_FQ1, 
        FQ2 = RMSK_FQ2
    

    shell:
        'samtools fastq -N@ {params.THREADS} {input.BAM} '\
        '-1 {output.FQ1} -2 {output.FQ2}'


rule extract_repeat_bam:
    input:
        BAM = PRIMARY_NO_SNC_BAM
    
    params:
        RMSK_BED = RMSK

    output:
        BAM = PRIMARY_REPEAT_BAM

    shell:
        'bedtools pairtobed -f 0.5 '\
        '-abam {input.BAM} -b {params.RMSK_BED} -type both '\
        '> {output.BAM}'
    

rule filter_snc_bam:
    input:
        BAM = PRIMARY_BAM
    
    params:
        SNC_RNA_BED = BEDPATH + '/sncRNA_rRNA_for_bam_filter.bed'

    output:
        BAM = PRIMARY_NO_SNC_BAM
    
    shell:
        'bedtools pairtobed -s -f 0.01 -abam {input.BAM} '\
        '-b {params.SNC_RNA_BED} -type neither '\
        '> {output.BAM} '


rule extract_snc_bam:
    input:
        BAM = PRIMARY_BAM
    
    params:
        BED = BEDPATH + '/sncRNA_no_tRNA.bed' 
    
    output:
        BAM = PRIMARY_SNC_BAM

    shell:
        BAM_EXTRACT_COMMAND

rule extract_rRNA_bam:
    input:
        BAM = PRIMARY_BAM
    
    params:
        BED = BEDPATH + '/rRNA_for_bam_filter.bed'
    
    output:
        BAM = PRIMARY_rRNA_BAM

    shell: 
        BAM_EXTRACT_COMMAND


rule extract_tRNA_bam:
    input:
        BAM = PRIMARY_BAM
    
    params:
        BED = BEDPATH + '/tRNA.bed'
    
    output:
        BAM = PRIMARY_tRNA_BAM
    
    shell:
        BAM_EXTRACT_COMMAND
    

rule make_primary_bam:
    input:
        MULTI = FILTERED_MULTI_BAM,
        HISAT = HISAT_UNIQUE_BAM,
        BOWTIE = BOWTIE_UNIQUE_BAM

    params:
        THREADS = THREADS,
        TEMP = COMBINED_FOLDER

    output:
        BAM = PRIMARY_BAM

    shell:
        'samtools cat {input.MULTI} {input.HISAT} {input.BOWTIE} '\
        '| filter_soft_clip.py -s 0.1 -b 0.2 -i - -o - --pe '\
        '| samtools sort -n@ {params.THREADS} -O bam -T {params.TEMP} -o {output.BAM}'


rule filter_multi:
    input:
        HISAT = HISAT_MULTI_BAM,
        BOWTIE = BOWTIE_MULTI_BAM

    output:
        BAM = MULTI_BAM,
        FILTERED_BAM = FILTERED_MULTI_BAM

    shell:
        'samtools cat {input.HISAT} {input.BOWTIE} ' \
        '| tee {output.BAM} '\
        '| reduce_multi_reads.py --infile - --outfile {output.FILTERED_BAM} '\
        '--bam_in --bam_out'

rule make_bowtie_multi:
    input:
        BAM = BOWTIE_BAM
    
    params:
        OUT_PREFIX = BOWTIE_FOLDER + '/bowtie'

    output:
        UNIQUE = BOWTIE_UNIQUE_BAM,
        MULTI = BOWTIE_MULTI_BAM
    
    shell:
        'split_uniq_bam.py -i {input.BAM} -o {params.OUT_PREFIX} -a bowtie2'


rule make_hisat_multi:
    input:
        BAM = HISAT_BAM
    
    params:
        OUT_PREFIX = HISAT_FOLDER + '/hisat'

    output:
        UNIQUE = HISAT_UNIQUE_BAM,
        MULTI = HISAT_MULTI_BAM
    
    shell:
        'split_uniq_bam.py -i {input.BAM} -o {params.OUT_PREFIX} -a bowtie2'


rule bowtie_align:
    input:
        FQ1 = HISAT_UNMAPPED_FQ1,
        FQ2 = HISAT_UNMAPPED_FQ2

    params:
        THREADS = THREADS,
        MULTI = MULTI,
        INDEX = BOWTIE2_INDEX

    output:
        BAM = BOWTIE_BAM

    shell:
        BOWTIE2 + ' -k {params.MULTI} -x {params.INDEX} '\
        '-1 {input.FQ1} -2 {input.FQ2} '\
        '| samtools view -@ {params.THREADS} -bS - '\
        '> {output.BAM}'

rule hisat_unmapped:
    input:
        BAM = HISAT_BAM

    params:
        THREADS = THREADS    

    output:
        FQ1 = HISAT_UNMAPPED_FQ1,
        FQ2 = HISAT_UNMAPPED_FQ2

    shell:
        'samtools fastq -f4 -@ {params.THREADS} '\
        ' -1 {output.FQ1} -2 {output.FQ2} {input.BAM}'

    
rule hisat_align:
    input:
        FQ1 = SNC_RNA_UNMAPPED_FQ1,
        FQ2 = SNC_RNA_UNMAPPED_FQ2
    
    params:
        THREADS = THREADS,
        MULTI = MULTI,
        INDEX = HISAT_INDEX,
        SPLICESITE = SPLICESITE,

    output:
        BAM = HISAT_BAM

    shell:
        HISAT2 + ' -k {params.MULTI} --known-splicesite-infile {params.SPLICESITE}' \
        '-x {params.INDEX} -1 {input.FQ1} -2 {input.FQ2} ' \
        '| samtools view -bS@ {params.THREADS} - '\
        '> {output.BAM}'

#### SMALL RNA FILTER ####
rule count_smRNA:
    input:
        SNC_RNA_MAPPED_BED

    output:
        SNC_RNA_COUNT

    run:
        count_bed(input[0], output[0])



rule make_smallRNA_BED:
    input:
        BAM = SNC_RNA_MAPPED_BAM

    output:
        BED = SNC_RNA_MAPPED_BED

    shell:
        RNA_BED_COMMAND

rule make_unalign_smallRNA_fq:
    input:
        BAM = SNC_RNA_MAPPED_BAM

    params:
        THREADS = THREADS,

    output:
        FQ1 = SNC_RNA_UNMAPPED_FQ1,
        FQ2 = SNC_RNA_UNMAPPED_FQ2

    shell:
        RNA_FILTER_FQ_COMMAND


rule smallRNA_align:
    input:
        FQ1 = mt_rRNA_UNMAPPED_FQ1,
        FQ2 = mt_rRNA_UNMAPPED_FQ2,
    
    params:
        THREADS = THREADS,
        INDEX = SMRNA_INDEX

    output:
        BAM = SNC_RNA_MAPPED_BAM

    shell:
        RNA_FILTER_ALIGN_COMMAND



#### MT rRNA FILTER ####
rule count_rRNA_align:
    input:
        BED = mt_rRNA_MAPPED_BED

    params:
        GENE_MODEL = BEDPATH + '/rRNA_mt.bed'

    output:
        COUNT = mt_rRNA_COUNT

    run:
        'bedtools coverage -a {params.GENE_MODEL} -b {input.bed} '\
        '-counts -s > {output.COUNT}'



rule make_rRNA_BED:
    input:
        BAM = mt_rRNA_MAPPED_BAM

    output:
        BED = mt_rRNA_MAPPED_BED

    shell:
        RNA_BED_COMMAND

rule make_unalign_rRNA_fq:
    input:
        BAM = mt_rRNA_MAPPED_BAM

    params:
        THREADS = THREADS

    output:
        FQ1 = mt_rRNA_UNMAPPED_FQ1,
        FQ2 = mt_rRNA_UNMAPPED_FQ2

    shell:
        RNA_FILTER_FQ_COMMAND


rule rRNA_align:
    input:
        FQ1 = UNIVEC_UNMAPPED_FQ1,
        FQ2 = UNIVEC_UNMAPPED_FQ2,
    
    params:
        THREADS = THREADS,
        INDEX = RRNA_MT_INDEX

    output:
        BAM = mt_rRNA_MAPPED_BAM

    shell:
        RNA_FILTER_ALIGN_COMMAND



#### UNIVEC FILTER ####
rule count_univec:
    input:
        UNIVEC_MAPPED_BED

    output:
        UNIVEC_COUNT

    run:
        count_bed(input[0], output[0])

rule make_univec_BED:
    input:
        BAM = UNIVEC_MAPPED_BAM

    output:
        BED = UNIVEC_MAPPED_BED

    shell:
        RNA_BED_COMMAND

rule make_unalign_univec_fq:
    input:
        BAM = UNIVEC_MAPPED_BAM

    params:
        THREADS = THREADS

    output:
        FQ1 = UNIVEC_UNMAPPED_FQ1,
        FQ2 = UNIVEC_UNMAPPED_FQ2

    shell:
        RNA_FILTER_FQ_COMMAND


rule univec_align:
    input:
        FQ1 = TRIMMED_FQ1,
        FQ2 = TRIMMED_FQ2,
    
    params:
        THREADS = THREADS,
        INDEX = UNIVEC_INDEX

    output:
        BAM = UNIVEC_MAPPED_BAM

    shell:
        RNA_FILTER_ALIGN_COMMAND


### TRIMMING ###
rule trim:
    input:
        FQ1 = FASTQ1,
        FQ2 = FASTQ2
    
    output:
        FQ1 = TRIMMED_FQ1,
        FQ2 = TRIMMED_FQ2

    run:
        trimming = fastp_trimming if config['fastp'] else atropos_trimming
        command = trimming(config, input, output)
        print(command)
        os.system(command)


def count_rRNA(RNA, start, end):
    gene = 'rDNA'
    RNA_5S = RNA == 'gi|23898|emb|X12811.1|' and end >= 274  and start <= 394
    RNA_18S = RNA== 'gi|555853|gb|U13369.1|HSU13369' and  end >= 3657  and start <= 5527
    RNA_58S = RNA == 'gi|555853|gb|U13369.1|HSU13369' and end >= 6623  and start <=  6779
    RNA_28S = RNA == 'gi|555853|gb|U13369.1|HSU13369' and end >= 7935 and start <= 12969
    if RNA_5S:
        gene = '5S_rRNA'
    elif RNA_18S:
        gene = '18S_rRNA'
    elif RNA_58S:
        gene = '5.8S_rRNA'
    elif RNA_28S:
        gene = '28S_rRNA'
    return gene


def count_bed(inbed, out_count):
    count_dict = Counter()

    with open(inbed, 'r') as inb:
        for line in inb:
            fields = line.rstrip().split('\t')
            RNA = fields[0]
            strand = fields[5]
            if strand == '+':
                if not RNA.startswith('gi'):
                    count_dict[RNA] += 1
                else:
                    gene = count_rRNA(RNA, int(fields[1]), int(fields[2]))
                    count_dict[gene] += 1

    DataFrame({'gene':list(count_dict.keys()),
                       'count': list(count_dict.values())}) \
        .filter(['gene','count'])\
        .to_csv(out_count, index=False, sep='\t', header=False)
    print('Written %s\n' %out_count, file=sys.stderr)

