from __future__ import division, print_function
import os
import sys
import time
from collections import defaultdict, deque
from functools import partial
from builtins import map, range
import re
import six
import pandas as pd
from tgirt_map.trim_function import fastp_trimming, atropos_trimming

class sample_object():
    def __init__(self, args):
        # read input
        self.fastq1 = args.fastq1
        self.fastq2 = args.fastq2
        self.outpath = args.outdir
        self.hisat_index = args.hisat_index
        self.bowtie2_index = args.bowtie2_index
        self.rRNA_mt_index = args.rRNA_mt_index
        self.rmsk_index = args.repeats_index
        self.univec_index = args.univec_index
        self.smRNA_index = args.smRNA_index
        self.bedpath = args.bedpath
        self.splicesite = args.splicesite
        self.threads = args.threads
        self.rmsk = args.repeats
        self.UMI = args.umi
        self.TTN = args.TTN
        self.trim_hard = args.trim_aggressive
        self.dry = args.dry
        self.count_all = args.count_all
        self.novel_splice = args.novel_splice
        self.polyA = args.polyA
        self.multi = args.multi
        self.use_fastp = args.fastp

        #### make folder
        self.trim_folder = self.outpath + '/Trim'
        self.count_folder= self.outpath + '/Counts'
        self.count_raw = self.count_folder + '/RAW'
        self.count_rmsk = self.count_folder + '/repeat_RAW'

        #define sample folder
        self.samplename = args.samplename
        self.sample_folder = self.outpath + '/' + self.samplename
        self.hisat_out = self.sample_folder + '/Hisat'
        self.rRNA_mt_out = self.sample_folder + '/rRNA_mt'
        self.bowtie_out = self.sample_folder + '/Bowtie'
        self.combined_out = self.sample_folder + '/Combined'
        self.repeat_out = self.sample_folder + '/repeats'
        self.univec_contaminants = self.sample_folder + '/UniVec'
        self.smRNA_out = self.sample_folder + '/smallRNA'

        #make output file names
        self.trimed1= '%s/%s.1.fq.gz' %(self.trim_folder, self.samplename)
        self.trimed2= self.trimed1.replace('.1.fq.gz','.2.fq.gz')

        if self.UMI == 0 or self.count_all:
                self.count_bam = self.combined_out + '/primary.bam'
        else:
                self.count_bam = self.combined_out + '/primary.dedup.bam'

        self.run_process = partial(system_run, args.dry, self.samplename)


        self.HISAT2 = 'hisat2 '\
                '--no-mixed --no-discordant '\
                ' --new-summary --dta --mp 4,2 '\
                '-p {threads} '.format(threads = self.threads) 
        
        self.BOWTIE2 = ' bowtie2  ' \
                '--very-sensitive-local ' \
                '-L 8  --mp 4,2 -N 1 '\
                '--no-mixed --no-discordant --dovetail '\
                '-p {threads}'.format(threads = self.threads)


    def make_result_dir(self):
        print('Checking output folders', file=sys.stderr)
        folders = [self.outpath, self.trim_folder, self.count_folder, self.count_raw,
                         self.sample_folder, self.hisat_out, self.rRNA_mt_out,
                        self.bowtie_out, self.combined_out, self.smRNA_out,
                        self.univec_contaminants]
        mf = deque(map(makeFolder, folders))
        if self.rmsk:
            makeFolder(self.repeat_out)
            makeFolder(self.count_rmsk)


    def trimming(self):
        trimming = fastp_trimming if self.use_fastp else atropos_trimming
        config, input, output, params = {}, {}, {}, {}
        config['TTN']  = self.TTN
        config['threads'] = self.threads
        config['umi'] = self.UMI
        config['trim_aggressive'] = self.trim_hard
        config['polyA'] = self.polyA

        input['FQ1'] = self.fastq1
        input['FQ2'] = self.fastq2

        output['FQ1'] = self.trimed1
        output['FQ2'] = self.trimed2
        command = trimming(config, input, output)
    
        self.run_process(command)


    def RNA_filter(self, RNA='univec'):

        if RNA=="univec":
            RNA_filter_out = self.univec_contaminants
            index = self.univec_index
        elif RNA == 'rRNA_mt':
            RNA_filter_out = self.rRNA_mt_out
            index = self.rRNA_mt_index
        elif RNA == 'smallRNA':
            RNA_filter_out = self.smRNA_out
            index = self.smRNA_index

        self.filtered_fq1 = RNA_filter_out + '/filtered.1.fq.gz'
        self.filtered_fq2 = RNA_filter_out + '/filtered.2.fq.gz'
        _input = '-1 {trimmed1} -2 {trimmed2}'.format(trimmed1 = self.trimed1, trimmed2 = self.trimed2)
        _out_bam = RNA_filter_out + '/aligned.bam'
        _out_bed = RNA_filter_out + '/aligned.bed'
        _out_count = RNA_filter_out + '/aligned.count'
        command = self.BOWTIE2 + \
                ' -k 1 -x {index} {input} '\
                '| samtools view -bS@{threads} - '\
                '> {out_bam} ' \
                '; samtools fastq -nf4 -1 {filtered_fq1} -2 {filtered_fq2} {out_bam}'\
                '; cat {out_bam} '\
                '| samtools view -bF2048 -F256 -F4 '\
                '| bam_to_bed.py -i - -o {out_bed} '\
                '-m 5 -M 10000 -p ' \
                .format(filtered_fq1 = self.filtered_fq1, 
                        filtered_fq2 = self.filtered_fq2,
                        index = index,
                        threads = self.threads,
                        input = _input,
                        out_bam = _out_bam,
                        out_bed = _out_bed)

        self.run_process(command)
        if not self.dry:
            count_bed(_out_bed, _out_count)


    def hisat_map(self):
        _input = '-1 {PREMAP_FASTQ1} -2 {PREMAP_FASTQ2} '\
                    .format(PREMAP_FASTQ1=self.trimed1, PREMAP_FASTQ2=self.trimed2)
        _split_option = ' '
        _unaligned = '| bamToFastq -i - -fq {bowtie_out}/unmapped.1.fq -fq2 {bowtie_out}/unmapped.2.fq'\
                    .format(bowtie_out=self.bowtie_out)
        _zip_command = 'gzip -f {bowtie_out}/unmapped.1.fq;gzip -f {bowtie_out}/unmapped.2.fq'\
                    .format(bowtie_out=self.bowtie_out)

        splice_option = ' '
        if self.novel_splice:
            splice_option = '--pen-canintronlen C,0,0 --pen-noncanintronlen C,1,0 ' \
                            '--pen-cansplice 0 --pen-noncansplice 2 --max-intronlen 1000000 '\
                            '--rna-strandness FR ' 

        # map reads
        command = self.HISAT2 +\
                ' -k {multi} --known-splicesite-infile {splicesite} {splice_option} '\
                '--novel-splicesite-outfile {hisat_out}/novelsite.txt -x {ref} {input}'\
                '| samtools view -bS -@ {threads} - > {hisat_out}/hisat.bam'\
                .format(threads=self.threads,
                        multi = self.multi,
                        splicesite=self.splicesite, 
                        splice_option = splice_option,
                        hisat_out=self.hisat_out, 
                        ref=self.hisat_index,
                        input=_input)
        self.run_process(command)

        #split to uniq and multimap
        uniq_command = 'split_uniq_bam.py -i {hisat_out}/hisat.bam '\
                '-o {hisat_out}/hisat -a hisat2 {option}'\
                .format(hisat_out=self.hisat_out, option=_split_option)
        self.run_process(uniq_command)

        #extract unaligned
        command = 'samtools view -@ {threads} -bf4 {hisat_out}/hisat.bam {unaligned_op}'\
                    .format(threads=self.threads, 
                            hisat_out=self.hisat_out, 
                            unaligned_op = _unaligned)
        self.run_process(command)
        self.run_process(_zip_command)

    def bowtie_map(self):

        _input = '-1 {fq_path}/unmapped.1.fq.gz '\
                '-2 {fq_path}/unmapped.2.fq.gz'\
                .format(fq_path = self.bowtie_out)
        _split_option = ' '

        # map reads
        command= self.BOWTIE2 + ' -k {multi} -x {index} {input} '\
                '| samtools view -@{threads} -bS - > {bowtie_out}/bowtie2.bam'\
                .format(threads=self.threads,
                        multi = self.multi, 
                        index=self.bowtie2_index,
                        input=_input,
                        bowtie_out=self.bowtie_out)
        self.run_process(command)

        # split to uniq and multimap
        command = 'split_uniq_bam.py '\
                '-i {bowtie_out}/bowtie2.bam '\
                '-o {bowtie_out}/bowtie '\
                '-a bowtie2 {option}'\
                .format(bowtie_out=self.bowtie_out, option=_split_option)
        self.run_process(command)

    def combined_aligned(self):
        command = 'samtools cat {hisat}/hisat.multi.bam '\
                        '{bowtie}/bowtie.multi.bam '\
                        ' > {out}/multi.bam' \
                .format(hisat = self.hisat_out, 
                        bowtie = self.bowtie_out,
                        out = self.combined_out)
        self.run_process(command)

        command = 'reduce_multi_reads.py --infile {combined_out}/multi.bam '\
                        '--outfile {combined_out}/multi_filtered.bam '\
                        ' --bam_in --bam_out '\
                   .format(combined_out = self.combined_out)
        self.run_process(command)

        command ='samtools cat {combined}/multi_filtered.bam {hisat}/hisat.unique.bam {bowtie}/bowtie.unique.bam' \
                '| filter_soft_clip.py -s 0.1 -b 0.2 -i - -o - --pe ' \
                '| samtools sort -n -@ {threads} -O bam -T {combined}/temp '\
                '> {combined}/primary.bam' \
                .format(combined = self.combined_out, 
                        hisat = self.hisat_out, 
                        bowtie = self.bowtie_out, 
                        threads = self.threads) 
        self.run_process(command)


    def dedup_bam(self):
        bam_file = self.combined_out + '/primary.bam'
        tag_bam = bam_file.replace('.bam','.tag.bam')
        sorted_bam = bam_file.replace('.bam','.sorted.bam')
        dedup_bam = bam_file.replace('.bam','.dedup.bam')
        umi_text = self.combined_out + '/primary.umi_metrics'
        duplicate_text = self.combined_out + '/primary.duplicate_metrics'

        umi_command = 'bam_umi_tag.py --in_bam {inbam} --out_bam - --tag RX '\
                        '| picard SortSam I=/dev/stdin O=/dev/stdout SORT_ORDER=queryname '\
                        '| picard FixMateInformation ADD_MATE_CIGAR=true ASSUME_SORTED=true '\
                            'INPUT=/dev/stdin OUTPUT=/dev/stdout '\
                        '> {outbam}'\
                        .format(inbam = bam_file, outbam = tag_bam)

        sort_command = 'samtools sort -@ %i -O bam -T %s/temp %s > %s' \
                            %(self.threads, self.combined_out, tag_bam, sorted_bam)
        index_command = 'samtools index %s' %(sorted_bam)
        dedup_command = 'picard UmiAwareMarkDuplicatesWithMateCigar UMI_METRICS_FILE=%s ' %(umi_text)+\
                        'MAX_EDIT_DISTANCE_TO_JOIN=1 TAG_DUPLICATE_SET_MEMBERS=true' +\
                        'UMI_TAG_NAME=RX INPUT=%s OUTPUT=%s ' %(sorted_bam, tag_bam) +\
                        'METRICS_FILE=%s REMOVE_DUPLICATES=false ASSUME_SORT_ORDER=coordinate' %(duplicate_text)
        resort_command = 'samtools sort -n@  %i -O bam -T %s/temp %s > %s ' \
                            %(self.threads, self.combined_out, tag_bam, dedup_bam)

        self.run_process(umi_command)
        self.run_process(sort_command)
        self.run_process(index_command)
        self.run_process(dedup_command)
        self.run_process(resort_command)

    def combined_filter(self):

        _verb = 'pairtobed'
        _option = '-type neither'


        ### filter out tRNA
        command = 'bedtools {verb} -s -f 0.01 -abam {combined_path}/primary.bam'\
                    ' -b {bed_path}/tRNA.bed '\
                    '> {out}/tRNA_primary.bam' \
                    .format(combined_path = self.combined_out, 
                            bed_path = self.bedpath, 
                            out=self.combined_out, 
                            verb = _verb)
        self.run_process(command)

        ### filter out rRNA
        command = 'bedtools {verb} -s -f 0.01 '\
                    '-abam {combined_path}/primary.bam '\
                    ' -b {bed_path}/rRNA_for_bam_filter.bed '\
                    '> {out}/rRNA_primary.bam' \
                    .format(combined_path = self.combined_out, 
                        bed_path = self.bedpath, 
                        out=self.combined_out, 
                        verb = _verb)
        self.run_process(command)

        ### filter out smallRNA
        command = 'bedtools {verb} -s -f 0.01 -abam {count_bam} '\
                        '-b {bed_path}/sncRNA_no_tRNA.bed '\
                        ' > {combined_path}/sncRNA.bam' \
                        .format(count_bam = self.count_bam, 
                                bed_path = self.bedpath, 
                                combined_path = self.combined_out, 
                                verb = _verb)
        self.run_process(command)

        ### filter out long RNA
        command = 'bedtools {verb} -s -f 0.01 {option} '\
                    '-abam {count_bam} '\
                    '-b {bed_path}/sncRNA_rRNA_for_bam_filter.bed '\
                    '> {combined_path}/primary_no_sncRNA_tRNA_rRNA.bam' \
                        .format(count_bam = self.count_bam, 
                                bed_path = self.bedpath, 
                                combined_path = self.combined_out, 
                                option=_option, 
                                verb = _verb)
        self.run_process(command)

        if self.rmsk:
            command = 'bedtools {verb} -f 0.01 {option} '\
                    '-abam {combined_path}/primary_no_sncRNA_tRNA_rRNA.bam '\
                    '-b {rmsk_bed} '\
                    '> {combined_path}/primary_no_sncRNA_tRNA_rRNA_repeats.bam' \
                .format(count_bam = self.count_bam, 
                        option=_option, 
                        bed_path = self.bedpath, 
                        combined_path = self.combined_out, 
                        verb = _verb, 
                        rmsk_bed = self.rmsk)
            self.run_process(command)

            command = 'bedtools {verb} -f 0.5 '\
                    '-abam {combined_path}/primary_no_sncRNA_tRNA_rRNA.bam '\
                    '-b {rmsk_bed} -type both '\
                    '> {combined_path}/repeats.bam' \
                .format(count_bam = self.count_bam, 
                        bed_path = self.bedpath, 
                        combined_path = self.combined_out, 
                        verb = _verb, 
                        rmsk_bed = self.rmsk)
            self.run_process(command)


    def make_alignment_bed(self):

        command = 'bam_to_bed.py -m 5 -M 1000000  '\
                '-i {combined_out}/sncRNA.bam > {combined_out}/sncRNA.bed'\
                .format(combined_out=self.combined_out)
        self.run_process(command)

        command = 'bam_to_bed.py '\
                '-i {combined_out}/primary_no_sncRNA_tRNA_rRNA.bam  '\
                '-m 5 -M 1000000 '\
                '> {combined_out}/primary_no_sRNAs.bed'\
                .format(combined_out=self.combined_out)
        self.run_process(command)



    def generate_all_count(self):
        command = 'bedtools coverage -s -counts -F 0.1 '\
                    '-a {bed_path}/sncRNA_no_tRNA.bed '\
                    '-b {combined}/sncRNA.bed '\
                    '> {combined}/sncRNA.counts'\
                    .format(combined=self.combined_out,
                            bed_path=self.bedpath)
        self.run_process(command)
        command = 'bedtools coverage -s -counts -F 0.1 '\
                    '-a {bed_path}/genes_no_sncRNA_rRNA_tRNA.bed '\
                    '-b {combined}/primary_no_sRNAs.bed '\
                    '> {combined}/non_sRNAs.counts'\
                    .format(combined=self.combined_out,
                            bed_path=self.bedpath)
        self.run_process(command)

        command = 'cat {combined}/non_sRNAs.counts '\
                    '{combined}/sncRNA.counts '\
                '> {count_path}/{samplename}.counts'\
                .format(count_path=self.count_raw, 
                        samplename = self.samplename,
                        combined=self.combined_out) 
        self.run_process(command)

    def generate_repeat_count(self):
        command = 'samtools fastq -N@ {threads} {combined_path}/repeats.bam '\
                '-1 {repeat_path}/repeats_1.fq.gz -2 {repeat_path}/repeats_2.fq.gz'\
                    .format(repeat_path=self.repeat_out,
                            combined_path =self.combined_out, 
                            threads = self.threads)
        self.run_process(command)
        fq_input = ' -1 {repeat_path}/repeats_1.fq.gz -2 {repeat_path}/repeats_2.fq.gz '\
                    .format(repeat_path = self.repeat_out)

        command = self.BOWTIE2 + \
                ' -x {repeat_index} {input} '\
                '| samtools view -bS@ {threads} - > {repeat_path}/repeat_remap.bam'\
                .format(threads=self.threads,
                        repeat_index=self.rmsk_index,
                        input = fq_input, 
                        repeat_path=self.repeat_out)
        command += '; filter_umi.py -i {repeat_path}/repeat_remap.bam '\
                '--consecutive_bases 3 '\
                '| bam_to_bed.py -i - -o {repeat_path}/repeat.bed '\
                '-m 5 -M 10000'\
                .format(repeat_path=self.repeat_out)
        self.run_process(command)

        repeat_count = defaultdict(lambda: defaultdict(int))
        repeat_bed = self.repeat_out + '/repeat.bed'
        repeat_count_file = self.count_rmsk + '/' + self.samplename + '.repeat'
        print('Reading from %s' %repeat_bed, file=sys.stderr)
        if not self.dry:
            with open(repeat_bed,'r') as bed:
                for line in bed:
                    fields = line.strip().split('\t')
                    repeat = fields[0]
                    strand = fields[5]
                    repeat_count[repeat][strand] += 1

            with open(repeat_count_file, 'w') as count_file:
                for repeat_name, strand_dict in six.iteritems(repeat_count):
                    for strand, value in six.iteritems(strand_dict):
                        print('%s\t%s\t%i' %(repeat_name, strand, value), file=count_file)
        print('Written %s' %repeat_count_file, file=sys.stderr)


def count_rRNA(RNA, start, end):
    gene = 'rDNA'
    RNA_5S = RNA == 'gi|23898|emb|X12811.1|' and end > 274  and start < 394
    RNA_18S = RNA== 'gi|555853|gb|U13369.1|HSU13369' and  end > 3657  and start < 5527
    RNA_58S = RNA == 'gi|555853|gb|U13369.1|HSU13369' and end > 6623  and start <  6779
    RNA_28S = RNA == 'gi|555853|gb|U13369.1|HSU13369' and end > 7935 and start < 12969
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
    count_dict = defaultdict(int)
    
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
    
    pd.DataFrame({'gene':list(count_dict.keys()),
                       'count': list(count_dict.values())}) \
        .filter(['gene','count'])\
        .to_csv(out_count, index=False, sep='\t', header=False)
    print('Written %s\n' %out_count, file=sys.stderr)


def system_run(dry, samplename, command):
    print('[%s] Running: %s' %(samplename, command), file=sys.stderr)
    if dry:
        return 0

    else:
        start = time.time()
        os.system(command)
        end = time.time() - start
        print('[%s] Used time %.3f min\n' %(samplename, end/60), file=sys.stderr)
        return 0


def makeFolder(folder):
    """
            Input a folder name and make a folder if it is non-existed
    """
    print('Creating %s....' %folder, file = sys.stderr)
    if os.path.isdir(folder):
        print('%s exists.' %folder, file = sys.stderr)
    else:
        os.mkdir(folder)
        print('Created %s.' %folder, file = sys.stderr)
    return 0
