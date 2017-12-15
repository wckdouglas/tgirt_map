from __future__ import division, print_function
import os
import sys
import time
from collections import defaultdict
from functools import partial
import re

class sample_object():
    def __init__(self, args):
        # read input
        self.fastq1 = args.fastq1
        self.fastq2 = args.fastq2
        self.outpath = args.outdir
        self.hisat_index = args.hisat_index
        self.bowtie2_index = args.bowtie2_index
        self.bedpath = args.bedpath
        self.splicesite = args.splicesite
        self.tRNA_index = args.tRNAindex
        self.rRNA_index = args.rRNAindex
        self.rRNA_tRNA_index = args.rRNA_tRNA_index
        self.threads = args.threads
        self.rmsk = args.repeats
        self.rmsk_index = args.repeats_index
        self.UMI = args.umi
        self.TTN = args.TTN
        self.dry = args.dry
        self.count_all = args.count_all

        self.single_end = not self.fastq2
        if self.single_end and self.UMI > 0:
            sys.exit('Cant do UMI for single end!!!')

        #### make folder
        self.trim_folder = self.outpath + '/Trim'
        self.count_folder= self.outpath + '/Counts'
        self.count_raw = self.count_folder + '/RAW'
        self.tRNA_raw = self.count_folder + '/tRNA_RAW'
        self.count_rmsk = self.count_folder + '/repeat_RAW'

        #define sample folder
        self.samplename = re.sub('.fastq|.gz|.fq','',self.fastq1.split('/')[-1])
        self.sample_folder = self.outpath + '/' + self.samplename
        self.hisat_out = self.sample_folder + '/Hisat'
        self.rRNA_tRNA_out = self.sample_folder + '/rRNA_tRNA_premap'
        self.bowtie_out = self.sample_folder + '/Bowtie'
        self.combined_out = self.sample_folder + '/Combined'
        self.tRNA_out = self.sample_folder + '/tRNA'
        self.rRNA_out = self.sample_folder + '/rRNA'
        self.repeat_out = self.sample_folder + '/repeats'

        #make output file names
        self.trimed1= '%s/%s.1.fq.gz' %(self.trim_folder, self.samplename)
        self.trimed2= self.trimed1.replace('.1.fq.gz','.2.fq.gz')

        self.premap_fastq1 = '%s/non_tRNA_rRNA.1.fq' %self.rRNA_tRNA_out
        self.premap_fastq2 = '%s/non_tRNA_rRNA.2.fq' %self.rRNA_tRNA_out

        self.tRNA_fastq1 = '%s/tRNA.1.fq' %self.rRNA_tRNA_out
        self.tRNA_fastq2 = '%s/tRNA.2.fq' %self.rRNA_tRNA_out

        self.rRNA_fastq1 = '%s/rRNA.1.fq' %self.rRNA_tRNA_out
        self.rRNA_fastq2 = '%s/rRNA.2.fq' %self.rRNA_tRNA_out

        if self.UMI == 0 or self.count_all:
                self.count_bam = self.combined_out + '/primary.bam'
        else:
                self.count_bam = self.combined_out + '/primary.dedup.bam'

        self.run_process = partial(system_run, args.dry, self.samplename)

    def make_result_dir(self):
        folders = [self.outpath, self.trim_folder, self.count_folder, self.count_raw,
                         self.tRNA_raw, self.sample_folder, self.hisat_out, self.rRNA_tRNA_out,
                        self.bowtie_out, self.combined_out, self.tRNA_out, self.rRNA_out]
        mf = map(makeFolder, folders)
        if self.rmsk:
            makeFolder(self.repeat_out)
            makeFolder(self.count_rmsk)


    def trimming(self):
        R2R = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
        R1R = 'GATCGTCGGACTGTAGAACTCTGAACGTGTAGA'
        R2 = 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTN'
        if self.TTN:
                option='-U 1'
                R2R = 'A' + R2R 
                R2 = R2.replace('TCTN','TCTTN')
        else:
                option = ''
        R2_frac = R2[-14:]
        R2R_frac = R2R[:14]


        single_end_adaptor = '-a {R2R} -b {R2_frac} -g {R2} '.format(R2R=R2R, R2 = R2, R2_frac=R2_frac)
        paired_end_adaptor = single_end_adaptor + '-A {R1R} -B {R2R_frac} -G {R2R}'.format(R2R=R2R, R1R=R1R, R2R_frac=R2R_frac)
        shared_options = '--nextseq-trim=20 --error-rate=0.2  -q 20 -m 15 -O 5 -n 3 '

        if self.UMI == 0:
            if not self.single_end:
                command = 'cutadapt {option} {adaptors} {shared_options} -o {trimed1} -p {trimed2} {file1} {file2}'\
                        .format(option=option, adaptors=paired_end_adaptor, shared_options=shared_options,
                                trimed1=self.trimed1, trimed2=self.trimed2,
                                file1= self.fastq1, file2= self.fastq2)
            else:
                command = 'cutadapt {option} {shared_options} --nextseq-trim=20 --error-rate=0.2 -q 20 -o {trimed1} {file1}'\
                        .format(adaptors=single_end_adaptor, option= option, shared_options=shared_options,
                                trimed1=self.trimed1,
                                file1= self.fastq1)
        elif self.UMI > 0:
            command = 'clip_fastq.py --fastq1={file1} --fastq2={file2} --idxBase={umi} '.format(file1= self.fastq1, file2= self.fastq2, umi=self.UMI*'X')+\
                        '--barcodeCutOff=20 --outputprefix=- --prefix_split=0 -r read1 '+\
                    '| cutadapt {option} {shared_options} --interleaved --quiet - '.format(option=option, adaptors=paired_end_adaptor, shared_options=shared_options)+\
                    '| deinterleave_fastq.py -i - -1 {trimed1} -2 {trimed2} '.format(trimed1=self.trimed1, trimed2=self.trimed2)
        self.run_process(command)

    def premap_tRNA_rRNA(self):
        '''
        premapping to tRNA and rRNA and make fastq files
        '''
        if not self.single_end:
            # paired end commands
            _input = '-1 {trimed1} -2 {trimed2} '.format(trimed1=self.trimed1, trimed2=self.trimed2)
            rRNA_extract = '| bamToFastq -fq {rRNA_FASTQ1} -fq2 {rRNA_FASTQ2} -i -'.format(rRNA_FASTQ1=self.rRNA_fastq1, rRNA_FASTQ2=self.rRNA_fastq2)
            tRNA_extract = '| bamToFastq -fq {TRNA_FASTQ1} -fq2 {TRNA_FASTQ2} -i -'.format(TRNA_FASTQ1=self.tRNA_fastq1, TRNA_FASTQ2=self.tRNA_fastq2)
            unmap_extract = '| bamToFastq -fq {PREMAP_FASTQ1} -fq2 {PREMAP_FASTQ2} -i - '.format(PREMAP_FASTQ1=self.premap_fastq1,
                                                                                    PREMAP_FASTQ2=self.premap_fastq2)
        else:
            # single end commands
            _input = '-U {trimed1} '.format(trimed1=self.trimed1)
            rRNA_extract = '| bamToFastq -fq {rRNA_FASTQ1} -i -'.format(rRNA_FASTQ1=self.rRNA_fastq1)
            tRNA_extract = '| bamToFastq -fq {TRNA_FASTQ1} -i -'.format(TRNA_FASTQ1=self.tRNA_fastq1)
            unmap_extract = '| bamToFastq -fq {PREMAP_FASTQ1} -i - '.format(PREMAP_FASTQ1=self.premap_fastq1)

        command =  'bowtie2 -p {threads} -D 20 -R 3 -N 0 -L 8 -i S,1,0.50 '.format(threads=self.threads)+\
                '--no-mixed --norc --no-discordant ' +\
                '-x {tRNA_rRNA_index} '.format(tRNA_rRNA_index = self.rRNA_tRNA_index) +\
                _input + \
                '| samtools view -bS@{threads} - '.format(threads=self.threads)+\
                '> {rRNA_tRNA_out}/tRNA_rRNA.bam'.format(rRNA_tRNA_out = self.rRNA_tRNA_out)
        self.run_process(command)


        ##extract tr/RNA reads
        command = 'samtools view -h -F4 {rRNA_tRNA_out}/tRNA_rRNA.bam '.format(rRNA_tRNA_out=self.rRNA_tRNA_out)+\
                "| awk '$1~\"^@\" || $2 ~/^(83|163|99|147|0|16)$/'" +\
                "| awk '$1~\"^@\" || $3~/gi\||rRNA/'" +\
                '| samtools view -b '+\
                rRNA_extract
        self.run_process(command)

        command = 'samtools view -h -F4 {rRNA_tRNA_out}/tRNA_rRNA.bam '.format(rRNA_tRNA_out=self.rRNA_tRNA_out)+\
                "| awk '$1~\"^@\" || $2 ~/^(83|163|99|147|0|16)$/'" +\
                "| awk '$1~\"^@\" || $3!~/gi\||rRNA/'" +\
                '| samtools view -b '+\
                tRNA_extract
        self.run_process(command)

        ##extract non tRNA/rRNA reads
        command = 'samtools view -bf4 {rRNA_tRNA_out}/tRNA_rRNA.bam'.format(rRNA_tRNA_out=self.rRNA_tRNA_out)+\
                unmap_extract
        self.run_process(command)


    def hisat_map(self):
        if not self.single_end:
            _input = '-1 {PREMAP_FASTQ1} -2 {PREMAP_FASTQ2} '.format(PREMAP_FASTQ1=self.premap_fastq1, PREMAP_FASTQ2=self.premap_fastq2)
            _split_option = ' '
            _unaligned = '| bamToFastq -i - -fq {bowtie_out}/unmapped.1.fq -fq2 {bowtie_out}/unmapped.2.fq'.format(bowtie_out=self.bowtie_out)
            _zip_command = 'gzip -f {bowtie_out}/unmapped.1.fq;gzip -f {bowtie_out}/unmapped.2.fq'.format(bowtie_out=self.bowtie_out)

        else:
            _input = '-U {PREMAP_FASTQ1} '.format(PREMAP_FASTQ1=self.premap_fastq1)
            _split_option = ' --single_end'
            _unaligned = '| bamToFastq -i - -fq {bowtie_out}/unmapped.1.fq '.format(bowtie_out=self.bowtie_out)
            _zip_command = 'gzip -f {bowtie_out}/unmapped.1.fq'.format(bowtie_out=self.bowtie_out)


        # map reads
        command = 'hisat2 -p {threads} -k 10 --no-mixed --no-discordant '.format(threads=self.threads)+\
                '--known-splicesite-infile {Splicesite} '.format(Splicesite=self.splicesite) +\
                '--novel-splicesite-outfile {hisat_out}/novelsite.txt -x {ref} '.format(hisat_out=self.hisat_out, ref=self.hisat_index)+\
                _input + \
                '| samtools view -bS - > {hisat_out}/hisat.bam'.format(hisat_out=self.hisat_out)
        self.run_process(command)

        #split to uniq and multimap
        uniq_command = 'split_uniq_bam.py -i {hisat_out}/hisat.bam -o {hisat_out}/hisat -a hisat2 {option}'.format(hisat_out=self.hisat_out, option=_split_option)
        self.run_process(uniq_command)

        #extract unaligned
        command = 'samtools view -@ {threads} -bf4 {hisat_out}/hisat.bam'.format(threads=self.threads, hisat_out=self.hisat_out) + \
                _unaligned
        self.run_process(command)
        self.run_process(_zip_command)

    def bowtie_map(self):

        if not self.single_end:
            _input = '-1 {fq_path}/unmapped.1.fq.gz -2 {fq_path}/unmapped.2.fq.gz'.format(fq_path = self.bowtie_out)
            _split_option = ' '

        else:
            _input = '-U {fq_path}/unmapped.1.fq.gz '.format(fq_path = self.bowtie_out)
            _split_option = ' --single_end'


        # map reads
        command = 'bowtie2 --mm  --local -D 20 -R 3 -N 0 -L 8 -i S,1,0.50 -p {threads} -k 10 '.format(threads=self.threads)+\
                '--no-mixed --no-discordant -x {index} '.format(index=self.bowtie2_index)+\
                _input +\
                '| samtools view -@{threads} -bS - > {bowtie_out}/bowtie2.bam'.format(threads=self.threads, bowtie_out=self.bowtie_out)
        self.run_process(command)

        # split to uniq and multimap
        command = 'split_uniq_bam.py -i {bowtie_out}/bowtie2.bam -o {bowtie_out}/bowtie -a bowtie2 {option}'.format(bowtie_out=self.bowtie_out, option=_split_option)
        self.run_process(command)

    def combined_aligned(self):
        _multi_option = '--single_end' if self.single_end else ' '
        _soft_clip_option = '--pe ' if not self.single_end else ' '
        command = 'samtools cat %s/hisat.multi.bam %s/bowtie.multi.bam ' %(self.hisat_out, self.bowtie_out)+\
                ' > %s/multi.bam' %(self.combined_out)
        self.run_process(command)

        command = 'reduce_multi_reads.py --infile {combined_out}/multi.bam --outfile {combined_out}/multi_filtered.bam '\
                                .format(combined_out = self.combined_out) +\
                        ' --bam_in --bam_out {option}'.format(option=_multi_option)
        self.run_process(command)

        command ='samtools cat %s/multi_filtered.bam %s/hisat.unique.bam %s/bowtie.unique.bam' %(self.combined_out, self.hisat_out, self.bowtie_out) +\
                '| filter_soft_clip.py -s 0.1 -b 0.2 -i - -o - %s' %_soft_clip_option + \
                '| samtools sort -n -@ %i -O bam -T %s/temp ' %(self.threads,self.combined_out) +\
                '> %s/primary.bam' %(self.combined_out)
        self.run_process(command)


    def dedup_bam(self):
        bam_file = self.combined_out + '/primary.bam'
        tag_bam = bam_file.replace('.bam','.tag.bam')
        sorted_bam = bam_file.replace('.bam','.sorted.bam')
        dedup_bam = bam_file.replace('.bam','.dedup.bam')
        umi_text = self.combined_out + '/primary.umi_metrics'
        duplicate_text = self.combined_out + '/primary.duplicate_metrics'

        umi_command = 'bam_umi_tag.py --in_bam %s --out_bam - --tag RX ' %bam_file + \
                        '| picard SortSam I=/dev/stdin O=/dev/stdout SORT_ORDER=queryname '+\
                        '| picard FixMateInformation ADD_MATE_CIGAR=true ASSUME_SORTED=true INPUT=/dev/stdin OUTPUT=/dev/stdout '   +\
                        '> %s' %(tag_bam)

        sort_command = 'samtools sort -@ %i -O bam -T %s/temp %s > %s' %(self.threads, self.combined_out, tag_bam, sorted_bam)
        index_command = 'samtools index %s' %(sorted_bam)
        dedup_command = 'picard UmiAwareMarkDuplicatesWithMateCigar UMI_METRICS_FILE=%s ' %(umi_text)+\
                        'MAX_EDIT_DISTANCE_TO_JOIN=1 TAG_DUPLICATE_SET_MEMBERS=true' +\
                        'UMI_TAG_NAME=RX INPUT=%s OUTPUT=%s ' %(sorted_bam, tag_bam) +\
                        'METRICS_FILE=%s REMOVE_DUPLICATES=false ASSUME_SORT_ORDER=coordinate' %(duplicate_text)
        resort_command = 'samtools sort -n@  %i -O bam -T %s/temp %s > %s ' %(self.threads, self.combined_out, tag_bam, dedup_bam)

        self.run_process(umi_command)
        self.run_process(sort_command)
        self.run_process(index_command)
        self.run_process(dedup_command)
        self.run_process(resort_command)

    def combined_filter(self):

        if not self.single_end:
            _verb = 'pairtobed'
            _option = '-type neither'

        else:
            _verb = 'intersect'
            _option = '-v'

        ### filter out tRNA
        command = 'bedtools {verb} -s -f 0.01 -abam {combined_path}/primary.bam -b {bed_path}/tRNA.bed > {tRNA_path}/tRNA_primary.bam' \
                        .format(combined_path = self.combined_out, bed_path = self.bedpath, tRNA_path=self.tRNA_out, verb = _verb)
        self.run_process(command)

        ### filter out rRNA
        command = 'bedtools {verb} -s -f 0.01 -abam {combined_path}/primary.bam -b {bed_path}/rRNA_for_bam_filter.bed > {rRNA_path}/rRNA_primary.bam' \
                        .format(combined_path = self.combined_out, bed_path = self.bedpath, rRNA_path=self.rRNA_out, verb = _verb)
        self.run_process(command)

        ### filter out tRNA
        command = 'bedtools {verb} -s -f 0.01 -abam {count_bam} -b {bed_path}/sncRNA_no_tRNA.bed > {combined_path}/sncRNA.bam' \
                        .format(count_bam = self.count_bam, bed_path = self.bedpath, combined_path = self.combined_out, verb = _verb)
        self.run_process(command)

        ### filter out long RNA
        command = 'bedtools {verb} -s -f 0.01 {option} -abam {count_bam} -b {bed_path}/sncRNA_rRNA_for_bam_filter.bed > {combined_path}/primary_no_sncRNA_tRNA_rRNA.bam' \
                        .format(count_bam = self.count_bam, bed_path = self.bedpath, combined_path = self.combined_out, option=_option, verb = _verb)
        self.run_process(command)

        if self.rmsk:
            command = 'bedtools {verb} -f 0.01 {option} -abam {combined_path}/primary_no_sncRNA_tRNA_rRNA.bam -b {rmsk_bed} > {combined_path}/primary_no_sncRNA_tRNA_rRNA_repeats.bam' \
                .format(count_bam = self.count_bam, option=_option, bed_path = self.bedpath, combined_path = self.combined_out, verb = _verb, rmsk_bed = self.rmsk)
            self.run_process(command)

            command = 'bedtools {verb} -f 0.01 -abam {combined_path}/primary_no_sncRNA_tRNA_rRNA.bam -b {rmsk_bed} > {combined_path}/repeats.bam' \
                .format(count_bam = self.count_bam, bed_path = self.bedpath, combined_path = self.combined_out, verb = _verb, rmsk_bed = self.rmsk)
            self.run_process(command)


    def make_alignment_bed(self):

        if not self.single_end:
            command = 'bam_to_bed.py -m 5 -M 1000000  -i {combined_out}/sncRNA.bam > {combined_out}/sncRNA.bed'.format(combined_out=self.combined_out)
            self.run_process(command)

            command = 'bam_to_bed.py -i {combined_out}/primary_no_sncRNA_tRNA_rRNA.bam  -m 5 -M 1000000 > {combined_out}/primary_no_sRNAs.bed'.format(combined_out=self.combined_out)
            self.run_process(command)
        else:
            command = 'bedtools bamtobed  -i {combined_out}/sncRNA.bam > {combined_out}/sncRNA.bed'.format(combined_out=self.combined_out)
            self.run_process(command)

            command = 'bedtools bamtobed -i {combined_out}/primary_no_sncRNA_tRNA_rRNA.bam  > {combined_out}/primary_no_sRNAs.bed'.format(combined_out=self.combined_out)
            self.run_process(command)


    def generate_tRNA_remap(self):
        # tRNA reads process
        if not self.single_end:
            command = 'bamToFastq -i {tRNA_path}/tRNA_primary.bam -fq {tRNA_path}/tRNA.1.fq -fq2 {tRNA_path}/tRNA.2.fq'.format(tRNA_path =self.tRNA_out)
            self.run_process(command)
            command = 'cat {tRNA_fastq1} {tRNA_path}/tRNA.1.fq | gzip > {tRNA_path}/tRNA.1.fq.gz'.format(tRNA_fastq1=self.tRNA_fastq1, tRNA_path=self.tRNA_out)
            self.run_process(command)
            command = 'cat {tRNA_fastq2} {tRNA_path}/tRNA.2.fq | gzip > {tRNA_path}/tRNA.2.fq.gz'.format(tRNA_fastq2=self.tRNA_fastq2, tRNA_path=self.tRNA_out)
            self.run_process(command)

            _input = '-1 {tRNA_path}/tRNA.1.fq.gz -2 {tRNA_path}/tRNA.2.fq.gz '.format(tRNA_path=self.tRNA_out)

        else:
            command = 'bamToFastq -i {tRNA_path}/tRNA_primary.bam -fq {tRNA_path}/tRNA.1.fq'.format(tRNA_path =self.tRNA_out)
            self.run_process(command)
            command = 'cat {tRNA_fastq1} {tRNA_path}/tRNA.1.fq | gzip > {tRNA_path}/tRNA.1.fq.gz'.format(tRNA_fastq1=self.tRNA_fastq1, tRNA_path=self.tRNA_out)
            self.run_process(command)

            _input = '-U {tRNA_path}/tRNA.1.fq.gz'.format(tRNA_path=self.tRNA_out)

        command = 'bowtie2 -p {threads} --local -D 20 -R 3 -N 0 -L 8 -i S,1,0.50 '.format(threads=self.threads)+\
                        '--norc --no-mixed --no-discordant -x {tRNA_index} '.format(tRNA_index=self.tRNA_index)+\
                        _input +\
                        '| samtools view -bS@ {threads} - > {tRNA_path}/tRNA_remap.bam'.format(tRNA_path=self.tRNA_out, threads=self.threads)
        self.run_process(command)

        if self.UMI > 0 and not self.count_all:
            command = ' bam_umi_tag.py --in_bam %s/tRNA_remap.bam --out_bam - --tag RX ' %(self.tRNA_out)+\
                    '| picard SortSam I=/dev/stdin O=/dev/stdout SORT_ORDER=queryname '+\
                    '| picard FixMateInformation ADD_MATE_CIGAR=true ASSUME_SORTED=true INPUT=/dev/stdin OUTPUT=/dev/stdout ' +\
                    '| samtools sort -@ {threads} -T {tRNA_path}/tRNA -O bam > {tRNA_path}/tRNA_remap.sort.bam'.format(threads=self.threads, tRNA_path=self.tRNA_out)
            self.run_process(command)
            command = 'samtools index {tRNA_path}/tRNA_remap.sort.bam'.format(tRNA_path=self.tRNA_out)
            self.run_process(command)
            dedup_command = 'picard UmiAwareMarkDuplicatesWithMateCigar UMI_METRICS_FILE=%s/tRNA.umi_metric ' %(self.tRNA_out)+\
                            'MAX_EDIT_DISTANCE_TO_JOIN=1 TAG_DUPLICATE_SET_MEMBERS=true ' +\
                            'UMI_TAG_NAME=RX INPUT={tRNA_path}/tRNA_remap.sort.bam OUTPUT=/dev/stdout '.format(tRNA_path=self.tRNA_out) +\
                            'METRICS_FILE=%s/tRNA.duplicate_metrics REMOVE_DUPLICATES=false ASSUME_SORT_ORDER=coordinate ' %(self.tRNA_out)   +\
                            '| samtools sort -n@ {threads} -T {tRNA_path}/tRNA -O bam - > {tRNA_path}/tRNA_remap.dedup.bam'.format(threads=self.threads, tRNA_path=self.tRNA_out)
            self.run_process(dedup_command)
            command = 'cat {tRNA_path}/tRNA_remap.dedup.bam | bam_to_bed.py -i - -o {tRNA_path}/tRNA.bed -m 5 -M 10000'.format(tRNA_path=self.tRNA_out)
            self.run_process(command)
        else:
            if not self.single_end:
                command = 'bam_to_bed.py -i {tRNA_path}/tRNA_remap.bam  -o {tRNA_path}/tRNA.bed -m 5 -M 10000'.format(tRNA_path=self.tRNA_out)
            else:
                command = 'bedtools bamtobed -i {tRNA_path}/tRNA_remap.bam  > {tRNA_path}/tRNA.bed'.format(tRNA_path=self.tRNA_out)
            self.run_process(command)

    def generate_tRNA_count(self):
        tRNA_count = defaultdict(int)
        tRNA_bed = self.tRNA_out + '/tRNA.bed'
        tRNA_count_file = self.tRNA_raw + '/' + self.samplename + '.tRNA'
        print('Reading from %s' %tRNA_bed, file=sys.stderr)
        if not self.dry:
            with open(tRNA_bed,'r') as bed:
                for line in bed:
                    tRNA=line.split('\t')[0]
                    tRNA_count[tRNA] += 1

            with open(tRNA_count_file, 'w') as count_file:
                for key, value in tRNA_count.iteritems():
                    count_file.write('%s\t%i\n' %(key, value))
        print('Written %s' %tRNA_count_file, file=sys.stderr)

    def generate_rRNA_count(self):

        if not self.single_end:
            command = 'bamToFastq -fq {rRna_path}/rRNA.1.fq -fq2 {rRna_path}/rRNA.2.fq -i {rRna_path}/rRNA_primary.bam'.format(rRna_path=self.rRNA_out)
            self.run_process(command)
            command = 'cat {rRNA_fastq1} {rRNA_path}/rRNA.1.fq | gzip > {rRNA_path}/rRNA.1.fq.gz'.format(rRNA_fastq1=self.rRNA_fastq1, rRNA_path=self.rRNA_out)
            self.run_process(command)
            command = 'cat {rRNA_fastq2} {rRNA_path}/rRNA.2.fq | gzip > {rRNA_path}/rRNA.2.fq.gz'.format(rRNA_fastq2=self.rRNA_fastq2, rRNA_path=self.rRNA_out)
            self.run_process(command)

            _input = '-1 {rRNA_path}/rRNA.1.fq.gz -2 {rRNA_path}/rRNA.2.fq.gz '.format(rRNA_path=self.rRNA_out)
        else:
            command = 'bamToFastq -fq {rRna_path}/rRNA.1.fq  -i {rRna_path}/rRNA_primary.bam'.format(rRna_path=self.rRNA_out)
            self.run_process(command)
            command = 'cat {rRNA_fastq1} {rRNA_path}/rRNA.1.fq | gzip > {rRNA_path}/rRNA.1.fq.gz'.format(rRNA_fastq1=self.rRNA_fastq1, rRNA_path=self.rRNA_out)
            self.run_process(command)

            _input = '-U {rRNA_path}/rRNA.1.fq.gz'.format(rRNA_path=self.rRNA_out)

        command = 'bowtie2 -p {threads} -D 20 -R 3 -N 0 -L 8 -i S,1,0.50 '.format(threads=self.threads)+\
                '--no-mixed --no-discordant -x {rRNA_index} '.format(rRNA_index=self.rRNA_index) +\
                _input +\
                '| samtools view -bS@ {threads} - > {rRNA_path}/rRNA_remap.bam'.format(rRNA_path=self.rRNA_out, threads=self.threads)
        self.run_process(command)

        if self.UMI > 0 and not self.count_all:
            command = ' bam_umi_tag.py --in_bam %s/rRNA_remap.bam --out_bam - --tag RX ' %(self.rRNA_out)+\
                    '| picard SortSam I=/dev/stdin O=/dev/stdout SORT_ORDER=queryname '+\
                    '| picard FixMateInformation ADD_MATE_CIGAR=true ASSUME_SORTED=true INPUT=/dev/stdin OUTPUT=/dev/stdout'+\
                    '| samtools sort -@ {threads} -T {rRNA_path}/rRNA -O bam > {rRNA_path}/rRNA_remap.sort.bam'.format(threads=self.threads, rRNA_path=self.rRNA_out)
            self.run_process(command)
            command = 'samtools index {rRNA_path}/rRNA_remap.sort.bam'.format(rRNA_path=self.rRNA_out)
            self.run_process(command)
            dedup_command = 'picard UmiAwareMarkDuplicatesWithMateCigar UMI_METRICS_FILE=%s/rRNA.umi_metric ' %(self.rRNA_out)+\
                            'MAX_EDIT_DISTANCE_TO_JOIN=1 TAG_DUPLICATE_SET_MEMBERS=true' +\
                            'UMI_TAG_NAME=RX INPUT={rRNA_path}/rRNA_remap.sort.bam OUTPUT=/dev/stdout '.format(rRNA_path=self.rRNA_out) +\
                            'METRICS_FILE=%s/rRNA.duplicate_metrics REMOVE_DUPLICATES=false ASSUME_SORT_ORDER=coordinate' %(self.rRNA_out)   +\
                            '| samtools sort -n@ {threads} -T {rRNA_path}/rRNA -O bam - > {rRNA_path}/rRNA_remap.dedup.bam'.format(threads=self.threads, rRNA_path=self.rRNA_out)
            self.run_process(dedup_command)
            command = 'cat {rRNA_path}/rRNA_remap.dedup.bam | bam_to_bed.py -i - -o {rRNA_path}/rRNA.bed -m 5 -M 10000'.format(rRNA_path=self.rRNA_out)
            self.run_process(command)
        else:
            if not self.single_end:
                command = 'bam_to_bed.py -i {rRNA_path}/rRNA_remap.bam  -o {rRNA_path}/rRNA.bed -m 5 -M 10000'.format(rRNA_path=self.rRNA_out)
            else:
                command = 'bedtools bamtobed -i {rRNA_path}/rRNA_remap.bam  > {rRNA_path}/rRNA.bed'.format(rRNA_path=self.rRNA_out)
            self.run_process(command)

        command = 'bedtools coverage -s -counts -F 0.1 -a {bed_path}/rRNA.bed -b {rRNA_path}/rRNA.bed > {rRNA_path}/rRNA.counts'.format(bed_path=self.bedpath, rRNA_path=self.rRNA_out)
        self.run_process(command)

    def generate_all_count(self):
        command = 'bedtools coverage -s -counts -F 0.1 -a {bed_path}/sncRNA_no_tRNA.bed '.format(bed_path=self.bedpath)+\
                '-b {combined}/sncRNA.bed > {combined}/sncRNA.counts'.format(combined=self.combined_out)
        self.run_process(command)
        command = 'bedtools coverage -s -counts -F 0.1 -a {bed_path}/genes_no_sncRNA_rRNA_tRNA.bed '.format(bed_path=self.bedpath)+\
                '-b {combined}/primary_no_sRNAs.bed > {combined}/non_sRNAs.counts'.format(combined=self.combined_out)
        self.run_process(command)

        command = 'cat {combined}/non_sRNAs.counts {combined}/sncRNA.counts {rRNA_path}/rRNA.counts '.format(combined=self.combined_out, rRNA_path=self.rRNA_out) +\
                '> {count_path}/{samplename}.counts'.format(count_path=self.count_raw, samplename = self.samplename)
        self.run_process(command)

    def generate_repeat_count(self):
        # repeat reads process
        if not self.single_end:
            command = 'samtools fastq -N@ {threads} {combined_path}/repeats.bam '\
                        .format(combined_path =self.combined_out, threads = self.threads)  +\
                '> {repeat_path}/repeats.fq'.format(repeat_path=self.repeat_out)
            self.run_process(command)
            _option=' --interleaved '

        else:
            command = 'samtools fastq -@ {threads} {combined_path}/repeats.bam '\
                        .format(combined_path =self.combined_out, threds=self.threds)  +\
                '> {repeat_path}/repeats.fq'.format(repeat_path=self.repeat_out)
            self.run_process(command)
            _option=' -U '

        command = 'bowtie2 -p {threads} --local -D 20 -R 3 -N 0 -L 8 -i S,1,0.50 '.format(threads=self.threads)+\
                        '--no-mixed --no-discordant -x {repeat_index} '.format(repeat_index=self.rmsk_index)+\
                        ' {option} {repeat_path}/repeats.fq'.format(option=_option, repeat_path=self.repeat_out)+\
                        '| samtools view -bS@ {threads} - > {repeat_path}/repeat_remap.bam'.format(repeat_path=self.repeat_out, threads=self.threads)
        self.run_process(command)

        if not self.single_end:
            command = 'bam_to_bed.py -i {repeat_path}/repeat_remap.bam  -o {repeat_path}/repeat.bed -m 5 -M 10000'.format(repeat_path=self.repeat_out)
        else:
            command = 'bedtools bamtobed -i {repeat_path}/repeat_remap.bam  > {repeat_path}/repeat.bed'.format(repeat_path=self.repeat_out)
        self.run_process(command)

        repeat_count = defaultdict(int)
        repeat_bed = self.repeat_out + '/repeat.bed'
        repeat_count_file = self.count_rmsk + '/' + self.samplename + '.repeat'
        print('Reading from %s' %repeat_bed, file=sys.stderr)
        if not self.dry:
            with open(repeat_bed,'r') as bed:
                for line in bed:
                    repeat=line.split('\t')[0]
                    repeat_count[repeat] += 1

            with open(repeat_count_file, 'w') as count_file:
                for key, value in repeat_count.iteritems():
                    count_file.write('%s\t%i\n' %(key, value))
        print('Written %s' %repeat_count_file, file=sys.stderr)
 

def system_run(dry, samplename, command):
    print('[%s] Running: %s' %(samplename, command), file=sys.stderr)
    if dry:
            return 0

    else:
            start = time.time()
            os.system(command)
            end = time.time() - start
            print('[%s] Used time %.3f min' %(samplename, end/60), file=sys.stderr)
            return 0


def makeFolder(folder):
    """
            Input a folder name and make a folder if it is non-existed
    """
    sys.stderr.write('Creating %s....\n' %folder)
    if os.path.isdir(folder):
            sys.stderr.write('%s exists.\n' %folder)
    else:
            os.mkdir(folder)
            sys.stderr.write('Created %s.\n' %folder)
    return 0
