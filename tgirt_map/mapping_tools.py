from __future__ import division
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
        self.UMI = args.umi
        self.TTN = args.TTN
        self.dry = args.dry

        #### make folder
        self.trim_folder = self.outpath + '/Trim'
        self.count_folder= self.outpath + '/Counts'
        self.count_raw = self.count_folder + '/RAW'
        self.tRNA_raw = self.count_folder + '/tRNA_RAW'

        #define sample folder
        self.samplename = re.sub('.fastq|.gz|.fq','',self.fastq1.split('/')[-1])
        self.sample_folder = self.outpath + '/' + self.samplename
        self.hisat_out = self.sample_folder + '/Hisat'
        self.rRNA_tRNA_out = self.sample_folder + '/rRNA_tRNA_premap'
        self.bowtie_out = self.sample_folder + '/Bowtie'
        self.combined_out = self.sample_folder + '/Combined'
        self.tRNA_out = self.sample_folder + '/tRNA'
        self.rRNA_out = self.sample_folder + '/rRNA'

        #make output file names
        self.trimed1= '%s/%s.1.fq.gz' %(self.trim_folder, self.samplename)
        self.trimed2= self.trimed1.replace('.1.fq.gz','.2.fq.gz')

        self.premap_fastq1 = '%s/non_tRNA_rRNA.1.fq' %self.rRNA_tRNA_out
        self.premap_fastq2 = '%s/non_tRNA_rRNA.2.fq' %self.rRNA_tRNA_out

        self.tRNA_fastq1 = '%s/tRNA.1.fq' %self.rRNA_tRNA_out
        self.tRNA_fastq2 = '%s/tRNA.2.fq' %self.rRNA_tRNA_out

        self.rRNA_fastq1 = '%s/rRNA.1.fq' %self.rRNA_tRNA_out
        self.rRNA_fastq2 = '%s/rRNA.2.fq' %self.rRNA_tRNA_out

        if self.UMI == 0:
                self.count_bam = self.combined_out + '/primary.bam'
        else:
                self.count_bam = self.combined_out + '/primary.dedup.bam'

        self.run_process = partial(system_run, args.dry, self.samplename)

    def make_result_dir(self):
        folders = [self.outpath, self.trim_folder, self.count_folder, self.count_raw,
                         self.tRNA_raw, self.sample_folder, self.hisat_out, self.rRNA_tRNA_out,
                        self.bowtie_out, self.combined_out, self.tRNA_out, self.rRNA_out]
        mf = map(makeFolder, folders)


    def trimming(self):
        R1 = 'AAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
        R2 = 'GATCGTCGGACTGTAGAACTCTGAACGTGTAGA'
        if self.TTN:
                option='-U 1'
        else:
                option = '' 

        if self.UMI == 0:
            command = 'cutadapt -m 15 -O 5 -n 3 {option} -q 20 -b {R1} -B {R2} -o {trimed1} -p {trimed2} {file1} {file2}'\
                    .format(R1=R1, R2=R2, option= option,
                            trimed1=self.trimed1, trimed2=self.trimed2,
                            file1= self.fastq1, file2= self.fastq2)
        else:
            command = 'clip_fastq.py --fastq1={file1} --fastq2={file2} --idxBase={umi} '.format(file1= self.fastq1, file2= self.fastq2, umi=self.UMI*'X')+\
                        '--barcodeCutOff=20 --outputprefix=- --prefix_split=0 -r read1 '+\
                    '| cutadapt -m 15 -O 5 -n 3 {option} -q 20 -b {R1} -B {R2} --interleaved --quiet - '.format(option=option, R1=R1,R2=R2)+\
                    '| deinterleave_fastq.py - {trimed1} {trimed2} '.format(trimed1=self.trimed1, trimed2=self.trimed2)
        self.run_process(command)

    def premap_tRNA_rRNA(self):
        command = 'bowtie2 -p {threads} -D 20 -R 3 -N 0 -L 8 -i S,1,0.50 '.format(threads=self.threads)+\
                '--no-mixed --norc --no-discordant ' +\
                '-x {tRNA_rRNA_index} -1 {trimed1} -2 {trimed2} '.format(tRNA_rRNA_index = self.rRNA_tRNA_index,
                                                                        trimed1=self.trimed1, trimed2=self.trimed2) +\
                '| samtools view -bS@{threads} - '.format(threads=self.threads)+\
                '> {rRNA_tRNA_out}/tRNA_rRNA.bam'.format(rRNA_tRNA_out = self.rRNA_tRNA_out)
        self.run_process(command)


        ##extract tr/RNA reads
        command = 'samtools view -h -F4 {rRNA_tRNA_out}/tRNA_rRNA.bam '.format(rRNA_tRNA_out=self.rRNA_tRNA_out)+\
                "| awk '$1~\"^@\" || $2 == 83 || $2 == 163 || $2 == 99 || $2 == 147'" +\
                "| awk '$1~\"^@\" || $3~/gi\||rRNA/'" +\
                '| samtools view -b '+\
                '| bamToFastq -fq {rRNA_FASTQ1} -fq2 {rRNA_FASTQ2} -i -'.format(rRNA_FASTQ1=self.rRNA_fastq1, rRNA_FASTQ2=self.rRNA_fastq2)
        self.run_process(command)

        command = 'samtools view -h -F4 {rRNA_tRNA_out}/tRNA_rRNA.bam '.format(rRNA_tRNA_out=self.rRNA_tRNA_out)+\
                "| awk '$1~\"^@\" || $2 == 83 || $2 == 163 || $2 == 99 || $2 == 147'" +\
                "| awk '$1~\"^@\" || $3!~/gi\||rRNA/'" +\
                '| samtools view -b '+\
                '| bamToFastq -fq {TRNA_FASTQ1} -fq2 {TRNA_FASTQ2} -i -'.format(TRNA_FASTQ1=self.tRNA_fastq1, TRNA_FASTQ2=self.tRNA_fastq2)
        self.run_process(command)

        ##extract non tRNA/rRNA reads
        command = 'samtools view -bf4 {rRNA_tRNA_out}/tRNA_rRNA.bam'.format(rRNA_tRNA_out=self.rRNA_tRNA_out)+\
                '| bamToFastq -fq {PREMAP_FASTQ1} -fq2 {PREMAP_FASTQ2} -i - '.format(PREMAP_FASTQ1=self.premap_fastq1,
                                                                                    PREMAP_FASTQ2=self.premap_fastq2)
        self.run_process(command)

    def hisat_map(self):
        # map reads
        command = 'hisat2 -p {threads} -k 10 --no-mixed --no-discordant '.format(threads=self.threads)+\
                '--known-splicesite-infile {Splicesite} '.format(Splicesite=self.splicesite) +\
                '--novel-splicesite-outfile {hisat_out}/novelsite.txt -x {ref} -1 {PREMAP_FASTQ1} -2 {PREMAP_FASTQ2} '\
                     .format(hisat_out=self.hisat_out, ref=self.hisat_index, PREMAP_FASTQ1=self.premap_fastq1, PREMAP_FASTQ2=self.premap_fastq2) +\
                '| samtools view -bS - > {hisat_out}/hisat.bam'.format(hisat_out=self.hisat_out)
        self.run_process(command)

        #split to uniq and multimap
        command = 'split_uniq_bam.py -i {hisat_out}/hisat.bam -o {hisat_out}/hisat -a hisat2'.format(hisat_out=self.hisat_out)
        self.run_process(command)

        #extract unaligned
        command = 'samtools view -@ {threads} -bf4 {hisat_out}/hisat.bam'.format(threads=self.threads, hisat_out=self.hisat_out) + \
                        '| bamToFastq -i - -fq {bowtie_out}/unmapped.1.fq -fq2 {bowtie_out}/unmapped.2.fq'.format(bowtie_out=self.bowtie_out)
        self.run_process(command)
        command = 'gzip -f {bowtie_out}/unmapped.1.fq;gzip -f {bowtie_out}/unmapped.2.fq'.format(bowtie_out=self.bowtie_out)
        self.run_process(command)

    def bowtie_map(self):
        # map reads
        command = 'bowtie2 --local -D 20 -R 3 -N 0 -L 8 -i S,1,0.50 -p {threads} -k 10 '.format(threads=self.threads)+\
                '--no-mixed --no-discordant -x {index} -1 {fq_path}/unmapped.1.fq.gz -2 {fq_path}/unmapped.2.fq.gz'\
                        .format(index=self.bowtie2_index, fq_path = self.bowtie_out) +\
                '| samtools view -@{threads} -bS - > {bowtie_out}/bowtie2.bam'.format(threads=self.threads, bowtie_out=self.bowtie_out)
        self.run_process(command)

        # split to uniq and multimap
        command = 'split_uniq_bam.py -i {bowtie_out}/bowtie2.bam -o {bowtie_out}/bowtie -a bowtie2'.format(bowtie_out=self.bowtie_out)
        self.run_process(command)

    def combined_aligned(self):
        command = 'samtools cat %s/hisat.multi.bam %s/bowtie.multi.bam ' %(self.hisat_out, self.bowtie_out)+\
                ' > %s/multi.bam' %(self.combined_out)
        self.run_process(command)

        command = 'reduce_multi_reads.py --infile {combined_out}/multi.bam --outfile {combined_out}/multi_filtered.bam '\
                                .format(combined_out = self.combined_out) +\
                        ' --bam_in --bam_out'
        self.run_process(command)

        command ='samtools cat %s/multi_filtered.bam %s/hisat.unique.bam %s/bowtie.unique.bam' %(self.combined_out, self.hisat_out, self.bowtie_out) +\
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

        ### filter out tRNA
        command = 'bedtools pairtobed -s -f 0.01 -abam {combined_path}/primary.bam -b {bed_path}/tRNA.bed > {tRNA_path}/tRNA_primary.bam' \
                        .format(combined_path = self.combined_out, bed_path = self.bedpath, tRNA_path=self.tRNA_out)
        self.run_process(command)

        ### filter out rRNA
        command = 'bedtools pairtobed -s -f 0.01 -abam {combined_path}/primary.bam -b {bed_path}/rRNA_for_bam_filter.bed > {rRNA_path}/rRNA_primary.bam' \
                        .format(combined_path = self.combined_out, bed_path = self.bedpath, rRNA_path=self.rRNA_out)
        self.run_process(command)

        ### filter out tRNA
        command = 'bedtools pairtobed -s -f 0.01 -abam {count_bam} -b {bed_path}/sncRNA_no_tRNA.bed > {combined_path}/sncRNA.bam' \
                        .format(count_bam = self.count_bam, bed_path = self.bedpath, combined_path = self.combined_out)
        self.run_process(command)

        ### filter out long RNA
        command = 'bedtools pairtobed -s -f 0.01 -type neither -abam {count_bam} -b {bed_path}/sncRNA_rRNA_for_bam_filter.bed > {combined_path}/primary_no_sncRNA_tRNA_rRNA.bam' \
                        .format(count_bam = self.count_bam, bed_path = self.bedpath, combined_path = self.combined_out)
        self.run_process(command)


    def make_alignment_bed(self):
        command = 'bam_to_bed.py -i {combined_out}/sncRNA.bam -o {combined_out}/sncRNA.bed -m 5 -M 1000000'.format(combined_out=self.combined_out)
        self.run_process(command)

        command = 'bam_to_bed.py -i {combined_out}/primary_no_sncRNA_tRNA_rRNA.bam  -o {combined_out}/primary_no_sRNAs.bed -m 5 -M 1000000'.format(combined_out=self.combined_out)
        self.run_process(command)

    def generate_tRNA_remap(self):
        # tRNA reads process
        command = 'bamToFastq -i {tRNA_path}/tRNA_primary.bam -fq {tRNA_path}/tRNA.1.fq -fq2 {tRNA_path}/tRNA.2.fq'.format(tRNA_path =self.tRNA_out)
        self.run_process(command)
        command = 'cat {tRNA_fastq1} {tRNA_path}/tRNA.1.fq | gzip > {tRNA_path}/tRNA.1.fq.gz'.format(tRNA_fastq1=self.tRNA_fastq1, tRNA_path=self.tRNA_out)
        self.run_process(command)
        command = 'cat {tRNA_fastq2} {tRNA_path}/tRNA.2.fq | gzip > {tRNA_path}/tRNA.2.fq.gz'.format(tRNA_fastq2=self.tRNA_fastq2, tRNA_path=self.tRNA_out)
        self.run_process(command)

        command = 'bowtie2 -p {threads} --local -D 20 -R 3 -N 0 -L 8 -i S,1,0.50 '.format(threads=self.threads)+\
                        '--norc --no-mixed --no-discordant -x {tRNA_index} '.format(tRNA_index=self.tRNA_index)+\
                        '-1 {tRNA_path}/tRNA.1.fq.gz -2 {tRNA_path}/tRNA.2.fq.gz '.format(tRNA_path=self.tRNA_out)+\
                        '| samtools view -bS@ {threads} - > {tRNA_path}/tRNA_remap.bam'.format(tRNA_path=self.tRNA_out, threads=self.threads)
        self.run_process(command)

        if self.UMI > 0:
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
                command = 'bam_to_bed.py -i {tRNA_path}/tRNA_remap.bam  -o {tRNA_path}/tRNA.bed -m 5 -M 10000'.format(tRNA_path=self.tRNA_out)
                self.run_process(command)

    def generate_tRNA_count(self):
        tRNA_count = defaultdict(int)
        tRNA_bed = self.tRNA_out + '/tRNA.bed'
        tRNA_count_file = self.tRNA_raw + '/' + self.samplename + '.tRNA'
        print >> sys.stderr, 'Reading from %s' %tRNA_bed
        if not self.dry:
            with open(tRNA_bed,'r') as bed:
                for line in bed:
                    tRNA=line.split('\t')[0]
                    tRNA_count[tRNA] += 1

            with open(tRNA_count_file, 'w') as count_file:
                for key, value in tRNA_count.iteritems():
                    count_file.write('%s\t%i\n' %(key, value))
        print >> sys.stderr, 'Written %s' %tRNA_count_file

    def generate_rRNA_count(self):
        command = 'bamToFastq -fq {rRna_path}/rRNA.1.fq -fq2 {rRna_path}/rRNA.2.fq -i {rRna_path}/rRNA_primary.bam'.format(rRna_path=self.rRNA_out)
        self.run_process(command)
        command = 'cat {rRNA_fastq1} {rRNA_path}/rRNA.1.fq | gzip > {rRNA_path}/rRNA.1.fq.gz'.format(rRNA_fastq1=self.rRNA_fastq1, rRNA_path=self.rRNA_out)
        self.run_process(command)
        command = 'cat {rRNA_fastq2} {rRNA_path}/rRNA.2.fq | gzip > {rRNA_path}/rRNA.2.fq.gz'.format(rRNA_fastq2=self.rRNA_fastq2, rRNA_path=self.rRNA_out)
        self.run_process(command)

        command = 'bowtie2 -p {threads} -D 20 -R 3 -N 0 -L 8 -i S,1,0.50 '.format(threads=self.threads)+\
                '--no-mixed --no-discordant -x {rRNA_index} '.format(rRNA_index=self.rRNA_index) +\
                '-1 {rRNA_path}/rRNA.1.fq.gz -2 {rRNA_path}/rRNA.2.fq.gz '.format(rRNA_path=self.rRNA_out)+\
                '| samtools view -bS@ {threads} - > {rRNA_path}/rRNA_remap.bam'.format(rRNA_path=self.rRNA_out, threads=self.threads)
        self.run_process(command)

        if self.UMI > 0:
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
                command = 'bam_to_bed.py -i {rRNA_path}/rRNA_remap.bam  -o {rRNA_path}/rRNA.bed -m 5 -M 10000'.format(rRNA_path=self.rRNA_out)
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

def system_run(dry, samplename, command):
    print >> sys.stderr, '[%s] Running: %s' %(samplename, command)
    if dry:
            return 0

    else:
            start = time.time()
            os.system(command)
            end = time.time() - start
            print >> sys.stderr, '[%s] Used time %.3f min' %(samplename, end/60)
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
