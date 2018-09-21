#!/bin/env  python
# This is the pair-end pipeline for tgirt sequencing
# Mapping with hisat + bowtie local
# and extract tRNA reads for reassigning counts


from __future__ import division, print_function
import time
import glob
import argparse
import sys
from tgirt_map.mapping_tools import sample_object

def getopt():
    parser = argparse.ArgumentParser(description='Pipeline for mapping and counting for TGIRT-seq paired end data')
    parser.add_argument('-1', '--fastq1', 
              help = 'pairedEnd fastq file (read1)', required=True)
    parser.add_argument('-2', '--fastq2', 
              help = 'pairedEnd fastq file (read2). Optional: only needed for paired-end', required=True)
    parser.add_argument('-o','--outdir', 
              help = 'result directory that all resulting/intermediate files will be stored\n' + \
                                         'will create 1. $resultpath/trimmed\n' + \
                                         '            2. $resultpath/hisat\n'  + \
                                         '            3. $resultpath/bowtie2\n' + \
                                         '            4. $resultpath/mergeBam (all useful result files)\n',
                        required=True)
    parser.add_argument('--hisat_index', 
              help = 'hisat2 index', required=True)
    parser.add_argument('--bowtie2_index', 
              help = 'bowtie2 index', required=True)
    parser.add_argument('--univec_index', 
              help = 'bowtie2 index for univec sequences')
    parser.add_argument('--smRNA_index', 
              help = 'bowtie2 index for small RNA sequences')
    parser.add_argument('--bedpath', 
              help = 'bed folder for gene counting', required=True)
    parser.add_argument('--splicesite', 
              help = 'splice site file generated by hisat', required=True)
    parser.add_argument('--rRNA_mt_index' , 
              help = 'bowtie2 index for rRNA and tRNA combined', required=True)
    parser.add_argument('-p', '--threads', default=1, type=int, 
              help = 'number of cores to be used for the pipeline (default:1)')
    parser.add_argument('--TTN', action='store_true',  
              help = 'used TTN primer')
    parser.add_argument('--novel_splice', action='store_true',  
              help = 'Tweak hisat2 to enable more novel splicings')
    parser.add_argument('--trim_aggressive', action='store_true',  
              help = 'trim all R2R to R2R junction')
    parser.add_argument('--polyA', action='store_true',  
              help = 'trim for polyA/T')
    parser.add_argument('--umi', default=0, type=int,
              help = "Number of UMI bases from 5' of R1 (default = 0)")
    parser.add_argument('--count_all', action='store_true',
              help = "Ignore UMI for counting, only evaluated with --umi option")
    parser.add_argument('--repeats', default=None,
              help = "Repeat mask BED file, recount repeat masks if a BED file is given (default: null)")
    parser.add_argument('--repeats_index', default=None,
              help = "Bowtie2 index of repeat mask fasta file (default: null)")
    parser.add_argument('--dry', action='store_true', help = "DEBUG: Dry run")
    parser.add_argument('--skip_trim', action='store_true',  
              help = 'DEBUG: skip trimming')
    parser.add_argument('--skip_univec', action='store_true',  
              help = 'DEBUG: skip univec filter')
    parser.add_argument('--skip_smRNA', action='store_true',  
              help = 'DEBUG: skip smallRNA filter')
    parser.add_argument('--skip_premap', action='store_true',  
              help = 'DEBUG: skip premapping tRNA and rRNA')
    parser.add_argument('--skip_hisat', action='store_true',  
              help = 'DEBUG: skip hisat')
    parser.add_argument('--skip_bowtie', action='store_true',  
              help = 'DEBUG: skip bowtie')
    parser.add_argument('--skip_post_process_bam', action='store_true',  
              help = 'DEBUG: skip combining BAM, multimap reassignment and BED file conversion')
    parser.add_argument('--skip_remap', action='store_true',  
              help = 'DEBUG: skip tRNA/rRNA remapping')
    parser.add_argument('--skip_count', action='store_true',  
              help = 'DEBUG: skip counting')
    #parser.add_argument('--hisat2', default='hisat2',  
    #          help = "PATH to Douglas's version of HISAT2, to allow dovetails")
    args = parser.parse_args()
    return args


def main():
    programname = sys.argv[0]
    start = time.time()

    args = getopt()
    process_sample = sample_object(args)

    process_sample.make_result_dir()

    if not args.skip_trim:
        process_sample.trimming()

    if args.univec_index:
        if not args.skip_univec:
            process_sample.RNA_filter(RNA='univec')
        process_sample.trimed1, process_sample.trimed2 = process_sample.filtered_fq1, process_sample.filtered_fq2

    if args.rRNA_mt_index:
        if not args.skip_premap:
            process_sample.RNA_filter(RNA='rRNA_mt')
        process_sample.trimed1, process_sample.trimed2 = process_sample.filtered_fq1, process_sample.filtered_fq2

    if args.smRNA_index:
        if not args.skip_smRNA:
            process_sample.RNA_filter(RNA='smallRNA')
        process_sample.trimed1, process_sample.trimed2 = process_sample.filtered_fq1, process_sample.filtered_fq2


    if not args.skip_hisat:
        process_sample.hisat_map()

    if not args.skip_bowtie:
        process_sample.bowtie_map()

    if not args.skip_post_process_bam:
        process_sample.combined_aligned()

        if args.umi > 0 and not args.count_all:
            process_sample.dedup_bam()

    if not args.skip_remap:
        process_sample.combined_filter()
        process_sample.make_alignment_bed()

    if not args.skip_count:
        process_sample.generate_all_count()
        if args.repeats and args.repeats_index:
            process_sample.generate_repeat_count()


    end = time.time()
    usedTime = end - start
    print('Finished: %s in %.3f hr\n' %(process_sample.samplename ,usedTime/3600), file=sys.stderr)
    return 0

if __name__ == '__main__':
    main()