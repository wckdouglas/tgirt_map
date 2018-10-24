#!/usr/bin/env Rscript

from __future__ import print_function
import pandas as pd
import sys


if len(sys.argv) != 2:
    sys.exit('[usage] python %s <ref path>' %(sys.argv[0]))
bed_path = sys.argv[1]#'/stor/work/Lambowitz/ref/RNASeqConsortium'

gene_bed = bed_path + '/genes.bed'
gene_bed = pd.read_table(gene_bed, 
                         names=['chrom','start','end','gene_name',
                                 'score','strand','bio_type','gene_id'],
                         dtype = {'chrom':'str'}) 


def filter_bed(type_pattern, filename):
    bed_file_name = '%s/%s.bed' %(bed_path, filename)
    if 'rRNA_for_bam_filter' == filename:
        gene_bed[gene_bed.bio_type.str.contains('rRNA|rDNA')] \
            .query('bio_type!="Mt_rRNA"') \
            .to_csv(bed_file_name, header=False, index=False ,sep='\t')
    elif 'genes_no' not in filename:
        gene_bed[gene_bed.bio_type.str.contains(type_pattern)] \
            .to_csv(bed_file_name, header=False, index=False ,sep='\t')
    else:
        gene_bed[~gene_bed.bio_type.str.contains(type_pattern)] \
            .to_csv(bed_file_name, header=False, index=False, sep='\t')
    print('Written: ', bed_file_name)

patterns = ['miRNA|misc_RNA|snoRNA|snRNA',
              'protein_coding',
             'miRNA|misc_RNA|snoRNA|snRNA|tRNA|scRNA',
             'rRNA|rDNA',
             '18S_rRNA|28S_rRNA|5.8S_rRNA|5S_rRNA|miRNA|misc_RNA|rRNA|snoRNA|snRNA|tRNA|scRNA|piRNA',
             'rRNA',
             'miRNA|misc_RNA|rRNA|snoRNA|snRNA|tRNA']
filenames = ['sncRNA_no_tRNA',
               'protein',
              'sncRNA_x_protein',
              'rDNA', 
              'genes_no_sncRNA_rRNA_tRNA',
              'rRNA_for_bam_filter',
              'sncRNA_rRNA_for_bam_filter']
for pattern, filename in zip(patterns, filenames):
    filter_bed(pattern, filename)
