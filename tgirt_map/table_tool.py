#!/usr/bin/env python

import re


ncRNA = ["sense_intronic","3prime_overlapping_ncRNA",'processed_transcript','TEC',
        'sense_overlapping','Other_lncRNA', 'macro_lncRNA','non_coding','known_ncrna', 'LRG_gene',
        'lincRNA','bidirectional_promoter_lncRNA', 'ribozyme','3prime_overlapping_ncrna']
smncRNA = ['misc_RNA','snRNA','piRNA','scaRNA','sRNA','scRNA']
large_rRNA = ['28S_rRNA','18S_rRNA']
small_rRNA = ['rRNA','5S_rRNA','58S_rRNA','5.8S_rRNA', 'rDNA']
protein_coding = ['protein_coding','TR','IG']
def change_gene_type(x):
    type = ''
    if x in ncRNA or re.search('pseudo',x):
        type = 'Other ncRNA'
    elif re.search('TR|IG|protein',x):
        type = 'Protein coding'
    elif x.startswith('Mt_'):
        type = 'Mt'
    elif x == 'tRNA':
        type = 'tRNA'
    elif x in small_rRNA or x in large_rRNA:
        type = 'rRNA'
    elif x in smncRNA:
        type = 'Other sncRNA'
    elif x =='antisense':
        type = 'Antisense'
    else:
        type = x
    return type

