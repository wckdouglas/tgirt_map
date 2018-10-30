#!/usr/bin/env python

import sys
from Bio import SeqIO
import re


tRNA_search = re.compile('^TR[A-Z]+-[AGCTN]{3}')
def categorize_gene(seq_name):
    gtype = ''

    if seq_name.startswith('MT-T'):
        gtype = 'tRNA'
    elif tRNA_search.search(seq_name):
        gtype = 'tRNA'
    elif 'hsa-mir' in seq_name or 'hsa-let' in seq_name:
        gtype = 'miRNA'
    elif seq_name.startswith('VTRNA'):
        gtype = 'vaultRNA'
    elif '7SL' in seq_name:
        gtype = '7SL'
    elif '7SK' in seq_name:
        gtype = '7SK'
    elif 'RNY' in seq_name:
        gtype = 'Y-RNA'
    return gtype




if len(sys.argv)!=2:
    sys.exit('[usage] python %s <smallRNA.fasta>' %sys.argv[0])

fasta_file = sys.argv[1]

bedline =  '{seq_name}\t0\t{seq_length}\t{gname}\t0\t+\t{seq_type}\t{seq_name}'
for record in SeqIO.parse(fasta_file, 'fasta'):
    gtype = categorize_gene(record.id)
    if gtype=='':
        sys.exit(record.id)

    if gtype == 'miRNA':
        gene_name = record.id.replace('hsa-','')

    else:
        gene_name = record.id

    outline = bedline.format(seq_name = record.id,
                   seq_length = len(record.seq),
                   gname = gene_name,
                   seq_type = gtype)
    print(outline)
