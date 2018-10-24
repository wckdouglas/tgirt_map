#!/usr/bin/env python

from __future__ import print_function
from Bio import Entrez, SeqIO
import ssl
ssl._create_default_https_context = ssl._create_unverified_context
import sys
import pysam
from xopen import xopen
from get_fa import get_sequences

if len(sys.argv) != 5:
    sys.exit('[usage] python %s <genome_fa> <in_bed> <out_bed_file> <out_fasta>' %sys.argv[0] )
genome = sys.argv[1]
in_bed = sys.argv[2]
out_bed = sys.argv[3]
out_fa = sys.argv[4]

rRNA_gene = ['gi|23898|emb|X12811.1| Human 5S DNA',
    'gi|555853|gb|U13369.1|HSU13369 Human ribosomal DNA complete repeating unit']
Entrez.email = 'wckdouglas@gmail.com'

rRNA_genes_bed = 'gi|23898|emb|X12811.1|\t274\t394\t5S_rRNA\t0\t+\t5S_rRNA\t5S_rRNA\n'\
            'gi|555853|gb|U13369.1|HSU13369\t3657\t5527\t18S_rRNA\t0\t+\t18S_rRNA\t18S_rRNA\n'\
            'gi|555853|gb|U13369.1|HSU13369\t6623\t6779\t5.8S_rRNA\t0\t+\t5.8S_rRNA\t5.8S_rRNA\n'\
            'gi|555853|gb|U13369.1|HSU13369\t7935\t12969\t28S_rRNA\t0\t+\t28S_rRNA\t28S_rRNA'

with xopen(out_fa, 'w') as fa, xopen(out_bed, 'w') as bed:
    for rRNA in rRNA_gene:
        id = rRNA.split('|')[1]
        handle = Entrez.efetch(db="nucleotide", id=id, 
                               rettype="fasta", retmode="text")
        record = handle.read()
        print('>' + rRNA +'\n'+ ''.join(record.split('\n')[1:]), file = fa)
    print(rRNA_genes_bed, file = bed)


    # MT
    genome_fa = pysam.Fastafile(genome)
    in_bed = open(in_bed)

    for gene in in_bed:
        fields = gene.strip().split('\t')
        gname = fields[3]
        gtype = fields[6]
        if gname.startswith('MT-') and not gtype.endswith('tRNA'):
            get_sequences(genome_fa, gene, bed, fa)

