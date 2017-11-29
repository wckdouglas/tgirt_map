#!/usr/bin/env python

import pandas as pd
import MySQLdb
import sys

def blocksize(a,b):
    starts = a.split(',')[:-1]
    ends = b.split(',')[:-1]
    sizes = map(lambda x,y: str(long(y)-long(x)), starts,ends)
    return ','.join(sizes) + ','

def adjust_start(a,b):
    exon_starts = b.split(',')[:-1]
    adjusted_starts = map(lambda x: str(long(x) - long(a)), exon_starts)
    return ','.join(adjusted_starts) + ','

outpath = '/corral-repl/utexas/2013lambowitz/Ref/hg19/genome'
out_bed = outpath + '/genes.bed'
user='genome' 
host='genome-mysql.soe.ucsc.edu'
ucsc= MySQLdb.connect(user=user,host=host, db='hg19')
#df = pd.read_sql('select chrom, txstart, txEnd, name, strand, geneSymbol, '+\
#                'from knownGene, kgXref where kgXref.kgId=knownGene.name',
#                con=ucsc)
df = pd.read_sql('select * from ensGene', con=ucsc) \
    .assign(rgb = '(255,0,0)') \
    .assign(block_size = lambda d: map(blocksize, d.exonStarts, d.exonEnds))\
    .assign(exonStarts = lambda d: map(adjust_start, d.txStart, d.exonStarts)) \
    .pipe(lambda d: d['chrom,txStart,txEnd,name2,score,strand,txStart,txEnd,rgb,exonCount,block_size,exonStarts'.split(',')])
df.to_csv(sys.stdout, sep='\t', index=False, header=False)
