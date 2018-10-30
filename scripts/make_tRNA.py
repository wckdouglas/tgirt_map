#!/usr/bin/env python

import fileinput
import sys
from itertools import islice
import re
import pandas as pd
import numpy as np



amino_acid = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
          'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
          'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
          'SEC': 'U', 'SUP':'SUP', 'UNDET':'UND', 'IMET': 'IMET',
          'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

def read_tRNA(tRNA_align):
    ln_count = 1
    for line in tRNA_align:
        if ln_count == 1:
            descriptor = line.strip()
        elif ln_count == 2:
            tRNA_info = line.strip()
        elif ln_count == 4:
            seq = line.strip()
        elif ln_count == 6:
            yield(descriptor, tRNA_info, seq)
            ln_count = 0
        ln_count += 1
        
def filter_file(infile):
    with open(infile, 'r') as in_f:
        lines = in_f.readlines()
        lines = filter(lambda x:  not 'Possible intron' in x and not 'HMM' in x, lines)
    return lines


def make_mature(sequence):
    sequence = re.sub('[actgn]','', sequence)
    sequence = sequence + 'CCAA'
    return sequence


def make_tRNA_info(filtered_lines):
    for tRNA_record in read_tRNA(filtered_lines):
        descriptor = tRNA_record[0]
        chrom = descriptor.split('.')[0]
        gene_range = re.search('\(([0-9]+-[0-9]+)\)',descriptor).group(1)
        start, end = gene_range.split('-')

        tRNA_info = tRNA_record[1]
        AA = tRNA_info.split(' ')[1].split('\t')[0]
        score = tRNA_info.split(' ')[-1]

        tRNA_detail = re.search('Anticodon: ([ACTGN]{3}) at ([0-9]+-[0-9]+)', tRNA_info)
        anticodon = tRNA_detail.group(1)
        anticodon_range = tRNA_detail.group(2)

        seq = tRNA_record[2].split(' ')[-1]
        if start > end:
            start, end = end, start
            strand = '-'
        else:
            strand = '+'
        
        AA_code = amino_acid[AA.upper()]
        name = 'TR%s-%s' %(AA_code, anticodon)
        mature_seq = make_mature(seq)
        yield(chrom, start, end, AA, score, mature_seq, anticodon, strand, name)


def make_tRNA_name(tRNA_gene_df):
    tRNA_gene_df = tRNA_gene_df\
        .assign(code = lambda d: d['score']\
                .astype('category')\
                .cat\
                .codes[::-1]) \
        .assign(dummy = 1) \
        .assign(tRNA_number = lambda d: d\
                                    .groupby('code')['dummy']\
                                    .transform(np.cumsum)) \
        .drop('dummy', axis=1)
    return tRNA_gene_df
    

def make_bed(tRNA_df, bed_file):
    bed = tRNA_df \
        .sort_values(['chrom','start']) \
        .pipe(lambda d: d[['chrom','start','end','name','score',
                            'strand','gene_type','gene_id']])
    bed.to_csv(bed_file, sep='\t', index=False, header=None)
    print('Written %s' %bed_file)


def make_fasta(tRNA_df, fa_file):
    with open(fa_file, 'w') as out_fa:
        for i, row in tRNA_df \
                .groupby(['anticodon','mature_seq'], as_index=False) \
                .agg({'gene_id': lambda x:  np.sort(x)[0]}) \
                .iterrows():
            print('>%s\n%s' %(row['gene_id'],
                                row['mature_seq']),
                    file = out_fa)
    print('Wrtieen %s' %fa_file)



def main():
    if len(sys.argv) != 4:
        sys.exit('[usage] python %s <hg38-tRNAs-detailed.ss> <out_tRNA_fa> <out_bed_file>' %sys.argv[0])

    infile = sys.argv[1]
    bed_file = sys.argv[2]
    fa_file = sys.argv[3]

    filtered_lines = filter_file(infile)
    lines = make_tRNA_info(filtered_lines)
    df = pd.DataFrame(lines,
                columns = 'chrom|start|end|AA|score|mature_seq'\
                        '|anticodon|strand|name'.split('|')) \
        .groupby(['name']) \
        .apply(make_tRNA_name) \
        .assign(gene_id = lambda d: d.name + \
                                d.code.astype(str) + \
                                '-' + \
                                d.tRNA_number.astype(str)) \
        .assign(gene_type = 'tRNA')
    make_bed(df, bed_file)
    make_fasta(df, fa_file)




if __name__ == '__main__':
    main()
