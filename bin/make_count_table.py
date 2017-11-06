#!/usr/bin/env python

import pandas as pd
import re
import numpy as np
import os
import glob
from functools import partial
from multiprocessing import Pool
import sys



ncRNA = ["sense_intronic","3prime_overlapping_ncRNA",'processed_transcript',
        'sense_overlapping','Other_lncRNA', 'macro_lncRNA','non_coding','known_ncrna', 'LRG_gene',
        'lincRNA','bidirectional_promoter_lncRNA', 'ribozyme','3prime_overlapping_ncrna']
smncRNA = ['misc_RNA','snRNA','piRNA','scaRNA','sRNA','scRNA']
large_rRNA = ['28S_rRNA','18S_rRNA']
small_rRNA = ['rRNA','5S_rRNA','58S_rRNA','5.8S_rRNA']
protein_coding = ['protein_coding','TR','IG']
def changeType(x):
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

def readDF(count_file_name):
    df = pd.read_table(count_file_name, header=None)  \
        .pipe(lambda d: d[[3,6,7,8]])
    df.columns = ['name','type','id','count']
    return df

def read_tRNA(count_file_name):
    df = pd.read_table(count_file_name, names= ['id','count'])  \
            .assign(name = lambda d: d['id']) \
            .assign(type = 'tRNA')
    return df


def readSample(count_file_path, tRNA_count_path, sample_id):    
    print 'Running %s' %sample_id
    df = readDF(count_file_path + '/' + sample_id + '.counts') 
    tRNA_df = read_tRNA(tRNA_count_path + '/' + sample_id + '.tRNA')
    df = pd.concat([df, tRNA_df],axis=0) \
        .assign(sample_name = sample_id.replace('-','_'))  \
        .assign(type = lambda d: d.type.map(changeType))
    return df

def main():
    if len(sys.argv) != 2:
        sys.exit('[usage] python %s <count_file_path>')
    count_file_path = sys.argv[1] 
    count_files = glob.glob(count_file_path + '/Counts/RAW/*counts')
    sample_ids = set(map(lambda x: x.split('/')[-1].split('.')[0], count_files))
    dfFunc = partial(readSample, count_file_path, tRNA_count_path)
    dfs = Pool(12).map(dfFunc, sample_ids)
    df = pd.concat(dfs, axis=0) \
        .query('sample_name != "try"') \
        .assign(count = lambda d: d['count'].astype(int)) \
        .pipe(pd.pivot_table,index = ['id','type','name'],  
            values = 'count' , columns = ['sample_name'],
            fill_value=0) \
        .reset_index() \
        .pipe(lambda d: d[~d.type.str.contains('ERCC')])
    df.iloc[:,3:] = df.iloc[:,3:].astype(int)
    tablename = count_file_path + '/Counts/combined_gene_count.tsv'
    df.to_csv(tablename, sep='\t', index=False)
    print 'Written %s' %tablename


if __name__ == '__main__':
    main()
