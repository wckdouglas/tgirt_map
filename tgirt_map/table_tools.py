#!/usr/bin/env python

from __future__ import print_function
import pandas as pd
import re
import numpy as np
import os
import glob
from functools import partial
from multiprocessing import Pool
import sys


ncRNA = ["sense_intronic","3prime_overlapping_ncRNA",'processed_transcript','TEC',
        'sense_overlapping','Other_lncRNA', 'macro_lncRNA','non_coding','known_ncrna', 'LRG_gene',
        'lincRNA','bidirectional_promoter_lncRNA', 'ribozyme','3prime_overlapping_ncrna']
smncRNA = ['misc_RNA','snRNA','piRNA','scaRNA','sRNA','scRNA','Y_RNA','Y RNA', 'vt RNA','7SK','7SL']
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

def readDF(count_file_name):
    df = pd.read_table(count_file_name, header=None)  \
        .pipe(lambda d: d[[3,6,7,8]])
    df.columns = ['name','type','id','count']
    return df


def assign_gene_type(gid):
    gt = ''
    if re.search('^TR[A-Z]|MT-T[A-Z]',gid):
        gt = 'tRNA'
    elif re.search('^RN7SK', gid):
        gt = '7SK'
    elif re.search('^RN7SL', gid):
        gt = '7SL'
    elif re.search('^RNY', gid):
        gt = 'Y RNA'
    elif re.search('hsa-mir|hsa-let', gid):
        gt = 'miRNA'
    elif re.search('^VT|vault', gid):
        gt = 'vt RNA'
    elif re.search('^MT-RNR', gid):
        gt = 'Mt-rRNA'
    elif re.search('^MT-', gid):
        gt = 'Mt-protein_coding'
    elif re.search('rRNA', gid):
        gt = 'rRNA'
    elif re.search('rDNA', gid):
        gt = 'rDNA'

    assert gt != '', gid
    return gt

def read_direct_counts(filename):
    df = pd.read_table(filename, names = ['id','count']) \
        .assign(type = lambda d: d['id'].map(assign_gene_type))\
        .assign(name = lambda d: d['id']) \
        .assign(name = lambda d: np.where((d.type == 'tRNA') & (~d.name.str.contains('MT')),
                                        d.name.str.extract('(TR[A-Za-z]+-[ACTGN]{3})', expand=False),
                                        d.name))
    return df


def readSample(project_path, sample_id):    
    print('Running %s' %sample_id)
    count_file_path = project_path + '/Counts/RAW'
    df = readDF(count_file_path + '/' + sample_id + '.counts') 

    smallRNA = project_path + '/' + sample_id + '/smallRNA/aligned.count'
    smallRNA_df = read_direct_counts(smallRNA)

    rRNA_mt = project_path + '/' + sample_id + '/rRNA_mt/aligned.count'
    rRNA_mt_df = read_direct_counts(rRNA_mt)

    df = pd.concat([df, smallRNA_df, rRNA_mt_df],axis=0, sort=True) \
        .assign(sample_name = sample_id.replace('-','_'))  \
        .assign(grouped_type = lambda d: d.type.map(change_gene_type))
    return df


def make_table(project_path):
    count_file_path = project_path + '/Counts/RAW'
    count_files = glob.glob(count_file_path + '/*counts')
    sample_ids = set(map(lambda x: x.split('/')[-1].split('.')[0], count_files))
    dfFunc = partial(readSample, project_path)
    dfs = Pool(12).map(dfFunc, sample_ids)

    df = pd.concat(dfs, axis=0, sort=True) \
        .assign(count = lambda d: d['count'].astype(int)) \
        .pipe(pd.pivot_table,
            index = ['id','grouped_type','type','name'],  
            values = 'count' , 
            columns = ['sample_name'],
            aggfunc = 'sum',
            fill_value=0) \
        .reset_index() 
    tablename = project_path + '/Counts/combined_gene_count.tsv'
    df.to_csv(tablename, sep='\t', index=False)
    print('Written %s' %tablename)

