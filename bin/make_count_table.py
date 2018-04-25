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
from tgirt_map.table_tool import change_gene_type


def readDF(count_file_name):
    df = pd.read_table(count_file_name, header=None)  \
        .pipe(lambda d: d[[3,6,7,8]])
    df.columns = ['name','type','id','count']
    return df

def read_tRNA(count_file_name):
    df = pd.read_table(count_file_name, names= ['id','count'])  \
        .assign(type = 'tRNA') \
        .assign(name = lambda d: np.where(d.id.str.contains('tRNA'),
                                        d.id.str.extract('(tRNA-[A-Za-z]{3,5}-[ACTGN]{3})', expand=False),
                                        d.id.str.extract('(MT-.*)', expand=False))) \
        .assign(name = lambda d: d.name.str.replace('^Homo_sapiens_|-[0-9]+$',''))\
        .assign(name = lambda d: d.name.str.replace('-[0-9]+$','')) \
        .assign(name = lambda d: d.name.str.replace('-[0-9]+$',''))  \
        .assign(name = lambda d: np.where(pd.isnull(d.name), d.id, d.name))
    return df


def readSample(count_file_path, tRNA_count_path, sample_id):    
    print('Running %s' %sample_id)
    df = readDF(count_file_path + '/' + sample_id + '.counts') 
    tRNA_df = read_tRNA(tRNA_count_path + '/' + sample_id + '.tRNA')
    df = pd.concat([df, tRNA_df],axis=0) \
        .assign(sample_name = sample_id.replace('-','_'))  \
        .assign(grouped_type = lambda d: d.type.map(change_gene_type))
    return df

def main():
    if len(sys.argv) != 2:
        sys.exit('[usage] python %s <output_folder_path>' %(sys.argv[0]))
    project_path = sys.argv[1] 
    count_file_path = project_path + '/Counts/RAW'
    tRNA_count_path = project_path + '/Counts/tRNA_RAW'
    count_files = glob.glob(count_file_path + '/*counts')
    sample_ids = set(map(lambda x: x.split('/')[-1].split('.')[0], count_files))
    dfFunc = partial(readSample, count_file_path, tRNA_count_path)
    dfs = Pool(12).map(dfFunc, sample_ids)
    df = pd.concat(dfs, axis=0) \
        .assign(count = lambda d: d['count'].astype(int)) \
        .pipe(pd.pivot_table,
            index = ['id','grouped_type','type','name'],  
            values = 'count' , 
            columns = ['sample_name'],
            fill_value=0) \
        .reset_index() 
    tablename = project_path + '/Counts/combined_gene_count.tsv'
    df.to_csv(tablename, sep='\t', index=False)
    print('Written %s' %tablename)


if __name__ == '__main__':
    main()
