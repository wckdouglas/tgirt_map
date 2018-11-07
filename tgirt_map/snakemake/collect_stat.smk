import glob
import re
import os
import pandas as pd
from tgirt_map.table_tools import read_flag_stat


PATH = config['PATH']
FOLDERS = glob.glob(PATH + '/*')
regex = config['regex']
SAMPLENAMES = map(os.path.basename,FOLDERS)
SAMPLENAMES = list(filter(lambda x: re.search(regex, x), SAMPLENAMES))
print('Running: %s' %'\n'.join(SAMPLENAMES))
FEATURES = ['bowtie','hisat','rRNA_mt','smallRNA', 'UniVec']
MAPPERS = ['bowtie','hisat']

# define varible
STAT_TABLE = PATH + '/mapping_stat.tsv'
SAMPLE_FOLDER_TEMPLATE = PATH + '/{SAMPLENAME}'
BAM_TEMPLATE = SAMPLE_FOLDER_TEMPLATE + '/{FEATURE}/aligned.bam'
FLAGSTAT_TEMPLATE = BAM_TEMPLATE.replace('.bam','.flagstat')
UNIQUE_BAM_TEMPLATE = SAMPLE_FOLDER_TEMPLATE + '/{MAPPER}/{MAPPER}.unique.bam'
UNIQUE_FLAGSTAT_TEMPLATE = UNIQUE_BAM_TEMPLATE.replace('.bam','.flagstat')
TRIMMED_FQ = SAMPLE_FOLDER_TEMPLATE + '/Trimmed/trim.1.fq.gz'
TRIMMED_TXT = TRIMMED_FQ.replace('.1.fq.gz','.count')
SAMPLE_STAT_TABLE = SAMPLE_FOLDER_TEMPLATE + '/stat_table.tsv'

rule all:
    input:
        STAT_TABLE

rule combine_table:
    input:
        TABLES = expand(SAMPLE_STAT_TABLE,
            SAMPLENAME = SAMPLENAMES)
    
    output:
        TABLE = STAT_TABLE

    run:
        def read_tab(fn):
            return pd.read_table(fn) \
                .assign(samplename = os.path.basename(os.path.dirname(fn)))
        pd.concat(map(read_tab, input['TABLES']))\
            .pipe(pd.pivot_table, index=['samplename'],
                columns = 'col_name',
                values = 'value') \
            .loc[:, map_cat] \
            .reset_index()\
            .assign(mapping_rate = lambda d: d.filter(['rRNA_mt','smallRNA','HISAT2','BOWTIE2']).sum(axis=1)/(d[map_cat[0]] - d[map_cat[1]]))\
            .to_csv(output['TABLE'], sep='\t', index=False)

    

rule make_table:
    input:
        FLAGSTATS = expand(FLAGSTAT_TEMPLATE \
                    .replace('{SAMPLENAME}','{{SAMPLENAME}}'), 
                FEATURE = FEATURES),
        UNIQUE_FLAGSTATS = expand(UNIQUE_FLAGSTAT_TEMPLATE\
                    .replace('{SAMPLENAME}','{{SAMPLENAME}}'),
                MAPPER = MAPPERS),
        TRIMMED = expand(TRIMMED_FQ \ 
                    .replace('{SAMPLENAME}','{{SAMPLENAME}}'))
    

    output:
        TABLE = SAMPLE_STAT_TABLE

    run:
        FLAGSTATS = input['FLAGSTATS']
        FLAGSTATS.extend(input['UNIQUE_FLAGSTATS'])
        fs_df = []
        #for FLAGSTAT in input['FLAGSTATS']:
        for FLAGSTAT in FLAGSTATS:
            stat_df = read_flag_stat(FLAGSTAT)
            in_read = stat_df['read1'] - stat_df['secondary']/2
            mapped_read = (stat_df['mapped'] - stat_df['secondary'] - stat_df['supplementary'])/2
            fs_df.append((in_read, mapped_read, FLAGSTAT))
        
        fs_df = pd.DataFrame(fs_df, columns = ['inread', 'mapped_read', 'flagstat_name'])\
            .pipe(pd.melt, id_vars = ['flagstat_name']) \
            .assign(col_name = lambda d: list(map(rename_variable, d.flagstat_name, d.variable))) \
            .query('col_name != ""') \
            .assign(col_name = lambda d: pd.Categorical(d.col_name, categories=map_cat, ordered=True))\
            .filter(['col_name','value']) \
            .sort_values('col_name')  
        
        fs_df.to_csv(output['TABLE'], index=False, sep='\t')


rule make_unique_stat:
    input:
        UNIQUE_BAM_TEMPLATE

    output:
        UNIQUE_FLAGSTAT_TEMPLATE

    shell:
        'samtools flagstat {input} > {output}'

    
rule make_stat:
    input:
        BAM_TEMPLATE

    output:
        FLAGSTAT_TEMPLATE

    shell:
        'samtools flagstat {input} > {output}'

rule trimmed_stat:
    input:
        TRIMMED_FQ

    output:
        TRIMMED_TXT
    
    shell:
        "zgrep -c '^+$' {input} > {output}"


map_cat = ['trimmed', 'UniVec_contam', 'rRNA_mt', 'smallRNA', 
        'HISAT2', 'HISAT2_unique', 'BOWTIE2','BOWTIE2_unique']
def rename_variable(fs_name, var_name):
    variable = ''
    if 'UniVec' in fs_name and var_name == 'inread':
        variable = map_cat[0]
    elif 'UniVec' in fs_name and var_name == 'mapped_read':
        variable = map_cat[1]
    elif 'rRNA_mt' in fs_name and var_name == 'mapped_read':
        variable = map_cat[2]
    elif 'smallRNA' in fs_name and var_name == 'mapped_read':
        variable = map_cat[3]
    elif 'hisat/aligned' in fs_name and var_name == 'mapped_read':
        variable = map_cat[4]
    elif 'hisat/hisat.unique' in fs_name and var_name == 'mapped_read':
        variable = map_cat[5]
    elif 'bowtie/aligned' in fs_name and var_name == 'mapped_read':
        variable = map_cat[6]
    elif 'bowtie/bowtie.unique' in fs_name and var_name == 'mapped_read':
        variable = map_cat[7]
    return variable
    