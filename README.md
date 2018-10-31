# TGIRT-map # 

[![Build Status](https://travis-ci.org/wckdouglas/tgirt_map.svg?branch=master)](https://travis-ci.org/wckdouglas/tgirt_map)


Pipeline for mapping and counting TGIRT RNA-seq


## Start ##

### Installation ###

Using conda to create an environment

```
conda config --add channels r
conda config --add channels bioconda

conda create -q -n tgirt_map python=python3.6 \  
        pandas biopython cython numpy networkx seaborn pyBigwig six pysam ujson \
        hisat2 bowtie2=2.2.5 bedtools samtools snakemake bedtools atropos seqkit \
        xopen ucsc-gtftogenepred ucsc-genePredToBed bedops xopen 
source activate tgirt_map
```

Download TGIRT-map
```
git clone git@github.com:wckdouglas/tgirt_map.git
cd tgirt_map
pip install .

# install sequencing tools
cd ../
git clone https://github.com/wckdouglas/sequencing_tools.git
cd sequencing_tools
pip install .
```

### prepare references ###

```
cd scripts
snakemake -s make_ref.smk -p -j2 
```


### Required pacakge ###
It is necessary to install [sequencing_tools](https://wckdouglas.github.io/sequencing_tools) to run the pipeline correctly as the pipeline extensively used some of the tools.
