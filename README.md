

# tgirt_map # [![Build Status](https://travis-ci.org/wckdouglas/tgirt_map.svg?branch=master)](https://travis-ci.org/wckdouglas/tgirt_map)

Pipeline for mapping and counting TGIRT RNA-seq

### Installation ###

```
pip install git+https://github.com/wckdouglas/tgirt_map.git
```


### Usage ###

```
usage: tgirt_count.py [-h] -1 FASTQ1 -2 FASTQ2 -o OUTDIR -x HISAT_INDEX -y
                      BOWTIE2_INDEX -b BEDPATH -s SPLICESITE -t TRNAINDEX -r
                      RRNAINDEX -e RRNA_TRNA_INDEX [-p THREADS] [--TTN]
                      [--umi UMI] [--dry] [--skip_trim] [--skip_premap]
                      [--skip_hisat] [--skip_bowtie] [--skip_post_process_bam]
                      [--skip_remap] [--skip_count]

Pipeline for mapping and counting for TGIRT-seq paired end data

optional arguments:
  -h, --help            show this help message and exit
  -1 FASTQ1, --fastq1 FASTQ1
                        pairedEnd fastq file (read1)
  -2 FASTQ2, --fastq2 FASTQ2
                        pairedEnd fastq file (read2)
  -o OUTDIR, --outdir OUTDIR
                        result directory that all resulting/intermediate files
                        will be stored will create 1. $resultpath/trimmed 2.
                        $resultpath/hisat 3. $resultpath/bowtie2 4.
                        $resultpath/mergeBam (all useful result files)
  -x HISAT_INDEX, --hisat_index HISAT_INDEX
                        hisat2 index
  -y BOWTIE2_INDEX, --bowtie2_index BOWTIE2_INDEX
                        bowtie2 index
  -b BEDPATH, --bedpath BEDPATH
                        bed folder for gene counting
  -s SPLICESITE, --splicesite SPLICESITE
                        splice site file generated by hisat
  -t TRNAINDEX, --tRNAindex TRNAINDEX
                        bowtie2 index for tRNA, for better tRNA counting
  -r RRNAINDEX, --rRNAindex RRNAINDEX
                        bowtie2 index for rRNA, for better rRNA counting
  -e RRNA_TRNA_INDEX, --rRNA_tRNA_index RRNA_TRNA_INDEX
                        bowtie2 index for rRNA and tRNA combined
  -p THREADS, --threads THREADS
                        number of cores to be used for the pipeline
                        (default:1)
  --TTN                 used TTN primer
  --umi UMI             Number of UMI bases from 5' of R1 (default = 0)
  --dry                 Dry run
  --skip_trim           DEBUG: skip trimming
  --skip_premap         DEBUG: skip premapping tRNA and rRNA
  --skip_hisat          DEBUG: skip hisat
  --skip_bowtie         DEBUG: skip bowtie
  --skip_post_process_bam
                        DEBUG: skip combining BAM, multimap reassignment and
                        BED file conversion
  --skip_remap          DEBUG: skip tRNA/rRNA remapping
  --skip_count          DEBUG: skip counting
```


### Required pacakge ###
It is necessary to install [sequencing_tools](https://wckdouglas.github.io/sequencing_tools) to run the pipeline correctly as the pipeline extensively used some of the tools.

### building reference ###

Coming soon....

