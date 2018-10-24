REF=$HOME/ref
mkdir -p $REF/genome
curl  http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chrX.fa.gz \
    | zcat \
    > $REF/hg19_genome.fa

bowtie2-build $REF/hg19_genome.fa $REF/hg19_genome
hisat2-build $REF/hg19_genome.fa $REF/hg19_genome


tgirt_count.py map \
    -1 test_map/test.1.fq.gz \
    -2 test_map/test.2.fq.gz \
    --outdir test_map/test_result \
    --samplename test --univec UniVec_Core \
    --hisat_index $REF/genome/hg19_genome --bowtie2_index $REF/genome/hg19_genome \
    --bedpath $REF/genes --splicesite $REF/genome/splicesites.tsv \
    --rRNA_mt_index $REF/genes/rRNA_mt --smRNA_index $REF/genes/smallRNA -p 24 \
    --trim_aggressive --repeats rmsk.bed.gz --umi 6 \
    --repeats_index all_rmsk_From_bed \
    --snakemake 
