REF=$HOME/ref
mkdir -p $REF/genome
for CHROM in chrX chrM
do
    curl  http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/${CHROM}.fa.gz \
        | zcat 
done > $REF/genome/hg19_genome.fa

bowtie2-build $REF/genome/hg19_genome.fa $REF/genome/hg19_genome
hisat2-build $REF/genome/hg19_genome.fa $REF/genome/hg19_genome

cd script
bash make_ref.sh $REF
cd ../


tgirt_count.py map \
    -1 test_map/test.1.fq.gz \
    -2 test_map/test.2.fq.gz \
    --outdir test_map/test_result \
    --samplename test --univec $REF/genes/UniVec_core \
    --hisat_index $REF/genome/hg19_genome --bowtie2_index $REF/genome/hg19_genome \
    --bedpath $REF/genes --splicesite $REF/genome/splicesites.tsv \
    --rRNA_mt_index $REF/genes/rRNA_mt --smRNA_index $REF/genes/smallRNA -p 3 \
    --trim_aggressive --repeats rmsk.bed.gz --umi 6 \
    --repeats_index all_rmsk_From_bed \
    --snakemake 
