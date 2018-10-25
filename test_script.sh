WORKING_DIR=$(pwd)
REF=$WORKING_DIR/test_map/ref
cd script
snakemake -s make_ref.smk -p -j 4 --config test=1 path=$REF
cd $WORKING_DIR


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
