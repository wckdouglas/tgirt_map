WORKING_DIR=$(pwd)
REF=$WORKING_DIR/test_map/ref
cd script
snakemake -s make_ref.smk -p -j 2 --config test=1 path=$REF
cd $WORKING_DIR


tgirt_count.py map \
    -1 $WORKING_DIR/data/test.1.fq.gz \
    -2 $WORKING_DIR/data/test.2.fq.gz \
    --outdir $WORKING_DIR/test_result \
    --samplename test --univec $REF/genes/UniVec_core \
    --hisat_index $REF/genome/hg19_genome --bowtie2_index $REF/genome/hg19_genome \
    --bedpath $REF/genes --splicesite $REF/genome/splicesites.tsv \
    --rRNA_mt_index $REF/genes/rRNA_mt --smRNA_index $REF/genes/smallRNA -p 3 \
    --trim_aggressive --repeats $REF/genes/rmsk.bed.gz --umi 6 \
    --repeats_index  $REF/genes/rmsk \
    --snakemake 
