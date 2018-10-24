REF_PATH=$REF/hg19
ANNOTATION_PATH=$REF_PATH/new_genes
GENOME_PATH=$REF_PATH/genome
GTF_LINK=ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/GRCh37_mapping/gencode.v28lift37.annotation.gtf.gz
tRNA_REF=http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi19/hg19-tRNAs.tar.gz
piRNA=http://www.regulatoryrna.org/database/piRNA/download/archive/v1.0/bed/piR_hg19_v1.0.bed.gz
MIR_LINK=ftp://mirbase.org/pub/mirbase/CURRENT/hairpin_high_conf.fa.gz

#annotationes
curl $GTF_LINK |zcat > $ANNOTATION_PATH/genes.gtf
cat $ANNOTATION_PATH/genes.gtf \
    | grep 'protein_coding' --color=no\
    | gtfToGenePred /dev/stdin /dev/stdout \
    | genePredToBed > $ANNOTATION_PATH/protein_coding.bed12
cat $ANNOTATION_PATH/genes.gtf \
    | grep 'protein_coding' --color=no \
    | awk '$3=="exon"' \
    | gtf2bed \
    | sort -k1,1 -k2,2n -k3,3n -k6,6 -u \
    > $ANNOTATION_PATH/exons.bed
hisat2_extract_splice_sites.py $ANNOTATION_PATH/genes.gtf > $ANNOTATION_PATH/splicesites.tsv
python gtf_to_bed.py $ANNOTATION_PATH/genes.gtf > $ANNOTATION_PATH/genes.bed


#piRNA
curl $piRNA \
    | zcat \
    | sort -k1,1 -k2,2n -k3,3n \
    | awk '{print $0, "piRNA","piRNA"}' OFS='\t' \
    | bgzip \
    > $ANNOTATION_PATH/piRNA.bed.gz
zcat $ANNOTATION_PATH/piRNA.bed.gz >> $ANNOTATION_PATH/genes.bed

#tRNA
curl $tRNA_REF > $ANNOTATION_PATH/tRNA.tar.gz
mkdir -p $ANNOTATION_PATH/tRNA
tar zxvf $ANNOTATION_PATH/tRNA.tar.gz --directory $ANNOTATION_PATH/tRNA
#python make_tRNA.py \
#    $ANNOTATION_PATH/tRNA/hg19-tRNAs-detailed.ss \
#    $ANNOTATION_PATH/tRNA.bed \
#    $ANNOTATION_PATH/tRNA/nucleo_tRNA.fa
seqkit  rmdup -s  $ANNOTATION_PATH/tRNA/hg19-mature-tRNAs.fa  \
    | python process_mature_tRNA.py \
    > $ANNOTATION_PATH/tRNA/nucleo_tRNA.fa
cat $ANNOTATION_PATH/tRNA.bed |cut -f1-8 >> $ANNOTATION_PATH/genes.bed
cat $ANNOTATION_PATH/genes.bed \
    | grep 'Mt_tRNA' \
    | bedtools getfasta  -fi $GENOME_PATH/hg19_genome.fa -bed - -s -name -tab \
    | tr ':' '\t' \
    | awk '{printf ">%s\n%s\n",$1,$NF}' \
    | sed 's/(-)//g' | sed 's/(+)//g' \
    > $ANNOTATION_PATH/tRNA/mt_tRNA.fa
cat $ANNOTATION_PATH/tRNA/mt_tRNA.fa $ANNOTATION_PATH/tRNA/nucleo_tRNA.fa > $ANNOTATION_PATH/tRNA.fa



#make rRNA
python get_rRNA.py $GENOME_PATH/hg19_genome.fa $ANNOTATION_PATH/genes.bed \
            $ANNOTATION_PATH/rRNA_mt.bed $ANNOTATION_PATH/rRNA_mt.fa 
cat $ANNOTATION_PATH/rRNA_mt.bed | awk '$7!~/Mt/' >> $ANNOTATION_PATH/genes.bed

python split_bed_for_count.py $ANNOTATION_PATH


#make tRNA filter
cat $ANNOTATION_PATH/tRNA.bed $ANNOTATION_PATH/rmsk_tRNA.bed $REF_PATH/genome/tRNA.bed \
    | bedtools sort \
    | bedtools merge -s -o first -c 4,5,6,7,8\
    > $ANNOTATION_PATH/tRNA_comprehensive.bed

#yRNA
cat $ANNOTATION_PATH/genes.bed \
    | grep --color=no 'RNY' \
    | awk '$4!~/[pP]/' \
    | python get_fa.py $GENOME_PATH/hg19_genome.fa $ANNOTATION_PATH/yRNA.bed $ANNOTATION_PATH/yRNA.fa

cat $ANNOTATION_PATH/genes.bed \
     | awk '$4~/.*7SK$|7SL[0-9]+$/' \
     | awk '{print $0, $3-$2}' OFS='\t' \
     | awk  '$NF~/299|330|296/' \
     | python get_fa.py $GENOME_PATH/hg19_genome.fa $ANNOTATION_PATH/srp.bed $ANNOTATION_PATH/srp.fa

cat $ANNOTATION_PATH/genes.bed \
    | grep --color=no 'VTRNA' \
    | awk '$4!~/[pP]/' \
    | python get_fa.py $GENOME_PATH/hg19_genome.fa $ANNOTATION_PATH/vaultRNA.bed  $ANNOTATION_PATH/vaultRNA.fa


curl $MIR_LINK \
    | seqkit grep -n -r -p 'Homo sapien' \
    | seqkit replace -s -p 'U' -r 'T' \
    > $ANNOTATION_PATH/miRNA_hairpin.fa 


#make rRNA tRNA
cat $ANNOTATION_PATH/tRNA.fa \
    $ANNOTATION_PATH/yRNA.fa\
    $ANNOTATION_PATH/srp.fa \
    $ANNOTATION_PATH/miRNA_hairpin.fa \
    $ANNOTATION_PATH/vaultRNA.fa  \
    > $ANNOTATION_PATH/smallRNA.fa 
python smallRNA_bed.py $ANNOTATION_PATH/smallRNA.fa > $ANNOTATION_PATH/smallRNA.bed

cat $ANNOTATION_PATH/tRNA.fa $ANNOTATION_PATH/yRNA.fa > $ANNOTATION_PATH/tRNA_yRNA.fa
cat $ANNOTATION_PATH/genes.bed | awk '$4~"RNY|Y_RNA"' > $ANNOTATION_PATH/yRNA.bed
cat $ANNOTATION_PATH/yRNA.bed $ANNOTATION_PATH/tRNA.bed > $ANNOTATION_PATH/tRNA_yRNA.bed

echo made tRNA_rRNA fasta
bowtie2-build $ANNOTATION_PATH/smallRNA.fa $ANNOTATION_PATH/smallRNA
bowtie2-build $ANNOTATION_PATH/rRNA_mt.fa $ANNOTATION_PATH/rRNA_mt
