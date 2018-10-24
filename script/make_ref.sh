#!/bin/bash

REF_PATH=${REF}/hg19/genome
GENOME_PREFIX=hg19_genome
mkdir -p  $REF_PATH/all_fasta

HG_VERSION=hg19
if echo $HG_VERSION | grep -v -q 'hg19\|hg38' 
then
	echo HG version has to be hg38 or hg19
	exit	
fi


tRNADB=Hsapi$(echo $HG_VERSION | cut -c3-)

#LINKS
CHROM=http://hgdownload.cse.ucsc.edu/goldenPath/$HG_VERSION/bigZips/chromFa.tar.gz
GENCODE_GTF=ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gz
tRNA_LINK=http://gtrnadb.ucsc.edu/genomes/eukaryota/$tRNADB/$HG_VERSION-tRNAs.tar.gz


#Download rRNA
python get_rRNA_fa.py > $REF_PATH/rRNA.fa
echo 'gi|23898|emb|X12811.1|  274     394     5S_rRNA 0       +       5S_rRNA 5S_rRNA
gi|555853|gb|U13369.1|HSU13369  3657    5527    18S_rRNA        0       +       18S_rRNA        18S_rRNA
gi|555853|gb|U13369.1|HSU13369  6623    6779    5.8S_rRNA       0       +       5.8S_rRNA       5.8S_rRNA
gi|555853|gb|U13369.1|HSU13369  7935    12969   28S_rRNA        0       +       28S_rRNA        28S_rRNA' \
		| awk '{print $1,$2,$3,$4,$5,$6,"rDNA",$8}' OFS='\t' \
			> $REF_PATH/rRNA.bed
echo 'Made rRNA'

# downolad reference
curl -o $REF_PATH/chroms.tar.gz $CHROM
tar zxvf $REF_PATH/chroms.tar.gz --directory $REF_PATH/all_fasta
cat $REF_PATH/all_fasta/*.fa $REF_PATH/rRNA.fa > $REF_PATH/${GENOME_PREFIX}.fa

#make index
hisat2-build $REF_PATH/${GENOME_PREFIX}.fa $REF_PATH/${GENOME_PREFIX}
bowtie2-build $REF_PATH/${GENOME_PREFIX}.fa $REF_PATH/${GENOME_PREFIX}
bwa index $REF_PATH/${GENOME_PREFIX}.fa
samtools faidx $REF_PATH/${GENOME_PREFIX}.fa
echo 'Made genome'

#downolaod annotation
curl  $GENCODE_GTF | zcat > $REF_PATH/gencode_genes.gtf
hisat2_extract_splice_sites.py $REF_PATH/gencode_genes.gtf > $REF_PATH/splicesite.tsv
echo 'Made splice sites'

#make genome file
python ~/ngs_qc_plot/make_genome.py -f $REF_PATH/${GENOME_PREFIX}.fa -o genome > $REF_PATH/${GENOME_PREFIX}.genome

#make gene bed
Rscript get_ensemble.R \
	| cat - $REF_PATH/rRNA.bed > $REF_PATH/genes.bed
python get_genes.py > $REF_PATH/genes.bed12
echo 'Downloaded genes.bed'

#make tRNA
tRNA_PATH=$REF_PATH/tRNA
mkdir -p $tRNA_PATH

FILENAME=$(basename $tRNA_LINK)
curl -o $tRNA_PATH/$FILENAME $tRNA_LINK
tar zxvf $tRNA_PATH/$FILENAME --directory $tRNA_PATH
cat $tRNA_PATH/hg19-tRNAs.bed \
	    | cut -f1-6 \
		| awk '{print $1, $2,$3,$4,$5,$6,"tRNA",$4}' OFS='\t' \
		> $tRNA_PATH/tRNA.bed
cat $tRNA_PATH/tRNA.bed >>  $REF_PATH/genes.bed
python make_tRNA_fasta.py $tRNA_PATH > $tRNA_PATH/nucleo_tRNA.fa
cat $REF_PATH/genes.bed \
	| grep 'Mt_tRNA' \
	| bedtools getfasta  -fi $REF_PATH/hg19_genome.fa -bed - -s -name -tab \
	| tr ':' '\t' \
	| awk '{printf ">%s\n%s\n",$1,$NF}' \
	> $tRNA_PATH/mt_tRNA.fa
cat $tRNA_PATH/mt_tRNA.fa $tRNA_PATH/nucleo_tRNA.fa > $REF_PATH/tRNA.fa
mv $tRNA_PATH/tRNA.bed $REF_PATH
echo 'Finished making tRNA'

# make index
cat  $REF_PATH/tRNA.fa $REF_PATH/rRNA.fa \
		> $REF_PATH/tRNA_rRNA.fa
bowtie2-build $REF_PATH/tRNA_rRNA.fa $REF_PATH/tRNA_rRNA
bowtie2-build $REF_PATH/tRNA.fa $REF_PATH/tRNA
bowtie2-build $REF_PATH/rRNA.fa $REF_PATH/rRNA
echo 'Indexed tRNA/rRNA'

#make gene count bed
python split_bed_for_count.py $REF_PATH

