from split_bed_for_count import define_pattern

REF_PATH = config['path']
try:
    test = config['test']
except KeyError:
    test = False

ANNOTATION_PATH = REF_PATH + '/genes'
GENOME_PATH = REF_PATH + '/genome'
GENOME_PREFIX = GENOME_PATH + '/hg19_genome'
GENOME_FA = GENOME_PREFIX + '.fa'
HISAT_GENOME_INDEX = expand(GENOME_PREFIX + '.{NUMBER}.ht2', NUMBER = range(1,9))
BOWTIE2_GENOME_INDEX = expand(GENOME_PREFIX + '.{NUMBER}.bt2', NUMBER = range(1,5))


GENE_BED = ANNOTATION_PATH + '/genes.bed'
GENE_GTF = ANNOTATION_PATH + '/genes.gtf'
SPLICE_SITE = ANNOTATION_PATH + '/splicesites.tsv'
EXONS = ANNOTATION_PATH + '/exons.bed'
BED12 = ANNOTATION_PATH + '/protein_coding.bed12'

UniVec_PREFIX = ANNOTATION_PATH + '/UniVec_core'
UniVec_FA = UniVec_PREFIX + '.fa'
BOWTIE2_UNIVEC_INDEX = expand(UniVec_PREFIX + '.{NUMBER}.bt2', NUMBER = range(1,5))

SMALL_RNA_PREFIX =  ANNOTATION_PATH + '/smallRNA'
SMALL_RNA_FA = SMALL_RNA_PREFIX + '.fa'
SMALL_RNA_BED = SMALL_RNA_PREFIX + '.bed'
BOWTIE2_SMRNA_INDEX = expand(SMALL_RNA_PREFIX + '.{NUMBER}.bt2', NUMBER = range(1,5))

mt_rRNA_PREFIX =  ANNOTATION_PATH + '/rRNA_mt'
mt_rRNA_FA = mt_rRNA_PREFIX + '.fa'
mt_rRNA_BED = mt_rRNA_PREFIX + '.bed'
BOWTIE2_mt_rRNA_INDEX = expand(mt_rRNA_PREFIX + '.{NUMBER}.bt2', NUMBER = range(1,5))

tRNA_FA = ANNOTATION_PATH + '/tRNA.fa'
tRNA_BED = ANNOTATION_PATH + '/tRNA.bed'
tRNA_FOLDER = ANNOTATION_PATH + '/tRNA'
tRNA_ZIP = ANNOTATION_PATH + '/tRNA.tar.gz'
GTRNA_FA = tRNA_FOLDER + '/hg19-mature-tRNAs.fa'
GTRNA_BED = tRNA_FOLDER + '/hg19-tRNAs.bed'
MT_tRNA_FA = tRNA_FOLDER + '/mt_tRNA.fa'
NUCLEO_tRNA_FA = tRNA_FOLDER + '/nucleo_tRNA.fa'

piRNA_BED = ANNOTATION_PATH + '/piRNA.bed.gz'
yRNA_FA = ANNOTATION_PATH + '/yRNA.fa'
yRNA_BED = ANNOTATION_PATH + '/yRNA.bed'
SRP_FA = ANNOTATION_PATH + '/srp.fa'
SRP_BED = ANNOTATION_PATH + '/srp.bed'
miRNA_FA = ANNOTATION_PATH + '/miRNA_hairpin.fa'
vaultRNA_FA = ANNOTATION_PATH + '/vaultRNA.fa'
vaultRNA_BED = ANNOTATION_PATH + '/vaultRNA.bed'
SPLIT_GENE_BED = expand(ANNOTATION_PATH + '/{NAME}.bed', 
                    NAME = list(define_pattern().keys()))

RMSK_PREFX = ANNOTATION_PATH + '/rmsk'
RMSK_BED = RMSK_PREFX + '.bed.gz'
RMSK_FA = RMSK_PREFX + '.fa'
BOWTIE2_RMSK_INDEX = expand(RMSK_PREFX + '.{NUMBER}.bt2', NUMBER = range(1,5))




GENOME_LINK  = 'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes'
GTF_LINK = 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/GRCh37_mapping/gencode.v28lift37.annotation.gtf.gz'
tRNA_REF = 'http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi19/hg19-tRNAs.tar.gz'
piRNA = 'http://www.regulatoryrna.org/database/piRNA/download/archive/v1.0/bed/piR_hg19_v1.0.bed.gz'
MIR_LINK = 'ftp://mirbase.org/pub/mirbase/CURRENT/hairpin_high_conf.fa.gz'
UNI_VEC_LINK = 'ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec_Core'
RMSK_LINK = 'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz'

def test_filter(test_bool):
    FILTER = ''
    if test_bool:
        FILTER = "| egrep --color=no 'chr[YM]'"
    return FILTER


rule all:
    input:
        HISAT_GENOME_INDEX,
        BOWTIE2_GENOME_INDEX,
        BOWTIE2_UNIVEC_INDEX,
        BOWTIE2_SMRNA_INDEX,
        BOWTIE2_mt_rRNA_INDEX,
        SMALL_RNA_BED,
        SPLIT_GENE_BED,
        piRNA_BED,
        EXONS,
        SPLICE_SITE,
        BOWTIE2_RMSK_INDEX,
        RMSK_BED

#### make RMSK ###
rule download_rmsk:
    input:

    params:
        LINK = RMSK_LINK,
        FILTER = test_filter(config['test'])

    output:
        RMSK_BED
    
    shell:
        'curl {params.LINK} '\
        '| zcat '\
        "| awk '{{print $6,$7,$8,$11,$NF, $10, $12, $13}}' OFS='\\t' "\
        "{params.FILTER}" \
        '| sort -k1,1 -k2,2n -k3,3n '\
        '| bgzip '\
        '> {output}'

rule make_rmsk_fa:
    input:
        FA = GENOME_FA,
        BED = RMSK_BED

    params:

    output:
        RMSK_FA
    
    shell:
        'bedtools getfasta -fi {input.FA} -bed {input.BED} -s -name '\
        '| seqkit rmdup -n '\
        '| seqkit rmdup -s '\
        '> {output}'

rule make_rmsk_index:
    input:
        FA = RMSK_FA
    
    params:
        PREFIX = RMSK_PREFX

    output:
        BOWTIE2_RMSK_INDEX

    shell:
        'bowtie2-build {input.FA} {params.PREFIX}'

#### make genes ####
rule download_GTF:
    input:

    params:
        LINK = GTF_LINK,
        FILTER = test_filter(config['test'])

    output:
        GENE_GTF
    
    shell:
        'curl {params.LINK} |zcat {params.FILTER}  > {output}'


rule make_bed12:
    input:
        GTF = GENE_GTF
    
    output:
        BED12

    shell:
        'cat {input.GTF} '\
        "| grep 'protein_coding' --color=no "\
        "| gtfToGenePred /dev/stdin /dev/stdout "\
        "| genePredToBed > {output}"


rule exon_bed:
    input:
        GTF = GENE_GTF

    output:
        EXONS = EXONS
    
    shell:
        'cat {input.GTF} '\
        "| grep 'protein_coding' --color=no "\
        "| awk '$3==\"exon\"' "\
        "| gtf2bed "\
        "| sort -k1,1 -k2,2n -k3,3n -k6,6 -u "\
        "> {output.EXONS}"



rule make_splice_site:
    input:
        GTF = GENE_GTF
    
    output:
        SPLICE_SITE
    
    shell:
        'hisat2_extract_splice_sites.py {input.GTF} > {output}'


rule make_gene_bed:
    input:
        GTF = GENE_GTF
    
    output:
        BED = GENE_BED
    
    shell:
        'python gtf_to_bed.py {input.GTF} > {output.BED}'


rule make_split_bed:
    input:
        BED = GENE_BED

    params:
        PATH = ANNOTATION_PATH
    
    output:
        SPLIT_GENE_BED

    shell:
        ' python split_bed_for_count.py {params.PATH} '



#### make rRNA ####
rule make_rRNA_mt_fa:
    input:
        FA = GENOME_FA,
        BED = GENE_BED
    
    output:
        FA = mt_rRNA_FA,
        BED = mt_rRNA_BED
    
    shell:
        'python get_rRNA.py {input.FA} {input.BED} {output.BED} {output.FA} '\
        "; cat {output.BED} | awk '$7!~/Mt/' >> {input.BED} "


rule make_rRNA_mt_index:
    input:
        FA = mt_rRNA_FA
    
    params:
        PREFIX = mt_rRNA_PREFIX

    output:
        BOWTIE2_mt_rRNA_INDEX
    
    shell:
        'bowtie2-build {input.FA} {params.PREFIX}'


##### small RNA ########
rule download_piRNA:
    input:
        GENE_BED

    params:
        LINK = piRNA


    output:
        piRNA_BED

    shell:
        'curl {params.LINK} '\
        '| zcat '\
        '| sort -k1,1 -k2,2n -k3,3n '\
        "| awk '{{print $0, \"piRNA\", \"piRNA\"}}' OFS='\\t'" \
        '| bgzip '\
        '> {output} ' \
        '; ccat {output}| zcat >> {input}'



rule download_gtRNA:
    input:
    
    params:
        LINK = tRNA_REF,
        DIR = tRNA_FOLDER

    output:
        BED = GTRNA_BED, 
        FA = GTRNA_FA, 
        ZIP = tRNA_ZIP

    shell:
        'curl {params.LINK} > {output.ZIP} '\
        '; mkdir -p {params.DIR} '\
        '; tar zxvf {output.ZIP} --directory {params.DIR}'


rule make_tRNA_bed:
    input:
        BED = GTRNA_BED,
        GENE_BED = GENE_BED

    output:
        BED = tRNA_BED
    
    shell:
        'cp {input.BED} {output.BED} '\
        ';cat {output.BED} | cut -f1-8 >> {input.GENE_BED}'

rule make_nucleo_tRNA:
    input:
        FA = GTRNA_FA
    
    output:
        FA = NUCLEO_tRNA_FA
    
    shell:
        'seqkit  rmdup -s  {input.FA} ' \
        '| python process_mature_tRNA.py ' \
        '> {output.FA}'


rule make_mt_tRNA:
    input:
        BED = GENE_BED,
        FA = GENOME_FA
    
    output:
        FA = MT_tRNA_FA
    
    shell:
        'cat {input.BED} '\
        "| grep 'Mt_tRNA' "\
        "| bedtools getfasta  -fi {input.FA} -bed - -s -name -tab "\
        "| tr ':' '\\t' "\
        "| awk '{{printf \">%s\\n%s\\n\",$1,$NF}}' "\
        "| sed 's/(-)//g' | sed 's/(+)//g' "\
        "> {output.FA}"

rule make_tRNA:
    input:
        MT_tRNA_FA,
        NUCLEO_tRNA_FA
    
    output:
        tRNA_FA
    
    shell:
        'cat {input} > {output}'

rule get_yRNA:
    input:
        BED = GENE_BED,
        FA = GENOME_FA
    
    output:
        FA = yRNA_FA,
        BED = yRNA_BED
    
    shell:
        'cat {input.BED} '
        "| grep --color=no 'RNY' "\
        "| awk '$4!~/[pP]/' "\
        "| python get_fa.py {input.FA} {output.BED} {output.FA}"

rule get_SRP:
    input:
        BED = GENE_BED,
        FA = GENOME_FA
    
    output:
        FA = SRP_FA,
        BED = SRP_BED
    
    shell:
        'cat {input.BED} '\
        "| awk '$4~/.*7SK$|7SL[0-9]+$/' "\
        "| awk '{{print $0, $3-$2}}' OFS='\\t' "\
        "| awk  '$NF~/299|330|296/'" \
        "| python get_fa.py {input.FA} {output.BED} {output.FA}"


rule get_vault_RNA:
    input:
        BED = GENE_BED,
        FA = GENOME_FA

    output:
        FA = vaultRNA_FA,
        BED = vaultRNA_BED

    shell:
        'cat {input.BED} '\
        "| grep --color=no 'VTRNA' "\
        "| awk '$4!~/[pP]/' "\
        "| python get_fa.py {input.FA} {output.BED} {output.FA}"




rule get_miRNA:
    input:

    params:
        LINK = MIR_LINK

    output:
        miRNA_FA


    shell:
        "curl {params.LINK} "\
        "| seqkit grep -n -r -p 'Homo sapien' " \
        "| seqkit replace -s -p 'U' -r 'T' " \
        '> {output}'

rule make_smallRNA_bed:
    input:
        SMALL_RNA_FA
    
    output:
        SMALL_RNA_BED
    
    shell:
        'python smallRNA_bed.py {input} > {output}'


rule make_smallRNA_fa:
    input:
        tRNA_FA,
        yRNA_FA,
        SRP_FA,
        miRNA_FA,
        vaultRNA_FA
    
    output:
        SMALL_RNA_FA
    
    shell:
        'cat {input} > {output}'


rule make_smallRNA_index:
    input:
        FA = SMALL_RNA_FA
    
    params:
        PREFIX = SMALL_RNA_PREFIX

    output:
        BOWTIE2_SMRNA_INDEX
    
    shell:
        'bowtie2-build {input.FA} {params.PREFIX}'


#### Univec ###################
rule download_univec:
    input:
    
    params:
        LINK = UNI_VEC_LINK
    
    output:
        FA = UniVec_FA
    
    shell:
        'curl {params.LINK} > {output.FA}'


rule make_univec_index:
    input:
        FA = UniVec_FA

    params:
        PREFIX = UniVec_PREFIX

    output:
        BOWTIE2_UNIVEC_INDEX,
    
    shell:
        'bowtie2-build {input.FA} {params.PREFIX}'
        

###### GENOME ######
rule make_hisat_genome_index:
    input:
        FA = GENOME_FA
    
    params:
        PREFIX = GENOME_PREFIX
    output:
        HISAT_GENOME_INDEX
    
    shell:
        'hisat2-build {input.FA} {params.PREFIX}'

rule make_bowtie_genome_index:
    input:
        FA = GENOME_FA
    
    params:
        PREFIX = GENOME_PREFIX
    
    output:
        BOWTIE2_GENOME_INDEX
    
    shell:
        'bowtie2-build {input.FA} {params.PREFIX}'


rule download_genome_fa:
    input:  

    params:
        LINK=GENOME_LINK.format(CHROM = '${CHROM}'),
        FILTER = test_filter(config['test'])

    output:
        GENOME_FA

    shell:
        "for CHROM in $(curl {params.LINK}/md5sum | awk '{{print $2}}'" \
            "| egrep --color=no 'fa.gz$' {params.FILTER}); "\
        'do '\
            'curl {params.LINK}/$CHROM ;'\
        'done |zcat '\
        '> {output}'
