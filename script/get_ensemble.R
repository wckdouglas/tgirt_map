#!/use/bin/env Rscript
suppressMessages(library(biomaRt))
suppressMessages(library(dplyr))
suppressMessages(library(readr))
suppressMessages(library(purrr))
suppressMessages(library(stringr))

ensembl<-useEnsembl(biomart="ensembl",GRCh=37,dataset = "hsapiens_gene_ensembl")
bed_fields <- c('chromosome_name',
                'start_position',
                'end_position',
                'external_gene_name',
                'strand',
                'gene_biotype',
                'ensembl_gene_id')
getBM(attributes = bed_fields, mart = ensembl) %>%
    set_names(c('chrom','start','end',
                'name','strand','biotype',
                'id')) %>%
    mutate(score=0) %>%
    dplyr::select(chrom,start,end,name,score,strand, biotype, id) %>%
    mutate(strand = ifelse(strand < 1, '-','+')) %>%
    filter(biotype!='TEC') %>%
	mutate(chrom = case_when(
            .$chrom=='MT'~ "chrM",
	    TRUE~str_c('chr',.$chrom))) %>%
    arrange(chrom,start,end) %>%
    write.table("", col.names=F,sep='\t', quote=F, row.names=F)
