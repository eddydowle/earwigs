#DESeq2 prelim analysis of Spotty 2018 RNA-seq gonad datasets against TrinitySI18_1 transcriptome.
#All samples: AI14A, SI16, SI18, SIP15 experiments
#95 samples as SI18_19G has 0 counts in table for unknown reason
#Stages selected based on my take on histo notes and Prelim_clustering analyses
#Exploratory plotting of some target genes - based on Oscar's scripts for generating early boxplots for blueahead

library(DESeq2)
library(org.Dm.eg.db)
library(ggpubr)
library(biomaRt)
library(stringr)
library(tidyverse)
setwd("/Users/edwinadowle/Documents/Earwigs/NugenRunsNov2018/hiseqruns/salmon_quants_june_gsa/genes_salmon/")

################################################
#sort out Trinotate annotations matrix to be useful as column records contain useless repetitive shit
#read in the table, ignore rando characters, replace '.' with 'NA'
annotations <- read.table('../trinotate_onessamplepergroup_annotation_report.xls', sep="\t", header=TRUE, row.names=NULL, quote="", na.strings=".")

colnames(annotations)[which(names(annotations) == "...trinotate.drosophila.dmel.all.translation.r6.27_BLASTX")] <- "Dmel_Top_BLASTX"
colnames(annotations)[which(names(annotations) == "...trinotate.elegans.Caenorhabditis_elegans.WBcel235.pep_BLASTX")] <- "Celegan_Top_BLASTX"
colnames(annotations)[which(names(annotations) == "...trinotate.drosophila.dmel.all.translation.r6.27_BLASTP")] <- "Dmel_Top_BLASTP"
colnames(annotations)[which(names(annotations) == "...trinotate.elegans.Caenorhabditis_elegans.WBcel235.pep_BLASTP")] <- "Celegan_Top_BLASTP"
colnames(annotations)
#successively isolate key info in separate columns (will be added as new columns at the end)
#use this for testing scripts https://regex101.com/

annotationsSimplified = annotations %>% mutate(sprot_BLASTX_gene=gsub('_.*',"",sprot_Top_BLASTX_hit))
annotationsSimplified$sprot_BLASTX_gene = gsub('"','',annotationsSimplified$sprot_BLASTX_gene)

annotationsSimplified = annotationsSimplified %>% mutate(sprot_BLASTX_name = str_extract(sprot_Top_BLASTX_hit, "Full=.*?;"))
annotationsSimplified$sprot_BLASTX_name = gsub('Full=','',annotationsSimplified$sprot_BLASTX_name)
annotationsSimplified$sprot_BLASTX_name = gsub(';','',annotationsSimplified$sprot_BLASTX_name)

annotationsSimplified = annotationsSimplified %>% mutate(sprot_BLASTP_gene=gsub('_.*',"",sprot_Top_BLASTP_hit))
annotationsSimplified$sprot_BLASTP_gene = gsub('"','',annotationsSimplified$sprot_BLASTP_gene)

annotationsSimplified = annotationsSimplified %>% mutate(sprot_BLASTP_name = str_extract(sprot_Top_BLASTP_hit, "Full=.*?;"))
annotationsSimplified$sprot_BLASTP_name = gsub('Full=','',annotationsSimplified$sprot_BLASTP_name)
annotationsSimplified$sprot_BLASTP_name = gsub(';','',annotationsSimplified$sprot_BLASTP_name)

			
annotationsSimplified = annotationsSimplified %>% mutate(Dmel_BLASTX=gsub('\\^.*',"",Dmel_Top_BLASTX))
annotationsSimplified$Dmel_BLASTX = gsub('"','',annotationsSimplified$Dmel_BLASTX)

annotationsSimplified = annotationsSimplified %>% mutate(Celegan_BLASTX=gsub('\\^.*',"",Celegan_Top_BLASTX))
annotationsSimplified$Celegan_BLASTX = gsub('"','',annotationsSimplified$Celegan_BLASTX)

annotationsSimplified = annotationsSimplified %>% mutate(Dmel_BLASTP=gsub('\\^.*',"",Dmel_Top_BLASTP))
annotationsSimplified$Dmel_BLASTP = gsub('"','',annotationsSimplified$Dmel_BLASTP)

annotationsSimplified = annotationsSimplified %>% mutate(Celegan_BLASTP=gsub('\\^.*',"",Celegan_Top_BLASTP))
annotationsSimplified$Celegan_BLASTP = gsub('"','',annotationsSimplified$Celegan_BLASTP)



#pull out and make a list of all unique Danio protein IDs, gsub the version (remove everything after '.'), remove 'NAs'
protList <- annotationsSimplified$Dmel_BLASTX
protList = gsub('\\..*','',protList)
protList <- unique(protList[!is.na(protList)])
protList

protList <- annotationsSimplified$Celegan_BLASTX
#protList = gsub('\\..*','',protList)
protList <- unique(protList[!is.na(protList)])
protList

################################################
#Use biomaRt to use protein IDs to add extra annotation into
#https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.html#how-to-build-a-biomart-query
listMarts()

#Select database and dataset for each fish species
ENSEMBL_Dm=useMart("ensembl", dataset="drerio_gene_ensembl")
ENSEMBL_Ce=useMart("ensembl", dataset="oniloticus_gene_ensembl")

####To query the databases
#attributes = those to retrieve, filters = input to query (i.e. gene names), values = values for filters (i.e. in gene_list) 
listFilters(ENSEMBL_Dr)
listFilters(ENSEMBL_On)
listFilters(ENSEMBL_Lb)
#useful: ensembl_gene_id  ensembl_peptide_id  external_gene_name  

listAttributes(ENSEMBL_Dr)
#useful: ensembl_gene_id  ensembl_peptide_id  external_gene_name  description

####RETRIEVE REQUIRED INFO:
####Use the Danio protein ID list 'protList' from Trinotate to generate table with extra data
annotations_with_Dr <- getBM(attributes=c('external_gene_name', 'ensembl_gene_id', 'ensembl_peptide_id', 'description'),
                             filters='ensembl_peptide_id',
                             values = protList,
                             mart = ENSEMBL_Dr)

#merge the extra Danio info with the Trinotate table
#first, gsub out protein version 
annotationsSimplified$Danio_blastx_ID = gsub('\\..*','',annotationsSimplified$Danio_blastx_ID)
annotationsSimplified_Dr <- merge(annotationsSimplified, annotations_with_Dr, by.x="Danio_blastx_ID", by.y="ensembl_peptide_id", all.x=TRUE)
head(annotationsSimplified_Dr)
#re-order columns and export
#annotationsSimplified_Dr <- annotationsSimplified_Dr[c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,1,34,33,35)]
write.table(annotationsSimplified_Dr, "../trinotate_onessamplepergroup_annotation_report.fullextras.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE) 

#To keep only some columns and export as new table
annotationsUseful_Dr <- select(annotationsSimplified_Dr, gene_id, transcript_id, sprot_BLASTX_gene, sprot_BLASTX_name, sprot_BLASTP_gene, sprot_BLASTP_name, Oreo_blastx_ID, Labrus_blastx_ID, Danio_blastp_ID, Oreo_blastp_ID, Labrus_blastp_ID, Danio_blastx_ID, ensembl_gene_id, external_gene_name, description)
write.table(annotationsUseful_Dr, "TrinitySI2018.1.cdh99_annotations_wDanioAnnots_Useful2.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE) 

annotationsUseful_Dr <-read.table('../Trinity_onepergroup_useful_annotations_transcripts.txt',header=T,row.names=NULL,sep='\t',quote='')
colnames(annotationsUseful_Dr)
colnames()
#Create new table with single entry per gene (so can merge with DESeq dds later) (dyplr - takes a while)

#maybe just to the ones for which there is multiple isoforms to reduce stress on R
annotationsUseful_Dr <-read.table('../Trinity_onepergroup_useful_annotations_transcripts.txt',header=T,row.names=NULL,sep='\t',quote='')
annotationsUseful_Dr_multiisoform<-annotationsUseful_Dr %>% 
  group_by(gene_id) %>% 
  filter(n()>1)


annotationsUseful_Dr_multiisoformsingle <- annotationsUseful_Dr_multiisoform %>% group_by(gene_id) %>%
  summarize(transcript_id_grp = paste(unique(transcript_id[!is.na(transcript_id)]), collapse = ","),
            sprot_BLASTX_name_grp = paste(unique(sprot_BLASTX_name[!is.na(sprot_BLASTX_name)]), collapse = ","),
            sprot_BLASTP_gene_grp = paste(unique(sprot_BLASTP_gene[!is.na(sprot_BLASTP_gene)]), collapse = ","),
            sprot_BLASTP_name_grp = paste(unique(sprot_BLASTP_name[!is.na(sprot_BLASTP_name)]), collapse = ","),
            Dmel_BLASTP_grp = paste(unique(Dmel_BLASTP[!is.na(Dmel_BLASTP)]), collapse = ","),
            Dmel_BLASTX_grp = paste(unique(Dmel_BLASTX[!is.na(Dmel_BLASTX)]), collapse = ","),
            Dmel_external_gene_name_grp = paste(unique(Dmel_external_gene_name[!is.na(Dmel_external_gene_name)]), collapse = ","),
            Dmel_ensembl_gene_id_grp = paste(unique(Dmel_ensembl_gene_id[!is.na(Dmel_ensembl_gene_id)]), collapse = ","),
            Dmel_external_gene_name_blastp_grp = paste(unique(Dmel_external_gene_name_blastp[!is.na(Dmel_external_gene_name_blastp)]), collapse = ","),
            Dmel_ensembl_gene_id_blastp_grp = paste(unique(Dmel_ensembl_gene_id_blastp[!is.na(Dmel_ensembl_gene_id_blastp)]), collapse = ","),
            Celegan_BLASTP_grp = paste(unique(Celegan_BLASTP[!is.na(Celegan_BLASTP)]), collapse = ","),
            Celegan_BLASTX_grp = paste(unique(Celegan_BLASTX[!is.na(Celegan_BLASTX)]), collapse = ","),
            Cele_external_gene_name_grp = paste(unique(Cele_external_gene_name[!is.na(Cele_external_gene_name)]), collapse = ","),
            Cele_ensembl_gene_id_grp = paste(unique(Cele_ensembl_gene_id[!is.na(Cele_ensembl_gene_id)]), collapse = ","),
            Cele_external_gene_name_blastp_grp = paste(unique(Cele_external_gene_name_blastp[!is.na(Cele_external_gene_name_blastp)]), collapse = ","),
            Cele_ensembl_gene_id_blastp_grp = paste(unique(Cele_ensembl_gene_id_blastp[!is.na(Cele_ensembl_gene_id_blastp)]), collapse = ","))

write.table(annotationsUseful_Dr_single, "../Trinity_onepergroup_useful_annotations_gene.txt", sep="\t", row.names=FALSE, header=TRUE, quote=FALSE) 



#Create new table with single entry per gene (so can merge with DESeq dds later) (aggregate - takes aaages, never finished - could try overnight)
#annotationsUseful_Dr_single <- aggregate(annotationsUseful_Dr[,2:15], list(annotationsUseful_Dr[,1]), function(x) paste0(unique(x[!is.na(x)])))


####For a list of candidate gene names in 'gene_list',  create table for those genes only
#Import a candidate (zebrafish) gene list for plotting (gene names must be correct as in database)
gene_list = data.frame(read.csv("gene_list.txt",header = FALSE,stringsAsFactors = FALSE))

#Create table 
#gene_list_Dr_annots <- getBM(attributes=c('external_gene_name', 'ensembl_gene_id', 'ensembl_peptide_id', 'description'),
      #filters='external_gene_name',
      #values = gene_list,
      #mart = ENSEMBL_Dr)


################################################
#OSCAR'S ORIGINAL WAY - adds annotations to expression results after analysis - less useful
##Import candidate gene list for plotting (gene names must be correct as in database)
gene_list = data.frame(read.csv("gene_list.txt",header = FALSE,stringsAsFactors = FALSE))
colnames(gene_list) = "SYMBOLS"

#Add gene IDs
#keytype=return results based on this column of db
ENSEMBL = mapIds (org.Dr.eg.db, gene_list$SYMBOLS, keytype="SYMBOL", column="ENSEMBL",multiVals = CharacterList)

table_SYMBOL_ENSEMBL <- data.frame()

for (i in 1:length(names(ENSEMBL))){
  temp <- unlist(ENSEMBL[i])
  SYMBOL_temp <- names(temp)
  ENSEMBL_temp <- unlist(temp[[1]])
  df_temp <- data.frame("SYMBOL" = rep(SYMBOL_temp,length(ENSEMBL_temp)),"ENSEMBL" = ENSEMBL_temp)
  table_SYMBOL_ENSEMBL <- rbind(table_SYMBOL_ENSEMBL, df_temp)
}

table_SYMBOL_ENSEMBL

################################################
###IMPORT RNA-seq DATA AND SET UP DESIGN
##Create and view data frame from counts matrix, indicating the file contains column headings and row headings. Filter genes with <10 total read count to reduce file size. View data.<p>

#import data and view
counts <- read.table('Spotty2018_matrix_29_04_2019.gene.counts.matrix.rounded.txt', header=TRUE, row.names=1)
head(counts)
dim(counts) #141205 95

#set-up counts matrix, remove rows with <10 counts total, to reduce size
counts = as.matrix(counts)
counts<-counts[apply(counts,1,sum)>=10,]
head(counts)
dim(counts) #140737 x 95
summary(rowSums(counts))
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#10       54      109     6676      356 42555888

#read in sample info table
coldata <- read.csv('coldata.csv', row.names=1)
head(coldata)
coldata <- coldata[,c("batch","condition")]
coldata

#Re-order coldata by stages
condOrder = c("CF", "NBF", "ET", "MT", "LT", "TP", "IP")
batchOrder = c("AI14A", "SIP15", "SI16", "SI18")
coldata$condition <- factor(coldata$condition, levels=condOrder)
coldata <- coldata[order(coldata$condition),]
coldata

#Re-order counts to match new coldata order
counts <- counts[, rownames(coldata)]
head(counts)

#check naming and sample order between counts and coldata is consistent
all(rownames(coldata) %in% colnames(counts))
#TRUE
all(rownames(coldata) == colnames(counts))
#TRUE

##MERGE annotations with couns into one table called counts_Dr
counts_Dr <- merge(counts, annotationsUseful_Dr_single, by.x="row.names", by.y="gene_id", all.x=TRUE)
head(counts_Dr)
names(counts_Dr)
write.table(counts_Dr, "TrinitySI2018_counts_wDr.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE) 


##create DESeq2 data set (dds) for all samples
#Specify design that includes 'batch' and 'condition' ('condition' last as is most important)
dds = DESeqDataSetFromMatrix(countData = counts, 
                             colData = coldata,
                             design = ~ batch + condition)
dds

##Add Danio feature data
featureData <- data.frame(counts_Dr[,c(1,97:110)],row.names=1)
mcols(dds) <- DataFrame(mcols(dds), featureData)
rownames(mcols(dds))

#normalize data
dds <- estimateSizeFactors(dds)
matrix_normalized <- as.data.frame(counts(dds, normalized=TRUE))
matrix_normalized_log <- log2(matrix_normalized + 1)


###############################
#PLOT GENES
library("ggplot2")
genes_selected <- read.csv("cyp19a1a_amh_contigs.csv",stringsAsFactors = FALSE)
colnames(genes_selected) <- c("symbol")

selected_contigs <- counts_Dr %>% filter(external_gene_name_grp %in% genes_selected$symbol) %>% select(Row.names,external_gene_name_grp) 

#example
d<-plotCounts(dds, gene="TRINITY_DN9355_c0_g1", intgroup="condition", returnData = TRUE)

ggplot(d, aes(x=condition, y=log(count), group=1)) + 
  geom_point(color="black",pch=19, size=1.5) + # Change pint size, size
  stat_summary(fun.y = "mean", geom="line", color="blue") + # Change color line, color
  stat_summary(fun.data = mean_sd, geom="errorbar", size=.5, width=.3) + #
  stat_summary(fun.y = "mean", geom="point", color="blue", shape=21, size=3, fill="white") +
  xlab("Sex change stage") + ylab("Normalized expression") +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  

#To plot in a grid of multiple genes
genes_selected <- read.csv("cyp19a1a_amh_contigs.csv",stringsAsFactors = FALSE)
colnames(genes_selected) <- c("symbol")

selected_contigs <- counts_Dr %>% filter(external_gene_name_grp %in% genes_selected$symbol) %>% select(Row.names,external_gene_name_grp) 

table_plot<-c()

for (i in 1:nrow(selected_contigs)){
  contig<-selected_contigs$Row.names[i] 
  d <- plotCounts(dds, gene=contig, intgroup="condition", returnData = TRUE)
  d$gene<-selected_contigs$external_gene_name_grp[i]
  d$contig<-contig
  if(length(table_plot)==0){
    table_plot<-d
  }
  else{table_plot<-rbind(table_plot,d)}
}
selected_contigs$Row.names[1]

library("ggplot2")
ggplot(table_plot, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  facet_wrap(~gene*contig)

ggplot(table_plot, aes(x= stage, y = value, group=1)) +
  geom_point(color="black",pch=19, size=1.5) + # Change pint size, size
  stat_summary(fun.y = "mean", geom="line", color="red") + # Change color line, color
  stat_summary(fun.data = mean_se, geom="errorbar", size=.5, width=.3) + #
  stat_summary(fun.y = "mean", geom="point", color="red", shape=21, size=3, fill="white") +
  xlab("Sex change stage") + ylab("Normalized expression") + 
  theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x = element_blank()) +  # Remove stages
  facet_grid(symbol ~ .) +
  theme(strip.text.y = element_text(size = 16))

library("ggplot2")
ggplot(table_plot, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  facet_wrap(~gene)

mcols(dds)
names(dds)
mcols(mcols(dds))
$external_gene_name_grp=="amh")
