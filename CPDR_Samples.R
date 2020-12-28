##########################
#
#
#Read example in CPDR
#
#
#
#
##########################

setwd("~/dn/Projects/TCGA_PRAD/")
library(psichomics)
library(pheatmap)
library(tidyverse)
library(readxl)
select=dplyr::select
filter=dplyr::filter
##prepareGeneQuant(Sys.glob("~/Documents/Aldrin/CPDR_RNA-seq/STAR/hg19/*Reads*"))
#prepareJunctionQuant(Sys.glob("~/Documents/Aldrin/CPDR_RNA-seq/STAR/hg19/*SJ.out.tab"))

data=loadLocalFiles("/home/dnousome/dn/Projects/TCGA_PRAD/Ethnicity_TCGA_PRAD/")

geneExpr      <- data[[1]]$`Gene expression`
junctionQuant <- data[[1]]$`Junction quantification`


hg19 <- listSplicingAnnotations(assembly="hg19")[[1]]
annotation <- loadAnnotation(hg19)
getSplicingEventTypes()

minReads <- 10 # default

##Quantify only known events
psi <- quantifySplicing(annotation, junctionQuant, 
                        minReads=minReads,
                        eventType =  c("SE", "MXE", "ALE", "AFE", "A3SS", "A5SS"))

##Since only few samples 4 of 17 
psi_miss=apply(psi,1,function(x)sum(is.na(x)))
psi_dt=psi[psi_miss<4,]

names(psi_dt)=sapply(strsplit(names(psi_dt),"_"),'[',1)
##20% of 17

##Ethnicity CPDR
covar=read_xlsx("~/Documents/G/Info_Mgmt/Bioinformatics-CPDR/Datasets/01_7AA_7CA_CaP_RNASeq/02 201210 ExpressionAnalytics RNASeq/EA RNA-Seq QC Data 101812/NextGen Seq RNA Specimens Oct12-2012.xlsx",sheet=3)



ethnicity=covar %>% 
  dplyr::filter(`Sample Type`=="TumorRNA") %>%
  arrange(match(SampleID,names(psi_dt)))

psi_dt1=psi_dt %>% dplyr::select(ethnicity$SampleID)
#psi_dt_lm=t(psi_dt)


psi_p=apply(psi_dt1,1,function(x){
  summary(lm(unlist(x)~ethnicity$Race))$coefficients[2,4]
})

ethnicity_dt=ethnicity %>% select(SampleID,Race) %>% column_to_rownames('SampleID') 
psi_q=p.adjust(psi_p,"fdr")



###Load in the TCGA results and see if these look different as well
tcga_splice=read_csv("TCGA_splicing.csv")


cpdr_splice=data.frame(cpdr_p=psi_p,cpdr_q=psi_q) %>% rownames_to_column("Splice")
splice_res=inner_join(tcga_splice,cpdr_splice,by="Splice")

top_tcgatovalidate=splice_res %>% dplyr::filter(TCGA_q<0.05) %>% filter(cpdr_p<0.1)


####################
sum(psi_p<0.05,na.rm=T)
sum(psi_q<0.05,na.rm=T)




psi_dt2=psi_dt1[psi_p<0.001 & !is.na(psi_p),]
psi_dt2=psi_dt1[order(psi_p)[1:25],]


rownames(psi_dt2)=make.unique(sapply(strsplit(rownames(psi_dt2),"_"),function(x)paste0(x[[1]],"_",x[length(x)])))

pheatmap(psi_dt3,annotation_col = ethnicity_dt,
         treeheight_row = 0, treeheight_col = 0,
        height = 8,width=5,filename = "CPDR.png",
         scale="row",
         fontsize = 8,clustering_method = "ward.D",show_colnames = F)



###########Using top to validate
psi_dt2=psi_dt1[top_tcgatovalidate$Splice,]


rownames(psi_dt2)=make.unique(sapply(strsplit(rownames(psi_dt2),"_"),function(x)paste0(x[[1]],"_",x[length(x)])))

pheatmap(psi_dt2,annotation_col = ethnicity_dt,
         treeheight_row = 0, treeheight_col = 0,
         height = 8,width=5,
         filename = "CPDR_validate.png",
         scale="row",
         fontsize = 8,clustering_method = "ward.D",show_colnames = F)




#################EXAMINE THE METHYL GENES
##load in the Kallisto data and then run it run it and run it
library(tximport)
library(EnsDb.Hsapiens.v86)
txdf <- transcripts(EnsDb.Hsapiens.v86, return.type="DataFrame")
#gdf=read_tsv("~/dn/Projects/LTF-Shashwat/geneid")
tx2gene <-  as.data.frame(txdf[,c("tx_id","gene_id")])

setwd("~/Documents/Aldrin/CPDR_RNA-seq//kallisto_quant/")
files=Sys.glob("*/*.h5")
files1=data.frame(file=sapply(strsplit(files,"_"),'[',1))

txi <- tximport(files[1:14], type = "kallisto", tx2gene = tx2gene,ignoreTxVersion = T)
covar=read_xlsx("~/Documents/G/Info_Mgmt/Bioinformatics-CPDR/Datasets/01_7AA_7CA_CaP_RNASeq/02 201210 ExpressionAnalytics RNASeq/EA RNA-Seq QC Data 101812/NextGen Seq RNA Specimens Oct12-2012.xlsx",sheet=3)

covar_sub=covar %>% dplyr::filter(!is.na(Race))
rownames(covar_sub)=covar_sub$SampleID
dds <- DESeqDataSetFromTximport(txi, covar_sub,~Race)
dds_res <- DESeq(dds, test = 'LRT', reduced = ~ 1, fitType = 'local')

library('org.Hs.eg.db')
gencode <- gsub('\\..*', '', names(recount_genes))

## Find the gene information we are interested in
gene_info <- select(org.Hs.eg.db, gencode, c('ENTREZID', 'GENENAME', 'SYMBOL',
                                             'ENSEMBL'), 'ENSEMBL')


methyl_geneinfo=gene_info %>% dplyr::filter(SYMBOL %in% methylgenes)
top_methyl=data.frame(results(dds_res)) %>% 
  rownames_to_column("GENENAME") %>% 
  separate(GENENAME,sep="[.]",into="GENE") %>%
  inner_join(.,methyl_geneinfo,by=c("GENE"="ENSEMBL")) %>%
  arrange(pvalue) %>%
  mutate(pfdr=p.adjust(pvalue,"fdr"))

write_csv(top_methyl,"TopRNAMethylation_cpdr.csv")  

