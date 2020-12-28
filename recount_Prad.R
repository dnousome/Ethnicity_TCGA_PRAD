##################################
#
#
#
#
#
##################################
library(readxl)
library(tidyverse)
library(recount)
library(DESeq2)


#load in library(the ethnicity in the PRD
s=read_xlsx("~/dn/Projects/PublicData/TCGA_PRAD_ethnicity.xlsx",skip=3)


##load the recount2 data
load("~/dn/Projects/PublicData/RNA-seq/rse_gene_prostate.Rdata")
rse <- scale_counts(rse_gene)


###
methylgenes=c("TRMT5","BUD23","TRMT61A","TRMT9B","MEPCE","TRMT1","TRMT2A","FBL","MRM1","DIMT1","TRMT11","TRMT10A",
  "METTL2A","MRM2","NOP2","ADAT2","CDA","PUS7","RPUSD2","RPUSD4","PUSL1","DKC1",
  "TRMT12","METTL27","TRMT61B","ALKBH8","BCDIN3D","TRMT1L","TRMT2B","FBLl1","MRM3","TARBP1",
  "TFB1M","TFB2M","THUMPD3","THUMPD2","TRMT10B","TRMT10C",
  "FTSJ1","FTSJ3","CMTR1","CMTR2","NSUN2","NSUN5","NSUN6","NSUN6","NSUN7","NSUN3",
  "ADAT3","ADAT1",'ADARB1',"ADAD1","ADAD2","ADARB2",
  "METTL6","METTL8","METTL2B","METTL2A","TRMT10C","TRMT10B","TRMT10A",
  'VIRMA', 'YTHDC2', 'YTHDF2', 'ALKBH5', 'METTL14', 'METTL3') 


##SUBSET TO THOSE WITH AA/EA NAMES
eigen=s %>% dplyr::select(ID="...1",EIGEN="...8") %>%
  dplyr::filter(EIGEN %in% c("EA","AA")) 

#rse$xml_bcr_patient_barcode
commonids=intersect(rse$xml_bcr_patient_barcode,eigen$ID)

eigen=eigen %>% dplyr::filter(ID %in% commonids) %>% arrange(match(ID,commonids))
rse_subset=rse[,rse$xml_bcr_patient_barcode %in% commonids]

rse_subset=rse_subset[,rse_subset$gdc_cases.samples.sample_type=="Primary Tumor"]
rse_subset=rse_subset[,order(rse_subset$xml_bcr_patient_barcode,rse_subset$reads_downloaded)]
rse_subset=rse_subset[,!duplicated(rse_subset$xml_bcr_patient_barcode )]
rse_subset=rse_subset[,order(match(rse_subset$xml_bcr_patient_barcode,commonids))]
rse_subset$EIGEN=eigen$EIGEN

dds <- DESeqDataSet(rse_subset, ~ EIGEN)
dds_res <- DESeq(dds, test = 'LRT', reduced = ~ 1, fitType = 'local')

library('org.Hs.eg.db')
gencode <- gsub('\\..*', '', names(recount_genes))

## Find the gene information we are interested in
gene_info <- select(org.Hs.eg.db, gencode, c('ENTREZID', 'GENENAME', 'SYMBOL',
                                             'ENSEMBL'), 'ENSEMBL')
methyl_geneinfo=gene_info %>% filter(SYMBOL %in% methylgenes)
top_methyl=data.frame(results(dds_res)) %>% 
  rownames_to_column("GENENAME") %>% 
  separate(GENENAME,sep="[.]",into="GENE") %>%
  inner_join(.,methyl_geneinfo,by=c("GENE"="ENSEMBL")) %>%
  arrange(pvalue) %>%
  mutate(pfdr=p.adjust(pvalue,"fdr"))

write_csv(top_methyl,"TopRNAMethylation_TCGA.csv")  
