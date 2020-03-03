#######################################
#
#
#Identify if there are splicing differences between AA/CA Mne
#
#
#
######################################
library(tidyverse)
library(psichomics)
library(readxl)
library(pheatmap)
# Set download folder
folder <- "~/dn/Projects/TCGA_PRAD/"
select=dplyr::select
filter=dplyr::filter
# Download and load most recent junction quantification and clinical data from
# TCGA/Firebrowse for Pca
data <- loadFirebrowseData(folder=folder,
                           cohort="PRAD",
                           data=c("clinical", "junction_quantification",
                                  "RSEM_genes"),
                           date="2016-01-28")


geneExpr=data[[1]]$`Gene expression`
filter <- filterGeneExpr(geneExpr)
geneExprFiltered <- geneExpr[filter, ]
junctionQuant <- data[[1]]$`Junction quantification (Illumina HiSeq)`
sampleInfo    <- data[[1]]$`Sample metadata`
clinical      <- data[[1]]$`Clinical data`

geneExprNorm <- normaliseGeneExpression(geneExprFiltered)


# Available alternative splicing annotation
annotList <- listSplicingAnnotations()
annotList

hg19 <- listSplicingAnnotations(assembly="hg19")[[1]]
annotation <- loadAnnotation(hg19)
getSplicingEventTypes()

minReads <- 10 # default

psi <- quantifySplicing(annotation, junctionQuant, minReads=minReads)


covar=data.frame(ID=toupper(clinical$patient.bcr_patient_barcode),
                 Gleason=clinical$patient.stage_event.gleason_grading.gleason_score,
                 AgeDx=clinical$patient.age_at_initial_pathologic_diagnosis)

ethnicity=read_xlsx("~/dn/Projects/PublicData/TCGA_PRAD_ethnicityEIGEN.xlsx",skip = 3) %>%
  select(ID=`...1`,SIR=`Self-identified race`,EIGEN=`...4`) %>%
  filter(EIGEN %in% c("AA",'EA')) %>%
  left_join(.,covar,by=c("ID")) %>%
  mutate(Gleason_bin=ifelse(Gleason==6,0,1))
  

id_psi=sapply(strsplit(colnames(psi),"-"),function(x)paste(x[1],x[2],x[3],sep="-"))
tumnorm_psi=substr(colnames(psi),14,15)

##01-09Tumor
##10-19 Normal
#table(tumnorm_psi)
norm_prad=tumnorm_psi %in% "01"& id_psi %in% ethnicity$ID
  
  
psi_dt=psi[,norm_prad]  

  


psi_miss=apply(psi_dt,1,function(x)sum(is.na(x)))
psi_dt=psi_dt[psi_miss<50,]


###########
psi_ids1=substr(colnames(psi_dt),1,12)

ethnicity1=ethnicity %>% filter(ID %in% psi_ids1) %>%
  arrange(match(ID,psi_ids1))

#psi_dt_lm=t(psi_dt)


psi_p=apply(psi_dt,1,function(x){
  summary(lm(unlist(x)~ethnicity1$EIGEN+ethnicity1$Gleason_bin+ethnicity1$AgeDx))$coefficients[2,4]
})

#x=psi_dt["A5SS_1_-_46655156_46655130_46655029_POMGNT1",]
#ggplot() +geom_violin(aes(x=ethnicity1$EIGEN,y=unlist(x)))

##########Plot heatmap first
ethnicity2=ethnicity1 %>% mutate(ID1=colnames(psi_top_dt)) %>%
  dplyr::select(ID1,Race=EIGEN) %>% column_to_rownames('ID1')

psi_top_dt=psi_dt[psi_p<0.00000005 & !is.na(psi_p),]
rownames(psi_top_dt)=make.unique(sapply(strsplit(rownames(psi_top_dt),"_"),function(x)paste0(x[[1]],"_",x[length(x)])))
pheatmap(psi_top_dt,annotation_col = ethnicity2,scale = "column",
         clustering_method = "complete",show_colnames = F)





psi_top=psi_p[psi_p<.05]

topforpath=unique(sapply(strsplit(names(psi_top),"_"),function(x)x[length(x)]))
topforpath=topforpath[!is.na(topforpath)]



########USE FDR
psi_fdr=p.adjust(psi_p,method = "fdr")
psi_fdr=psi_fdr[psi_fdr<.1]

topforpath=unique(sapply(strsplit(names(psi_fdr),"_"),function(x)x[length(x)]))
topforpath=topforpath[!is.na(topforpath)]


#################Clsuterprofiler
library(clusterProfiler)
library(EnsDb.Hsapiens.v86)
edb=EnsDb.Hsapiens.v86
edb1=data.frame(genes(edb)) %>% dplyr::select("gene_id","gene_name","entrezid")


##Convert to ENTREZ ID
path_entrez=genes(edb, filter = ~ gene_name %in% topforpath)$entrezid
path_entrez=unique(unlist(path_entrez))


kk <- enrichMKEGG(gene         = path_entrez,
                 organism     = 'hsa',
                 pvalueCutoff = .1)


kkplot=setReadable(kk,'org.Hs.eg.db','ENTREZID')

cnet_kk=cnetplot(kkplot,showCategory = 5)
cnet_kk

barplot(kk, showCategory=20)
emapplot(kk)

library(ReactomePA)
path_ptah=enrichPathway(path_entrez)


#formulas=lapply(paste0(colnames(psi_dt_lm),"ethnicity1$EIGEN"),formula)
#lm_res=lapply(formulas,function(x)glm(x,data=dt,family=binomial))
#lapply(lm_res,tidy)



psi_sd=apply(psi_dt,1,function(x)mad(x,na.rm = T))

psi_dt1=psi_dt[order(-psi_sd),]
psi_dt1=psi_dt1[1:50,]

pheatmap(psi_dt1,show_colnames = F)


# Check the identifier of the splicing events in the resulting table
events <- rownames(psi)
head(events)
# Group by normal and tumour samples
types  <- createGroupByAttribute("Sample types", sampleInfo)
stages <- createGroupByAttribute(
  "patient.stage_event.gleason_grading.gleason_score", clinical)
groups <- list()




for (i in c("6", "9")) {
  stage <- Reduce(union,
                  stages[grep(sprintf("stage %s[a|b|c]{0,1}$", i), names(stages))])
  # Include only tumour samples
  stageTumour <- names(getSubjectFromSample(tumour, stage))
  elem <- list(stageTumour)
  names(elem) <- paste("Tumour Stage", toupper(i))
  groups <- c(groups, elem)
}
groups <- c(groups, Normal=list(normal))