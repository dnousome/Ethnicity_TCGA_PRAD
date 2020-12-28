####################################3
#
#
#TIL and Somatic Mutations
#
#
#####################################


library(tidyverse)
library(survival)

setwd("~/Downloads")
s1=read_tsv("TME/xCell_data_mRNA_seq_fpkm_capture_xCell_0818110619.txt")
s2=read_tsv("TME/xCell_data_mRNA_seq_fpkm_polya_xCell_0825110619.txt")


#names(tme)[!names(tme) %in% unique(mutations_dt$Tumor_Sample_Barcode)]

tme=cbind(s2,s1[,-1])
tme_comb=tme[,!duplicated(names(tme))]
tme_dt=data.frame(t(tme_comb[,-1])) %>% rownames_to_column("ID")


mutations1=read_tsv("TME/data_mutations_extended.txt")
mutations2=read_tsv("TME/data_mutations_mskcc.txt")

###Remove if fewer than 5
mutations_dt=mutations1 %>% 
  mutate(Tumor_Sample_Barcode=gsub("-",".",Tumor_Sample_Barcode)) %>%
  mutate(Tumor_Sample_Barcode=gsub("SUTC","SU2C",Tumor_Sample_Barcode)) 


tmb=mutations_dt %>% 
  group_by(Tumor_Sample_Barcode) %>% 
  tally() %>% 
  ungroup() %>% 
  mutate(Tumor_Sample_Barcode=gsub("-",".",Tumor_Sample_Barcode)) %>%
  mutate(Tumor_Sample_Barcode=gsub("SUTC","SU2C",Tumor_Sample_Barcode)) 

####39640 mutations with 5 or more samples taht have it
prad_samps=unique(mutations_dt$Tumor_Sample_Barcode)
prad_samps[!prad_samps %in% tme_dt$ID]

mutations_dt=mutations_dt %>%
  group_by(Hugo_Symbol) %>%
  filter(Tumor_Sample_Barcode %in% tme_dt$ID)%>%
  filter(n() >=10) %>%
  ungroup()



x=split(mutations_dt,mutations_dt$Hugo_Symbol)[[1]]
x=mutations_dt %>% filter(Hugo_Symbol=="IL1R2")

muts_tme=lapply(split(mutations_dt,mutations_dt$Hugo_Symbol),function(x){
    print(x$Hugo_Symbol[1])
    if(nrow(x)>0){
     x1=x %>% mutate(MUT=1) %>% 
     select(Tumor_Sample_Barcode,MUT) %>%
     ###IF ONE MUT
      distinct(.,.keep_all = T) %>%
      full_join(tme_dt,.,by=c("ID"="Tumor_Sample_Barcode")) %>%
      mutate(MUT=ifelse(MUT %in% 1,1,0))
     x2=select(x1,starts_with("X"))
     x2_mut=x1$MUT
      p.genes=apply(x2,2,function(y){
        summary(lm(y~x2_mut))$coefficients[2,4]
      })
    
    
    data.frame(tme$X1,p.genes)
    }
    })


muts_tme_res=lapply(names(muts_tme),function(x){
  muts_tme[[x]] %>% 
    mutate(FDR=p.adjust(p.genes,method="fdr")) %>%
    filter(FDR<0.1) %>% mutate(Gene=x)
}) %>% bind_rows()

muts_tme_res_top=muts_tme_res %>% arrange(p.genes)


x="BRCA1"
y='Keratinocytes'

pgenes<-c('SPOP',"TP53","PTEN","FOXA1","CDKN1B","RB1","CHD1","FOXA1","BRCA2",
          "ATM","PIK3CA","CDK12","AR","AKT1","BRAF","BRCA1","CTNNB1","ERCC2",   
          "THSD7B","MED12","NIPA2","PIK3CA","SCN11A")
muts_tme_res_top %>% filter(Gene %in% pgenes)

dt=filter(mutations_dt,Hugo_Symbol==x)
dt=dt %>% mutate(MUT=1) %>% 
  select(Tumor_Sample_Barcode,MUT) %>%
  distinct(.,.keep_all = T) %>%
  full_join(tme_dt,.,by=c("ID"="Tumor_Sample_Barcode")) %>%
  mutate(MUT=ifelse(MUT %in% 1,1,0)) %>%
  select(which(y==tme$X1)+1,MUT)
ggplot(dt,aes(as.factor(MUT),dt[,1])) + geom_violin()

summary(lm(X41~MUT,dt))

names(s1)
names(s2)


lapply(split(mutations1$Entrez_Gene_Id))


intersect(names(s1),names(s2))



###Read in clinical
covar=read_tsv("TME/data_clinical_patient.txt",comment='#')
covar=read_tsv("TME/data_clinical_sample.txt",comment = "#") %>% 
  full_join(.,covar,by=c("PATIENT_ID")) %>%
  mutate(Gleason=ifelse(GLEASON_SCORE=="UNK",NA,as.numeric(GLEASON_SCORE)))

covar_surv=covar %>% filter(!is.na(OS_STATUS)) %>% 
  mutate(SAMPLE_ID=gsub("-",".",SAMPLE_ID)) %>%
  mutate(SAMPLE_ID=gsub("SUTC","SU2C",SAMPLE_ID)) %>%
  select(PATIENT_ID,SAMPLE_ID,OS_STATUS,OS_MONTHS)

covar_surv$SurvObj <- with(covar_surv, Surv(OS_MONTHS, OS_STATUS == "DECEASED"))


tmb_surv=tmb %>% 
  filter(Tumor_Sample_Barcode %in% covar_surv$SAMPLE_ID) %>%
  arrange(match(Tumor_Sample_Barcode,covar_surv$SAMPLE_ID)) %>%
  mutate(TMB_MED=ifelse(n>median(n),1,0))


km.one <- survfit(SurvObj ~ 1, data = covar_surv, conf.type = "log-log")
plot(km.one)

km.tmb <- survfit(SurvObj ~ tmb_surv$TMB_MED, data = covar_surv, conf.type = "log-log")


plot(km.tmb)
res.cox1 <- coxph(SurvObj ~ tmb_surv$TMB_MED, data =  covar_surv)
res.cox1 <- coxph(SurvObj ~ tmb_surv$n, data =  covar_surv)






tmb_tme=tmb %>% 
  filter(Tumor_Sample_Barcode %in% tme_dt$ID) %>%
  arrange(match(Tumor_Sample_Barcode,covar_surv$SAMPLE_ID))

tmb_tme_res=apply(tme_dt[,-1],2,function(x){
  summary(lm(x~log2(tmb_tme$n)))$coefficients[2,4]
})


library(caret)
tmb_rf=data.frame(tme_dt[,-1],tmb=log(tmb_tme$n))
rf1=train(tmb~.,data=tmb_rf,method="rf")

names(tmb_tme_res)=tme$X1
tmb_tme_res[order(tmb_tme_res)]

ggplot(tmb_tme,aes(x=log2(n),y=tme_dt$X64)) + geom_point()

##150 with both

##213 with rrna and 271 with polya