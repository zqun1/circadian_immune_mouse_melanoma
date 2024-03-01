#!!!Coxph analysis should follow one in ten rule to avoid overfitting!!!!#
rm(list = ls())
input_dir <- "./input/"
if(!dir.exists(input_dir)){dir.create(input_dir,recursive = T)} 
output_dir <- "./output/mMel_in/DEGs/"
if(!dir.exists(output_dir)){dir.create(output_dir,recursive = T)} 
project <- 'mMel_in'
marker_dir <- "./marker/"
if(!dir.exists(marker_dir)){dir.create(marker_dir,recursive = T)} 
TCGA_dir<- "./TCGA_data/"
if(!dir.exists(TCGA_dir)){dir.create(TCGA_dir,recursive = T)} 
library(TCGAbiolinks)
library(survival)
library(survminer)
library(tidyverse)
library(DT)
library(SummarizedExperiment)
load(paste(input_dir,project, ".RData", sep = ""))

#######download TCGA data#######
# Only run once
# FPKM2TPM <- function(fpkm){
#   exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
# }
# request_cancer=c("SKCM")
# for (i in request_cancer) {
#   cancer_type=paste("TCGA",i,sep="-")
#   print(cancer_type)
#   #下载临床数据
#   clinical <- GDCquery_clinic(project = cancer_type, type = "clinical")
#   clinical <- clinical[!is.na(clinical$primary_diagnosis),]
#   write.csv(clinical,file = paste0(TCGA_dir,cancer_type,"_clinical.csv"))
# 
#   #下载rna-seq的counts数据
#   query <- GDCquery(project = cancer_type,
#                     data.category = "Transcriptome Profiling",
#                     data.type = "Gene Expression Quantification",
#                     workflow.type = "STAR - Counts")
# 
#   GDCdownload(query, method = "api") #, files.per.chunk = 100
#   expdat <- GDCprepare(query = query)
#   TPM_matrix=assays(expdat)$tpm_unstrand #assay(expdat)%>% apply(., 2, FPKM2TPM)
#   write.csv(TPM_matrix,file = paste0(TCGA_dir,cancer_type,"_TPM.csv"))
# }
#read in data
SKCM.clinical=read.csv(file = paste0(TCGA_dir,"TCGA-SKCM_clinical.csv"),row.names = 1) #377 69
SKCM.clinical$submitter_id= gsub("-",".",SKCM.clinical$submitter_id)
rownames(SKCM.clinical)=SKCM.clinical$submitter_id
SKCM.tpm=read.csv(file = paste0(TCGA_dir,"TCGA-SKCM_TPM.csv"),row.names = 1) %>% #60660 
  .[,as.numeric(substr(colnames(.),14,15)) < 10] #remove normal tissue
SKCM.tpm[1:4,1:4]

#To remove duplicated sample, e.g. "TCGA.A6.6650.01A.11R.1774.07" "TCGA.A6.6650.01B.02R.A277.07" "TCGA.A6.6650.01A.11R.A278.07":
#01B means FFPE vial-> do not use. For transcriptome choose H > R. last 7 digits mean plate and center. Use larger portion and plate number: later.
RemoveDuplicate<- function(barcode){
  keep= rep(TRUE,length(barcode))
  id= substr(barcode,1,12)
  vial= substr(barcode,16,16)
  Analyte= substr(barcode,20,20)
  plate= substr(barcode,22,25)
  pp= substr(barcode,18,25)
  dup= id[duplicated(id)] %>% unique()
  for (i in dup) {
    index= which(id==i)
    keep[index[vial[index] %in% 'B']]= FALSE
    keep[index[vial[index] %in% 'C']]= FALSE
    keep[index[Analyte[index] %in% 'T']]= FALSE
    if ('H' %in% Analyte[index]) {
      keep[index[Analyte[index]!='H']]= FALSE
    } else {
      index= index[keep[index]]
      maxi= pp[index] %>% max()
      index2= index[pp[index]!= maxi]
      keep[index2]=FALSE
    }
  }
  return(barcode[keep])
}

SKCM.tpm= SKCM.tpm %>% colnames() %>% RemoveDuplicate() %>% SKCM.tpm[,.] #60660 470
#use first 12 characters (submitter_id) as the colnames.
colnames(SKCM.tpm) %>% substr(1,12)  -> colnames(SKCM.tpm) 
SKCM.tpm= intersect(colnames(SKCM.tpm),rownames(SKCM.clinical)) %>% SKCM.tpm[,.] #60660 469
SKCM.clinical= colnames(SKCM.tpm) %>% SKCM.clinical[.,] #469 69
write.csv(SKCM.clinical,file = paste0(TCGA_dir,"TCGA_SKCM_clinical.csv"))

#Ensembl Id to gene symbol
ensembl.genes <- rownames(SKCM.tpm) %>% gsub('(.*)\\..*','\\1',.)
library(EnsDb.Hsapiens.v86) #GRCh38.p7. got problem. E.g. ENSG00000165507 converts to C10orf10, an alias of DEPP1.
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= ensembl.genes, keytype = "GENEID", columns = c("GENEID","SYMBOL")) #57116

library(AnnotationHub)
hub <- AnnotationHub()
dm <- query(hub, c("EnsDb", "sapiens")) #AH78783 AH79689 AH83216 AH89180 AH89426 | Ensembl 103 EnsDb for Homo sapiens
geneIDs2= dm[['AH109606']] #v109, only 56424 mappings
# geneIDs1= dm[['AH78783']] #v99, only 56499 mappings
# geneIDs1= dm[['AH73881']] #v97, the last version of GRCh38.p12. only 56537 mappings
geneIDs= select(geneIDs2, keys=ensembl.genes, columns=c("GENEID","SYMBOL"), keytype="GENEID") #60432
intersect(geneIDs1$GENEID,geneIDs$GENEID) %>% length() #56961

# library(biomaRt) #also has problem. Only mapped 56383/56602 ensembl_gene_ids from the data.
# ensembl= useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# geneIDs1 <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), mart = ensembl)
# geneIDs1= geneIDs1[geneIDs1$ensembl_gene_id %in% ensembl.genes,]

# library(msigdbr) #only 38905 ensembl ID mapped???
# geneIDs1= msigdbr(species = "human")
# geneIDs1= geneIDs1[geneIDs1$ensembl_gene %in% ensembl.genes,c('ensembl_gene','gene_symbol')] %>% unique()

# #only 12405 genes?
# geneIDs1= read.table('/Users/u0143346/Desktop/Peter Carmeliet/Relavant work/HCC/R HCC/Human_ENSEMBL_Gene_ID_MSigDB.v7.4.chip',sep = '\t',header = T)
# geneIDs1= geneIDs1[geneIDs1$Probe.Set.ID %in% ensembl.genes,c('Probe.Set.ID','Gene.Symbol')] %>% unique() 
ensembl.genes.version= data.frame(row.names = rownames(SKCM.tpm),ENSEMBL= ensembl.genes,SYMBOL= rownames(SKCM.tpm))
ensembl.genes.version$SYMBOL= geneIDs[match(ensembl.genes.version$ENSEMBL,geneIDs$GENEID),]$SYMBOL
SKCM.tpm$SYMBOL= ensembl.genes.version$SYMBOL

#ENSG00000157654 PALM2-AKAP2 was renamed as PALM2AKAP2
SKCM.tpm$SYMBOL= sub('PALM2AKAP2','PALM2-AKAP2',SKCM.tpm$SYMBOL)
write.csv(SKCM.tpm,file = paste0(TCGA_dir,"TCGA_SKCM_tpm.csv"))
# save(SKCM.clinical,SKCM.tpm,file = paste0(TCGA_dir,"TCGA_SKCM.RData"))
# load(file = paste0(TCGA_dir,"TCGA_SKCM.RData"))

#######DEGs between ZT01 and ZT13############
# tmp=mMel_in %>% subset(subclustering_round1.5 %in% c('CD8','CD4','T_mitotic','Treg'))
# DEGs= FindMarkers(tmp,ident.1 = 'ZT13',ident.2 = 'ZT01',group.by = 'timepoint',min.pct = 0.0001,logfc.threshold = 0)
# write.csv(DEGs,file = paste0(output_dir,project,'_abT_DEGs.csv'))
DEGs= read.csv(file = paste0(output_dir,project,'_abT_DEGs.csv'),row.names = 1)
DEGs= DEGs %>%  rownames_to_column('gene')
DEGs$cluster= 'ns'
DEGs$cluster[DEGs$avg_log2FC >0.5 &DEGs$p_val_adj <=0.05] = 'ZT13'
DEGs$cluster[DEGs$avg_log2FC < -0.5 &DEGs$p_val_adj <=0.05] = 'ZT01'

filename_prefix= 'ZT13vsZT01' #'ZT_signature' #

SKCM.tpm= read.csv(row.names = 1,file = paste0(TCGA_dir,"TCGA_SKCM_tpm.csv"))
SKCM.clinical= read.csv(row.names = 1,file = paste0(TCGA_dir,"TCGA_SKCM_clinical.csv"))

ZT13= DEGs %>% dplyr::filter(cluster=='ZT13') %>% pull(gene) %>% toupper() %>% .[.%in% SKCM.tpm$SYMBOL]
ZT01= DEGs %>% dplyr::filter(cluster=='ZT01') %>% pull(gene) %>% toupper() %>% .[.%in% SKCM.tpm$SYMBOL]
all_gene= c(ZT13,ZT01) %>% unique() #
SKCM=SKCM.tpm[SKCM.tpm$SYMBOL %in%all_gene,]
rownames(SKCM)=SKCM$SYMBOL
SKCM=SKCM[,colnames(SKCM)!='SYMBOL' & colnames(SKCM)!='GENEID']


##############Survival analysis of signature scores###########
data= cbind(SKCM.clinical,t(SKCM)%>%as.data.frame())
notDead <- is.na(data$days_to_death)
if (any(notDead == TRUE)) {
  data[notDead, "days_to_death"] <- data[notDead, "days_to_last_follow_up"]
}
data$s <- grepl("dead|deceased", data$vital_status, ignore.case = TRUE) # change alive/dead to FALSE/TRUE
data$days_to_death=as.numeric(data$days_to_death)
table(data$days_to_death %>% is.na(),data$vital_status) #still have 1 patient with unknown days_to_last_follow_up"
write.csv(data,file = paste0(TCGA_dir,'SKCM_filtered_RNAseq_clinical.csv'))


SKCM.tpm= read.csv(row.names = 1,file = paste0(TCGA_dir,"TCGA_SKCM_tpm.csv"))
SKCM.clinical= read.csv(row.names = 1,file = paste0(TCGA_dir,"TCGA_SKCM_clinical.csv"))
ZT01_ENSEMBL= SKCM.tpm %>% dplyr::filter(SYMBOL %in% ZT01) %>% rownames() #%>% c('',.) %>% data.frame(ZT01_signature=.)
ZT13_ENSEMBL= SKCM.tpm %>% dplyr::filter(SYMBOL %in% ZT13) %>% rownames() #%>% c('',.) %>% data.frame(ZT13_signature=.)
Self_geneset= list(ZT01_signature=ZT01_ENSEMBL,ZT13_signature=ZT13_ENSEMBL)
library(GSVA)
signature_results= GSVA::gsva(as.matrix(SKCM.tpm[,colnames(SKCM.tpm)!='SYMBOL']),gset.idx.list = Self_geneset,method='ssgsea',ssgsea.norm=T)

#cbind to the clinical data
data= cbind(SKCM.clinical,t(signature_results)%>%as.data.frame())
notDead <- is.na(data$days_to_death)
if (any(notDead == TRUE)) {
  data[notDead, "days_to_death"] <- data[notDead, "days_to_last_follow_up"]
}
data$s <- grepl("dead|deceased", data$vital_status, ignore.case = TRUE) # change alive/dead to FALSE/TRUE
data$days_to_death=as.numeric(data$days_to_death)
table(data$days_to_death %>% is.na(),data$vital_status) #still have 1 patient with unknown days_to_last_follow_up"

#COX analysis for each single gene
all_gene= c('ZT01_signature','ZT13_signature')
COX_single= data.frame()
for (gene in all_gene) {
  fit <- coxph(as.formula(paste('Surv(days_to_death, s)~', gene)), #Has to be an explicitly defined formula
               data=data)
  coefficients= fit %>% summary() %>% .$coefficients
  COX_single=rbind(COX_single,coefficients)
} #
colnames(COX_single)[5]='pval'
COX_single= COX_single %>% slice_min(.,order_by = pval,n=nrow(.))
tmp= COX_single[COX_single$pval<0.05 ,] #ADAMTS5 = 1.16 

threetile=function(x){ifelse(x>quantile(x,0.6666),"high",ifelse(x>quantile(x,0.3333),"medium","low"))}
twotile=function(x){ifelse(x>quantile(x,0.5),"high","low")}
data1=data
for (gene in rownames(tmp)) {
  group= paste0(gene,'_levels')
  data1[,group]<-twotile(data1[,gene])
  TCGAanalyze_survival(data1,clusterCol = group,conf.int = F,risk.table = F,legend = gene,
                       ggtheme=theme_survminer(font.x = 24,font.y = 24,font.legend = 24),color = c('#BCBBDD','#F19B9B'),
                       filename = paste0(TCGA_dir,gene,'_SKCM_survival.pdf'))
}

