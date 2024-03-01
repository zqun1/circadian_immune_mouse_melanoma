###################integrate all bulk to train the model################
rm(list = ls())
library(Seurat)
library(data.table)
library(plyr)
library(tidyverse)
library(ggplot2)
library(dittoSeq)
addTaskCallback(function(...) {set.seed(2023);TRUE}) #this will keep the seed 2023 used for the whole session! 
input_dir <- "./input/"
if(!dir.exists(input_dir)){dir.create(input_dir,recursive = T)} 
output_dir= './output/mMel_in/zeitZeiger/corrected_bulk/'
if(!dir.exists(output_dir)){dir.create(output_dir,recursive = T)} 
marker_dir <- "./marker/mMel_in/"
if(!dir.exists(marker_dir)){dir.create(marker_dir,recursive = T)} 
project <- 'mMel_in'
source('Modified function.R')
load(paste(input_dir,project, ".RData", sep = ""))

library('data.table')
library('ggplot2')
library('zeitzeiger')

doParallel::registerDoParallel(cores = 2) # register a parallel backend to minimize runtimev

#randomly subset each object to create more training and test data
dir1= '../Collaboration_Lydia/public_dermal_GSE223109/GSE223109_RAW/'
load(file = paste0(dir1,'mSkin_pub.RData'))
mSkin_pub@meta.data$timepoint= case_match(mSkin_pub@meta.data$timepoint,'ZT10'~10,'ZT14'~14,'ZT18'~18,'ZT2'~2,'ZT22'~22,'ZT6'~6)
mSkin_pub@meta.data$mouseID= mSkin_pub$sampleID

load( "../Collaboration_Lydia/Skin/input/mSkin_in.RData")
mSkin_in@meta.data$all='all'

coefficient_list= list()
error_list= list()
# load(file = paste0(output_dir,project,'_zeitZeiger_prediction_model.RData'))
# length(bulk_sample1) and length(bulk_sample2) = 14

###################prepare bulk data##############
other_cells= list(mMel_in,mSkin_in,mSkin_pub)
names(other_cells)= c('mMel_in','mSkin_in','mSkin_pub')
genes_shared= intersect(rownames(mMel_in),rownames(mSkin_in)) %>% intersect(.,rownames(mSkin_pub))
other_cells_sample= list()
time_sample= list()

# bulk skin: GSE115104
library(GEOquery)
options(timeout = max(300, getOption("timeout")))
options(download.file.method.GEOquery = "wget")
Skin_GSE115104= getGEO('GSE115104',GSEMatrix = T) 
Skin_GSE115104= pData(Skin_GSE115104[[1]])
colnames(Skin_GSE115104)[colnames(Skin_GSE115104)=='title']= 'sampleID'
Skin_GSE115104$timepoint= sub('.*timepoint(.*)re.*','\\1',Skin_GSE115104$description) %>% as.numeric()
Skin_GSE115104$timepoint=4* (Skin_GSE115104$timepoint-1)
# Skin_GSE115104= getGEOSuppFiles('GSE115104')
tmp= list()
for(i in list.files('GSE115104/',pattern = '.gz')){
  tmp[[i]] = read.table(paste0('GSE115104/',i),row.names = 1)
}
Skin_GSE115104_raw= do.call(cbind,tmp)
colnames(Skin_GSE115104_raw)= list.files('GSE115104/',pattern = '.gz') %>% str_split(.,'_',simplify = T) %>% .[,1]
Skin_GSE115104_raw= Skin_GSE115104_raw[,Skin_GSE115104$geo_accession]
#Ensembl Id to gene symbol
ensembl.genes <- rownames(Skin_GSE115104_raw) 
library(AnnotationHub)
hub <- AnnotationHub()
dm <- query(hub, c("EnsDb", "mmusculus")) #AH78783 AH79689 AH83216 AH89180 AH89426 | Ensembl 103 EnsDb for Homo sapiens
geneIDs2= dm[['AH109655']] #
geneIDs= AnnotationDbi::select(geneIDs2, keys=ensembl.genes, columns=c("GENEID","SYMBOL"), keytype="GENEID") #60432
table(duplicated(geneIDs$SYMBOL))
table(geneIDs$SYMBOL=='') #has duplicated IDs and blank IDs
Skin_GSE115104_raw$ID= geneIDs[match(rownames(Skin_GSE115104_raw),geneIDs$GENEID),]$SYMBOL
#use mean as the duplicated gene expression levels
Skin_GSE115104_raw <- Skin_GSE115104_raw%>%
  dplyr::filter(ID!='') %>% 
  dplyr::group_by(ID) %>%
  dplyr::summarise_all(mean)  %>%
  as.data.frame() %>% 
  column_to_rownames('ID')

#convert counts to tpm
id_length= lengthOf(geneIDs2, of = "gene") %>% as.data.frame()
colnames(id_length)= 'length'
geneIDs= genes(geneIDs2) %>% as.data.frame() %>% .[rownames(id_length),c('gene_id','symbol')] %>% cbind(id_length)
write.csv(geneIDs,file = 'mouse_ensembl_length.csv')
count2tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}
gene_length= geneIDs[match(rownames(Skin_GSE115104_raw),geneIDs$symbol),'length']
Skin_GSE115104_tpm= apply(Skin_GSE115104_raw,2,count2tpm,gene_length)

#select the WT mice
ids= grep('WT',Skin_GSE115104$description)
Skin_GSE115104_tpm=Skin_GSE115104_tpm[,ids]
Skin_GSE115104=Skin_GSE115104[ids,]


# bulk skin: restricted feeding
library(GEOquery)
options(timeout = max(300, getOption("timeout")))
options(download.file.method.GEOquery = "wget")
Skin_GSE83855= getGEO('GSE83855',GSEMatrix = T) 
Skin_GSE83855= pData(Skin_GSE83855[[1]])
colnames(Skin_GSE83855)[colnames(Skin_GSE83855)=='title']= 'sampleID'
Skin_GSE83855$timepoint= sub('.* ZT(.*)','\\1',Skin_GSE83855$characteristics_ch1.7) %>% str_split(.,pattern = '\\(',simplify = T) %>% .[,1] %>% as.numeric()
# Skin_GSE83855_raw= getGEOSuppFiles('GSE83855')
Skin_GSE83855_raw= read_tsv(paste0('GSE83855/GSE83855_rpkm_by_cqn.tsv.gz'))
Skin_GSE83855_raw= Skin_GSE83855_raw %>% dplyr::filter(transcript=='exon') %>% .[,-c(1,3)]%>% 
  dplyr::filter(geneSymbol!='') %>% 
  dplyr::group_by(geneSymbol) %>%
  dplyr::summarise_all(mean)  %>%
  as.data.frame() %>% 
  column_to_rownames('geneSymbol')
Skin_GSE83855_raw= Skin_GSE83855_raw[,Skin_GSE83855$sampleID]


fpkmToTpm <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
Skin_GSE83855_tpm= fpkmToTpm(Skin_GSE83855_raw)


#########################training on both pseudobulk and bulk#############
other_cells= list(mMel_in,mSkin_in,mSkin_pub)
names(other_cells)= c('mMel_in','mSkin_in','mSkin_pub')
genes_shared= intersect(rownames(mMel_in),rownames(mSkin_in)) %>% intersect(.,rownames(mSkin_pub))
other_cells_sample= list()
time_sample= list()
for(seed in 1:5){
  other_cells_sample[[seed]]=list()
  time_sample[[seed]]=list()
  addTaskCallback(function(...) {set.seed(seed+2008);TRUE}) #this will keep the seed 2023 used for the whole session! 
  for(j in names(other_cells)){
    cells_to_sample= sample(rownames(other_cells[[j]]@meta.data),round(nrow(other_cells[[j]]@meta.data)*0.6))
    other_cells_sample[[seed]][[j]]= other_cells[[j]] %>% subset(cells=cells_to_sample,features= genes_shared) %>% 
      AverageExpression(assays = 'RNA',return.seurat = F,group.by = 'mouseID') %>% 
      .$RNA %>% t()
    meta= table(other_cells[[j]]$timepoint,other_cells[[j]]$mouseID) %>% reshape2::melt() %>% dplyr::filter(value!=0) %>%.[,1:2] %>% 
      .[match(rownames(other_cells_sample[[seed]][[j]]),.$Var2),]
    colnames(meta)=c('timepoint','Cell')
    if(j!='mSkin_pub'){meta$timepoint=  case_match(meta$timepoint,'ZT01'~1,'ZT07'~7,'ZT13'~13,'ZT19'~19)}else{meta$timepoint=as.numeric(meta$timepoint)}
    time_sample[[seed]][[j]]= meta$timepoint
  }
}
pseudo_cell= other_cells_sample[1:4] %>% unlist(recursive = F) 
pseudo_cell[c(2,3,5,6,9)]=NULL #mSKin_pub has 18 samples while mSKin_in only 8. To balance the samples, mix at 1:2 ratio. melanoma is underpresented. Therefore, mix with total healthy samples at 1:2
pseudo_cell= pseudo_cell%>% do.call(rbind,.) *100 #tmp like
time= time_sample[1:4]%>% unlist(recursive = F) 
time[c(2,3,5,6,9)]= NULL
time= time %>% unlist()/24

pseudo_cell_others= other_cells_sample[[5]] %>% do.call(rbind,.)*100 #tmp like
timeTest= time_sample[[5]]%>% unlist() /24


#add in bulk data
genes_shared= intersect(rownames(mMel_in),rownames(mSkin_in)) %>% intersect(.,rownames(mSkin_pub))%>% 
  intersect(.,rownames(Skin_GSE115104_tpm))%>% intersect(.,rownames(Skin_GSE83855_tpm))
bulk_sample1= sample(1:ncol(Skin_GSE115104_tpm),round(ncol(Skin_GSE115104_tpm)*0.6))
bulk_sample2= sample(1:ncol(Skin_GSE83855_tpm),round(ncol(Skin_GSE83855_tpm)*0.6))

pseudo_cell= rbind(pseudo_cell[,genes_shared],t(Skin_GSE115104_tpm[genes_shared,bulk_sample1]),t(Skin_GSE83855_tpm[genes_shared,bulk_sample2]))
time= c(time,Skin_GSE115104$timepoint[bulk_sample1]/24,Skin_GSE83855$timepoint[bulk_sample2]%%24/24)

pseudo_cell_others= rbind(pseudo_cell_others[,genes_shared],t(Skin_GSE115104_tpm[genes_shared,-bulk_sample1]),t(Skin_GSE83855_tpm[genes_shared,-bulk_sample2]))
timeTest= c(timeTest,Skin_GSE115104$timepoint[-bulk_sample1]/24,Skin_GSE83855$timepoint[-bulk_sample2]%%24/24)

x=pseudo_cell
xTest= pseudo_cell_others
nObs_others= nrow(xTest)
nObs= nrow(x)

#correct batch effect
library(sva)
batch= c(rep('mMel_in1',8),
         rep('mMel_in4',8),
         rep('mMel_in2',8),
         rep(c('mSkin_in21','mSkin_in22'),4),
         rep('mMel_in3',8),
         rep(c('mSkin_in31','mSkin_in32'),4),
         rep('mSkin_pub3',18),
         
         rep('Skin_GSE115104_tpm',length(bulk_sample1)),
         rep('Skin_GSE83855_tpm',length(bulk_sample2)),
         
         rep('mMel_in5',8),
         rep(c('mSkin_in51','mSkin_in52'),4),
         rep('mSkin_pub5',18),
         rep('Skin_GSE115104_tpm',ncol(Skin_GSE115104_tpm)-length(bulk_sample1)),
         rep('Skin_GSE83855_tpm',ncol(Skin_GSE83855_tpm)-length(bulk_sample2)))
combat_edata = ComBat(dat= cbind(t(pseudo_cell),t(pseudo_cell_others)), batch=batch, mean.only = F) #mean.only=F work better than mean.only=T
x= t(combat_edata[,1:94])
xTest= t(combat_edata[,95:ncol(combat_edata)])


sumabsv = c(1, 1.5, 2, 3)
nSpc = 1:4
nFolds = 10
timeTrain = time
xTrain= x
foldid = sample(rep(1:nFolds, length.out = nObs))
fitResultList = zeitzeigerFitCv(x, time, foldid)
spcResultList = list()
for (ii in seq_len(length(sumabsv))) {
  spcResultList[[ii]] = zeitzeigerSpcCv(fitResultList, sumabsv = sumabsv[ii])}

predResultList = list()
for (ii in seq_len(length(sumabsv))) {
  predResultList[[ii]] = zeitzeigerPredictCv(
    x, time, foldid, spcResultList[[ii]], nSpc = nSpc)}

#cross validation results
timePredList = lapply(predResultList, function(a) a$timePred)

cvResult = data.table(
  do.call(rbind, timePredList),
  timeObs = rep(timeTrain, length(sumabsv)),
  sumabsv = rep(sumabsv, each = nObs),
  obs = rep(1:nObs, length(sumabsv)))

cvResultMelt = melt(
  cvResult, id.vars = c('obs', 'timeObs', 'sumabsv'), variable.name = 'nSpc',
  value.name = 'timePred', variable.factor = FALSE)
cvResultMelt[, nSpc := as.integer(substr(nSpc, 2, 2))]
cvResultMelt[, sumabsv := factor(sumabsv)]
cvResultMelt[, timeError := getCircDiff(timePred, timeObs)]

cvResultMeltGroup =
  cvResultMelt[, .(medae = median(abs(timeError))), by = .(sumabsv, nSpc)]
ggplot(cvResultMeltGroup) +
  geom_point(aes(x = nSpc, y = medae, shape = sumabsv, color = sumabsv), size = 2) +
  labs(x = 'Number of SPCs', y = 'Median absolute error') +
  theme_bw() + theme(legend.position = c(0.7, 0.7))

#get the best sumabsv and SPCs
SPCs_optimal= cvResultMeltGroup %>% arrange(medae) %>% pull(nSpc) %>% .[1]
sumabsv_optimal= cvResultMeltGroup %>% arrange(medae) %>% pull(sumabsv) %>% .[1] %>% as.character() %>% as.numeric()

#Train a model on the full test dataset
fitResultFinal = zeitzeigerFit(x, time)
spcResultFinal = zeitzeigerSpc(
  fitResultFinal$xFitMean, fitResultFinal$xFitResid, sumabsv = 3) #but the optimal sumabsv_optimal is 3

#test
predResult = zeitzeigerPredict(x, time, xTest, spcResultFinal, nSpc = 4) #but the optimal SPCs_optimal is 4
#alculate the difference between predicted time and observed time
dfTest = data.frame(
  timeObs = timeTest, timePred = predResult$timePred,
  timeError = getCircDiff(predResult$timePred, timeTest))

dfTest= 24*dfTest
avg_error= (sum(abs(dfTest$timeError))/nrow(dfTest)) %>% round(.,digits = 2)
p4=ggplot(dfTest,aes(x = as.factor(timeObs), y = timeError)) +
  geom_boxplot(outlier.color=NA) + #to avoid plotting outlier as dot
  geom_jitter()+
  labs(x = 'Observed time (h)', y = 'Error (h)',caption = paste0('Average error: ',avg_error,'h')) +
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text = element_text(size=32),
                   panel.border = element_blank(),axis.line = element_line(colour = "black"))
p4 #all big errors are from mMel_in
ggsave(paste0(output_dir,'_pseudobulk+bulk_prediction_batch_corrected.pdf'),p4,width = 10,height = 10)


#Plot the coefficients of the features for the SPCs
for(i in 1:length(spcResultList)){
  coefficient_list[[i]]= list()
  for(j in 1:length(spcResultList[[i]])){
    v = data.frame(spcResultList[[i]][[j]]$v[, 1:3])
    colnames(v) = c('SPC 1', 'SPC 2', 'SPC 3')
    rownames(v)= colnames(x)
    v = v[apply(v, 1, function(r) any(r != 0)), ]
    v[v == 0] = NA
    v = v[do.call(order, v), ]
    v$feature = rownames(v)
    vMelt = melt(setDT(v), id.vars = 'feature', variable.name = 'spc',
                 value.name = 'Coefficient')
    vMelt[, feature := factor(feature, rev(v$feature))]
    
    coefficient_list[[i]][[j]]= vMelt
  }
}

res= lapply(coefficient_list,rbindlist) %>% lapply(.,function(x){x %>% dplyr::filter(!is.na(Coefficient))})
res_shared_features= lapply(res, function(x){table(x$feature) %>% sort(decreasing = T) %>%.[.>1] %>% names }) #shared coefficient
genes_to_use_for_prediction= lapply(res, function(x){table(x$feature) %>%.[.>6] %>% names }) %>% unlist() %>% unique()
# genes_to_use_for_prediction_more_mMel_training= c('Nr1d2','Dbp','Ttc21b','Tomm6','Gpsm3','Per3','Arid4b','Rev1','Tef',
#                                                   'Arntl','Cul4b','Ofd1','Fancf','Cry1','D330050G23Rik','Fyb','Mapk8',
#                                                   'Fam76a','Rbm3','Prr13','Dnajc30','Brca1','Guk1')
# 
# genes_to_use_for_prediction_less_mMel_training= c("Cry1","Dbp","Nr1d2","Tef","Arfgap3","Rorc","Per3","Fam76a","Arntl",        
#                                                   "A930015D03Rik","Pafah1b3","Lonp2","Lpar6","Mthfd1l","Per2","Npas2",
#                                                   "Polb","Gja1","Lrig1","Rbfox2","Zfp78","Prr13","Gapdh","Ppia","Psmb2",
#                                                   "Klf9","Nr1d1"  ,"Cry2","Per1")

#refer to the end of this file for calculating seed-conserved coefficients
genes_to_use_for_prediction_final_less= c('Nr1d2','Ucp2','Cry1','Dbp','Arntl','Tef','Polq','Rorc','Arfgap3','Per3','Htr2b','Pdgfrl',
                                          'Lrig1','Per2','Fanci','Gja1','Fam76a','Polb','Rev1','Tppp3','Pafah1b3','Npas2','Lonp2','Scaper',
                                          'Mthfd1l','Rsad1','Guk1','Leo1','Has2','Nr1d1','Kdelr3', #herer >150
                                          'Lama3','Ift80','Col6a1','Prr13','Bscl2','Mrpl45','D330050G23Rik','Zfp493','Cstf1','Col5a2','Col1a2',
                                          'Nr1i3','Secisbp2','Rad54b','Wrap53','Ikbip','Pcolce','Steap4','Tubd1','Ceacam19','Hlf','Ino80','Ahsa1')
genes_to_use_for_prediction_final_more= c('Nr1d2','D330050G23Rik','Tomm6','Dbp','Tef','Per3','Prr13','Cstf1','Ttc21b','Arid4b','Arntl',
                                          'Rev1','Cry1','Gpsm3','Guk1','Ucp2','Rbm3','Sh2d1b1','Brca1','Snrpd3','Thop1','Enho','Orai3',
                                          'Vrk3','Polb','Tppp3','Dnajc30','Fyb','Fam76a','Rrp7a','Cacybp','Pi4k2b','Hnrnpa1','Enox1', #here >150
                                          'Per2','Ofd1','Incenp','Malat1','Pdcl3','Prrc1','Slc25a13','Fzr1','Htr2b','Rad23a','Zfp511',
                                          'Gnai2','Cul4b','Ndufaf2','Tial1')

genes_to_use_for_prediction= c(genes_to_use_for_prediction_final_less, genes_to_use_for_prediction_final_more) %>% unique()

#Plot the behavior of the SPCs over time
z = x %*% spcResultFinal$v[, 1:4]
colnames(z) = c('SPC 1', 'SPC 2', 'SPC 3', 'SPC 4')
zMelt = melt(
  data.table(z, obs = 1:nObs, Time = time, check.names = FALSE),
  id.vars = c('obs', 'Time'), variable.name = 'SPC', value.name = 'Abundance')
ggplot(zMelt) +
  facet_grid(vars(SPC), scales = 'free_y') +
  geom_point(aes(x = Time, y = Abundance), size = 2, shape = 1) + theme_bw()


#Plot the coefficients of the features for the SPCs
v = data.frame(spcResultFinal$v[, 1:3])
colnames(v) = c('SPC 1', 'SPC 2', 'SPC 3')
rownames(v)= colnames(x)
v = v[apply(v, 1, function(r) any(r != 0)), ]
v[v == 0] = NA
v = v[do.call(order, v), ]
v$feature = rownames(v)
vMelt = melt(setDT(v), id.vars = 'feature', variable.name = 'spc',
             value.name = 'Coefficient')
vMelt[, feature := factor(feature, rev(v$feature))]

ggplot(vMelt) +
  facet_wrap(vars(spc), nrow = 1) +
  geom_bar(aes(x = feature, y = Coefficient), stat = 'identity') +
  labs(x = 'Feature') + coord_flip() +
  theme_bw() + theme(panel.spacing = unit(1.2, 'lines'))



#test on pseudobulk of all cells in scRNA-seq
other_cells= list(mMel_in,mSkin_in,mSkin_pub)
pseudobulk_all= other_cells %>% lapply(.,subset,features=genes_shared)%>% lapply(.,AverageExpression,group.by='mouseID') %>% lapply(.,function(x){x$RNA *100})
pseudobulk_all= do.call(cbind,pseudobulk_all)
time_all= time[33:66]
#correct batch effect
library(sva)
batch= c(rep('mMel_in1',8),
         rep('mMel_in4',8),
         rep('mMel_in2',8),
         rep(c('mSkin_in21','mSkin_in22'),4),
         rep('mMel_in3',8),
         rep(c('mSkin_in31','mSkin_in32'),4),
         rep('mSkin_pub3',18),
         
         rep('Skin_GSE115104_tpm',length(bulk_sample1)),
         rep('Skin_GSE83855_tpm',length(bulk_sample2)),
         
         rep('mMel_in5',8),
         rep(c('mSkin_in51','mSkin_in52'),4),
         rep('mSkin_pub5',18))
combat_edata = ComBat(dat= cbind(t(pseudo_cell),pseudobulk_all), batch=batch, mean.only = F)
x= t(combat_edata[,1:94])
pseudobulk_all= t(combat_edata[,95:ncol(combat_edata)])

#Train a model on the full test dataset
fitResultFinal = zeitzeigerFit(x, time)
spcResultFinal = zeitzeigerSpc(
  fitResultFinal$xFitMean, fitResultFinal$xFitResid, sumabsv = 3) #but the optimal sumabsv_optimal is 3

#test
predResult = zeitzeigerPredict(x, time, pseudobulk_all, spcResultFinal, nSpc = 4) #but the optimal SPCs_optimal is 4
#alculate the difference between predicted time and observed time
dfTest = data.frame(
  timeObs = time_all, timePred = predResult$timePred,
  timeError = getCircDiff(predResult$timePred, time_all))

dfTest= 24*dfTest
avg_error= (sum(abs(dfTest$timeError))/nrow(dfTest)) %>% round(.,digits = 2)
p4=ggplot(dfTest,aes(x = as.factor(timeObs), y = timeError)) +
  geom_boxplot(outlier.color=NA) + #to avoid plotting outlier as dot
  geom_jitter()+
  labs(x = 'Observed time (h)', y = 'Error (h)',caption = paste0('Average error: ',avg_error,'h')) +
  theme_bw()+theme(text = element_text(size=18))
p4 
ggsave(paste0(output_dir,'_pseudobulk_of_all_cells_prediction_batch_corrected.pdf'),p4,width = 10,height = 10)


####################only use Coefficient genes for training and predicting########
#if only use genes_to_use_for_prediction_more_mMel_training, the performance for predicting subsets of mMel_in is better, 3.68h difference/8 predictions= 0.46h/prediction
#but predicting healthy skin is worse: 0.3h/2predictions= 0.15 h/prediction
#To prevent overfitting, I would suggest to use genes_to_use_for_prediction_more_mMel_training+ genes_to_use_for_prediction_less_mMel_training 

#!!!!!using only coefficients to predict testing data and pseudobulk of all scRNA cells reduced error by 0.19h/prediction.
#!!!!!However, increased the error of predicting pseudobulk of scRNA subsets by 1.18h/prediction
genes_to_use_for_prediction_final_less= c('Nr1d2','Ucp2','Cry1','Dbp','Arntl','Tef','Polq','Rorc','Arfgap3','Per3','Htr2b','Pdgfrl',
                                          'Lrig1','Per2','Fanci','Gja1','Fam76a','Polb','Rev1','Tppp3','Pafah1b3','Npas2','Lonp2','Scaper',
                                          'Mthfd1l','Rsad1','Guk1','Leo1','Has2','Nr1d1','Kdelr3'#, #herer >150
                                          # 'Lama3','Ift80','Col6a1','Prr13','Bscl2','Mrpl45','D330050G23Rik','Zfp493','Cstf1','Col5a2','Col1a2',
                                          # 'Nr1i3','Secisbp2','Rad54b','Wrap53','Ikbip','Pcolce','Steap4','Tubd1','Ceacam19','Hlf','Ino80','Ahsa1'
                                          )
genes_to_use_for_prediction_final_more= c('Nr1d2','D330050G23Rik','Tomm6','Dbp','Tef','Per3','Prr13','Cstf1','Ttc21b','Arid4b','Arntl',
                                          'Rev1','Cry1','Gpsm3','Guk1','Ucp2','Rbm3','Sh2d1b1','Brca1','Snrpd3','Thop1','Enho','Orai3',
                                          'Vrk3','Polb','Tppp3','Dnajc30','Fyb','Fam76a','Rrp7a','Cacybp','Pi4k2b','Hnrnpa1','Enox1'#, #here >150
                                          # 'Per2','Ofd1','Incenp','Malat1','Pdcl3','Prrc1','Slc25a13','Fzr1','Htr2b','Rad23a','Zfp511',
                                          # 'Gnai2','Cul4b','Ndufaf2','Tial1'
                                          )

genes_to_use_for_prediction= c(genes_to_use_for_prediction_final_less, genes_to_use_for_prediction_final_more) %>% unique()
genes_to_use_for_prediction= c(genes_to_use_for_prediction,'Klf9','Nr1d1','Cry2','Per1','Gapdh', 'Hprt1', 'Ppia', 'Psmb2') %>% unique() %>% # house keeping genes
  .[.%in% genes_shared] %>% .[.%notin%c('Brca1','Dnajc30')] #'Brca1','Dnajc30': sequence are <70% identical to human
x= pseudo_cell[,genes_to_use_for_prediction]
xTest= pseudo_cell_others[,genes_to_use_for_prediction]

#correct batch effect
library(sva)
batch= c(rep('mMel_in1',8),
         rep('mMel_in4',8),
         rep('mMel_in2',8),
         rep(c('mSkin_in21','mSkin_in22'),4),
         rep('mMel_in3',8),
         rep(c('mSkin_in31','mSkin_in32'),4),
         rep('mSkin_pub3',18),
         
         rep('Skin_GSE115104_tpm',length(bulk_sample1)),
         rep('Skin_GSE83855_tpm',length(bulk_sample2)),
         
         rep('mMel_in5',8),
         rep(c('mSkin_in51','mSkin_in52'),4),
         rep('mSkin_pub5',18),
         rep('Skin_GSE115104_tpm',ncol(Skin_GSE115104_tpm)-length(bulk_sample1)),
         rep('Skin_GSE83855_tpm',ncol(Skin_GSE83855_tpm)-length(bulk_sample2)))
combat_edata = ComBat(dat= cbind(t(x),t(xTest)), batch=batch, mean.only = F) #mean.only=F work better than mean.only=T
x= t(combat_edata[,1:94])
xTest= t(combat_edata[,95:ncol(combat_edata)])


#Train a model on the full test dataset
fitResultFinal = zeitzeigerFit(x, time)
spcResultFinal = zeitzeigerSpc(
  fitResultFinal$xFitMean, fitResultFinal$xFitResid, sumabsv = 3) #but the optimal sumabsv_optimal is 3

#test
predResult = zeitzeigerPredict(x, time, xTest, spcResultFinal, nSpc = 4) #but the optimal SPCs_optimal is 4
#alculate the difference between predicted time and observed time
dfTest = data.frame(
  timeObs = timeTest, timePred = predResult$timePred,
  timeError = getCircDiff(predResult$timePred, timeTest))

dfTest= 24*dfTest
avg_error= (sum(abs(dfTest$timeError))/nrow(dfTest)) %>% round(.,digits = 2)
p4=ggplot(dfTest,aes(x = as.factor(timeObs), y = timeError)) +
  geom_boxplot(outlier.color=NA) + #to avoid plotting outlier as dot
  geom_jitter()+
  labs(x = 'Observed time (h)', y = 'Error (h)',caption = paste0('Average error: ',avg_error,'h')) +
  theme_bw()+theme(text = element_text(size=18))
p4 #all big errors are from mMel_in
ggsave(paste0(output_dir,'_pseudobulk+bulk_prediction_batch_corrected_only_coEfficient_genes.pdf'),p4,width = 10,height = 10)


#Plot the behavior of the SPCs over time
z = x %*% spcResultFinal$v[, 1:3]
colnames(z) = c('SPC 1', 'SPC 2', 'SPC 3')
zMelt = melt(
  data.table(z, obs = 1:nObs, Time = time, check.names = FALSE),
  id.vars = c('obs', 'Time'), variable.name = 'SPC', value.name = 'Abundance')
ggplot(zMelt) +
  facet_grid(vars(SPC), scales = 'free_y') +
  geom_point(aes(x = Time, y = Abundance), size = 2, shape = 1) + theme_bw()


#Plot the coefficients of the features for the SPCs
v = data.frame(spcResultFinal$v[, 1:3])
colnames(v) = c('SPC 1', 'SPC 2', 'SPC 3')
rownames(v)= colnames(x)
v = v[apply(v, 1, function(r) any(r != 0)), ]
v[v == 0] = NA
v = v[do.call(order, v), ]
v$feature = rownames(v)
vMelt = melt(setDT(v), id.vars = 'feature', variable.name = 'spc',
             value.name = 'Coefficient')
vMelt[, feature := factor(feature, rev(v$feature))]

ggplot(vMelt) +
  facet_wrap(vars(spc), nrow = 1) +
  geom_bar(aes(x = feature, y = Coefficient), stat = 'identity') +
  labs(x = 'Feature') + coord_flip() +
  theme_bw() + theme(panel.spacing = unit(1.2, 'lines'))



#test on pseudobulk of all cells in scRNA-seq
other_cells= list(mMel_in,mSkin_in,mSkin_pub)
pseudobulk_all= other_cells %>% lapply(.,subset,features=genes_to_use_for_prediction)%>% lapply(.,AverageExpression,group.by='mouseID') %>% lapply(.,function(x){x$RNA *100})
pseudobulk_all= do.call(cbind,pseudobulk_all)
time_all= time[33:66]
x= pseudo_cell[,genes_to_use_for_prediction]
#correct batch effect
library(sva)
batch= c(rep('mMel_in1',8),
         rep('mMel_in4',8),
         rep('mMel_in2',8),
         rep(c('mSkin_in21','mSkin_in22'),4),
         rep('mMel_in3',8),
         rep(c('mSkin_in31','mSkin_in32'),4),
         rep('mSkin_pub3',18),
         
         rep('Skin_GSE115104_tpm',length(bulk_sample1)),
         rep('Skin_GSE83855_tpm',length(bulk_sample2)),
         
         rep('mMel_in5',8),
         rep(c('mSkin_in51','mSkin_in52'),4),
         rep('mSkin_pub5',18))
combat_edata = ComBat(dat= cbind(t(x),pseudobulk_all), batch=batch, mean.only = F)
x= t(combat_edata[,1:94])
pseudobulk_all= t(combat_edata[,95:ncol(combat_edata)])

#Train a model on the full test dataset
fitResultFinal = zeitzeigerFit(x, time)
spcResultFinal = zeitzeigerSpc(
  fitResultFinal$xFitMean, fitResultFinal$xFitResid, sumabsv = 3) #but the optimal sumabsv_optimal is 3

#test
predResult = zeitzeigerPredict(x, time, pseudobulk_all, spcResultFinal, nSpc = 4) #but the optimal SPCs_optimal is 4
#alculate the difference between predicted time and observed time
dfTest = data.frame(
  timeObs = time_all, timePred = predResult$timePred,
  timeError = getCircDiff(predResult$timePred, time_all))

dfTest= 24*dfTest
avg_error= (sum(abs(dfTest$timeError))/nrow(dfTest)) %>% round(.,digits = 2)
p4=ggplot(dfTest,aes(x = as.factor(timeObs), y = timeError)) +
  geom_boxplot(outlier.color=NA) + #to avoid plotting outlier as dot
  geom_jitter()+
  labs(x = 'Observed time (h)', y = 'Error (h)',caption = paste0('Average error: ',avg_error,'h')) +
  theme_bw()+theme(text = element_text(size=18))
p4 
ggsave(paste0(output_dir,'_pseudobulk_of_all_cells_prediction_batch_corrected_only_coEfficient_genes.pdf'),p4,width = 10,height = 10)

####################predict the time of mouse bulk skin#############
#use all cells for training model
other_cells= list(mMel_in,mSkin_in,mSkin_pub)
names(other_cells)= c('mMel_in','mSkin_in','mSkin_pub') #use all cells from mSKin_pub
genes_shared= intersect(rownames(mMel_in),rownames(mSkin_in)) %>% intersect(.,rownames(mSkin_pub))
other_cells_sample= list()
time_sample= list()
for(seed in 1:5){
  other_cells_sample[[seed]]=list()
  time_sample[[seed]]=list()
  addTaskCallback(function(...) {set.seed(seed+2008);TRUE}) #this will keep the seed 2023 used for the whole session! 
  for(j in names(other_cells)){
    cells_to_sample= sample(rownames(other_cells[[j]]@meta.data),round(nrow(other_cells[[j]]@meta.data)*0.75)) #use higher percentage than before
    other_cells_sample[[seed]][[j]]= other_cells[[j]] %>% subset(cells=cells_to_sample,features= genes_shared) %>% 
      AverageExpression(assays = 'RNA',return.seurat = F,group.by = 'mouseID') %>% 
      .$RNA %>% t()
    meta= table(other_cells[[j]]$timepoint,other_cells[[j]]$mouseID) %>% reshape2::melt() %>% dplyr::filter(value!=0) %>%.[,1:2] %>% 
      .[match(rownames(other_cells_sample[[seed]][[j]]),.$Var2),]
    colnames(meta)=c('timepoint','Cell')
    if(j!='mSkin_pub'){meta$timepoint=  case_match(meta$timepoint,'ZT01'~1,'ZT07'~7,'ZT13'~13,'ZT19'~19)}else{meta$timepoint=as.numeric(meta$timepoint)}
    time_sample[[seed]][[j]]= meta$timepoint
  }
}
pseudo_cell= other_cells_sample[1:4] %>% unlist(recursive = F) 
pseudo_cell[c(2,3,5,6,9)]=NULL #mSKin_pub has 18 samples while mSKin_in only 8. To balance the samples, mix at 1:2 ratio. melanoma is underpresented. Therefore, mix with total healthy samples at 1:2
pseudo_cell[['mSkin_pub']]= mSkin_pub %>% subset(features= genes_shared) %>%  #use all cells from mSkin_pub
  AverageExpression(assays = 'RNA',return.seurat = F,group.by = 'mouseID') %>% 
  .$RNA %>% t()

pseudo_cell= pseudo_cell%>% do.call(rbind,.) *100 #tmp like
time= time_sample[1:4]%>% unlist(recursive = F) 
time[c(2,3,5,6,9)]= NULL
time= time %>% unlist()/24

#add in bulk data
genes_shared= intersect(rownames(mMel_in),rownames(mSkin_in)) %>% intersect(.,rownames(mSkin_pub))%>% 
  intersect(.,rownames(Skin_GSE115104_tpm))%>% intersect(.,rownames(Skin_GSE83855_tpm))
bulk_sample1= sample(1:ncol(Skin_GSE115104_tpm),round(ncol(Skin_GSE115104_tpm)*0.6))
bulk_sample2= sample(1:ncol(Skin_GSE83855_tpm),round(ncol(Skin_GSE83855_tpm)*0.6))

pseudo_cell= rbind(pseudo_cell[,genes_shared],t(Skin_GSE115104_tpm[genes_shared,bulk_sample1]),t(Skin_GSE83855_tpm[genes_shared,bulk_sample2]))
time= c(time,Skin_GSE115104$timepoint[bulk_sample1]/24,Skin_GSE83855$timepoint[bulk_sample2]%%24/24)
save(pseudo_cell,time,genes_to_use_for_prediction,file = paste0(output_dir,project,'_zeitZeiger_prediction_model.RData'))

dfTest_list=list()


# bulk skin: GSE114943 --> DD instead of LD as in GSE115104 (training dataset)
library(GEOquery)
Skin_GSE114943= getGEO('GSE114943',GSEMatrix = T) 
Skin_GSE114943= pData(Skin_GSE114943[[1]])
colnames(Skin_GSE114943)[colnames(Skin_GSE114943)=='title']= 'sampleID'
Skin_GSE114943$timepoint= sub('.*timepoint(.*)re.*','\\1',Skin_GSE114943$sampleID) %>% as.numeric()
Skin_GSE114943$timepoint=4* (Skin_GSE114943$timepoint-1)
# Skin_GSE114943= getGEOSuppFiles('GSE114943')
tmp= list()
for(i in list.files('GSE114943/',pattern = '.gz')){
  tmp[[i]] = read.table(paste0('GSE114943/',i),row.names = 1,header = T)
}
Skin_GSE114943_raw= do.call(cbind,tmp)
colnames(Skin_GSE114943_raw)= list.files('GSE114943/',pattern = '.gz') %>% str_split(.,'_',simplify = T) %>% .[,1]
Skin_GSE114943_raw= Skin_GSE114943_raw[,Skin_GSE114943$geo_accession]
#Ensembl Id to gene symbol
geneIDs=read.csv(file = 'mouse_ensembl_length.csv',row.names = 1)
Skin_GSE114943_raw$ID= geneIDs[match(rownames(Skin_GSE114943_raw),geneIDs$gene_id),]$symbol
#use mean as the duplicated gene expression levels
Skin_GSE114943_raw <- Skin_GSE114943_raw%>%
  dplyr::filter(ID!='') %>% 
  dplyr::group_by(ID) %>%
  dplyr::summarise_all(mean)  %>%
  as.data.frame() %>% 
  column_to_rownames('ID')

#convert counts to tpm
count2tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}
gene_length= geneIDs[match(rownames(Skin_GSE114943_raw),geneIDs$symbol),'length']
Skin_GSE114943_tpm= apply(Skin_GSE114943_raw,2,count2tpm,gene_length)

pseudo_cell_others= Skin_GSE114943_tpm %>% t
timeTest= Skin_GSE114943$timepoint%%24/24

x= pseudo_cell[,genes_to_use_for_prediction %>% .[.%in%rownames(Skin_GSE114943_raw)]]
xTest= pseudo_cell_others[,genes_to_use_for_prediction %>% .[.%in%rownames(Skin_GSE114943_raw)]]
nObs_others= nrow(xTest)
nObs= nrow(x)

#correct batch effect
library(sva)
batch= c(rep('mMel_in1',8),
         rep('mMel_in4',8),
         rep('mMel_in2',8),
         rep(c('mSkin_in21','mSkin_in22'),4),
         rep('mMel_in3',8),
         rep(c('mSkin_in31','mSkin_in32'),4),
         rep('mSkin_pub3',18),
         
         rep('Skin_GSE115104_tpm',length(bulk_sample1)),
         rep('Skin_GSE83855_tpm',length(bulk_sample2)),
         
         rep('Skin_GSE114943',nObs_others/3),
         rep('Skin_GSE114943',nObs_others/3),
         rep('Skin_GSE1149432',nObs_others/3)
)
combat_edata = ComBat(dat= cbind(t(x),t(xTest)), batch=batch, mean.only = F) #mean.only=F work better than mean.only=T
x= t(combat_edata[,1:94])
xTest= t(combat_edata[,95:ncol(combat_edata)])

#Train a model on the full test dataset
fitResultFinal = zeitzeigerFit(x, time)
spcResultFinal = zeitzeigerSpc(
  fitResultFinal$xFitMean, fitResultFinal$xFitResid, sumabsv = 3) #but the optimal sumabsv_optimal is 3

#test
predResult = zeitzeigerPredict(x, time, xTest, spcResultFinal, nSpc = 4) #but the optimal SPCs_optimal is 4
#alculate the difference between predicted time and observed time
dfTest = data.frame(
  timeObs = timeTest, timePred = predResult$timePred,
  timeError = getCircDiff(predResult$timePred, timeTest))

dfTest= 24*dfTest
dfTest$genotype= Skin_GSE114943$`genotype/variation:ch1`
avg_error= (sum(abs(dfTest$timeError))/nrow(dfTest)) %>% round(.,digits = 2)
p4=ggplot(dfTest,aes(x = as.factor(timeObs), y = timeError)) +
  geom_boxplot(outlier.color=NA) + #to avoid plotting outlier as dot
  geom_jitter(aes(color= genotype))+
  labs(x = 'Observed time (h)', y = 'Error (h)',caption = paste0('Average error: ',avg_error,'h')) +
  theme_bw()+theme(text = element_text(size=18))
p4 #all big errors are from mMel_in
ggsave(paste0(output_dir,'_Skin_GSE114943_prediction_batch_corrected_only_coEfficient_genes.pdf'),p4,width = 10,height = 10)
dfTest_list[['Skin_GSE114943']]= dfTest




#peritoneal Mac in young and old age
library(GEOquery)
Mac_age= getGEO('GSE128830',GSEMatrix = T) 
Mac_age= pData(Mac_age[[1]])
colnames(Mac_age)[colnames(Mac_age)=='title']= 'sampleID'
Mac_age$timepoint= sub('(.*)_.*','\\1',Mac_age$sampleID) %>% gsub('Young|Old','',.) %>% as.numeric()
Mac_age$age= sub('(.*)([[:digit:]]+)_.*','\\1',Mac_age$sampleID) %>% sub('(.*)([[:digit:]]+)','\\1',.) 
# Mac_age_raw= getGEOSuppFiles('GSE128830')
Mac_age_young= read.table('GSE128830/RNA_young_counts_n47643x21.txt',row.names = 1)
Mac_age_old= read.table('GSE128830/RNA_old_counts_n47643x21.txt',row.names = 1) 
Mac_age_raw= cbind(Mac_age_young,Mac_age_old)
#Ensembl Id to gene symbol
geneIDs= read.csv(file = 'mouse_ensembl_length.csv',row.names = 1)
count2tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}
Mac_age_raw= Mac_age_raw[rownames(Mac_age_raw) %in% geneIDs$gene_id,]
gene_length= geneIDs[match(rownames(Mac_age_raw),geneIDs$gene_id),'length']
Mac_age_tpm= apply(Mac_age_raw,2,count2tpm,gene_length) %>% as.data.frame()
Mac_age_tpm$ID= geneIDs[match(rownames(Mac_age_tpm),geneIDs$gene_id),]$symbol
#use mean as the duplicated gene expression levels
pseudo_cell_others <- Mac_age_tpm%>%
  dplyr::filter(ID!='') %>% 
  dplyr::group_by(ID) %>%
  dplyr::summarise_all(mean)  %>%
  as.data.frame() %>% 
  column_to_rownames('ID') %>%t()
timeTest= Mac_age$timepoint%%24/24

x= pseudo_cell[,genes_to_use_for_prediction %>% .[.%in%Mac_age_tpm$ID]]
xTest= pseudo_cell_others[,genes_to_use_for_prediction %>% .[.%in%Mac_age_tpm$ID]]
nObs_others= nrow(xTest)
nObs= nrow(x)

#correct batch effect
library(sva)
batch= c(rep('mMel_in1',8),
         rep('mMel_in4',8),
         rep('mMel_in2',8),
         rep(c('mSkin_in21','mSkin_in22'),4),
         rep('mMel_in3',8),
         rep(c('mSkin_in31','mSkin_in32'),4),
         rep('mSkin_pub3',18),
         
         rep('Skin_GSE115104_tpm',length(bulk_sample1)),
         rep('Skin_GSE83855_tpm',length(bulk_sample2)),
         
         rep('Mac_age',nObs_others/2),
         rep('Mac_age3',nObs_others/2)
)
combat_edata = ComBat(dat= cbind(t(x),t(xTest)), batch=batch, mean.only = F) #mean.only=F work better than mean.only=T
x= t(combat_edata[,1:94])
xTest= t(combat_edata[,95:ncol(combat_edata)])

#Train a model on the full test dataset
fitResultFinal = zeitzeigerFit(x, time)
spcResultFinal = zeitzeigerSpc(
  fitResultFinal$xFitMean, fitResultFinal$xFitResid, sumabsv = 3) #but the optimal sumabsv_optimal is 3

#test
predResult = zeitzeigerPredict(x, time, xTest, spcResultFinal, nSpc = 4) #but the optimal SPCs_optimal is 4
#alculate the difference between predicted time and observed time
dfTest = data.frame(
  timeObs = timeTest, timePred = predResult$timePred,
  timeError = getCircDiff(predResult$timePred, timeTest))
dfTest= 24*dfTest
dfTest$age= Mac_age$age
avg_error= (sum(abs(dfTest$timeError))/nrow(dfTest)) %>% round(.,digits = 2)
p4=ggplot(dfTest,aes(x = as.factor(timeObs), y = timeError)) +
  geom_boxplot(outlier.color=NA) + #to avoid plotting outlier as dot
  geom_jitter(aes(color=age))+
  labs(x = 'Observed time (h)', y = 'Error (h)',caption = paste0('Average error: ',avg_error,'h')) +
  theme_bw()+theme(text = element_text(size=18))
p4 #all big errors are from mMel_in
ggsave(paste0(output_dir,'_Mac_age_prediction_batch_corrected_only_coEfficient_genes.pdf'),p4,width = 10,height = 10)
dfTest_list[['Mac_age']]=dfTest



#BMDM bulk
library(GEOquery)
options(timeout = max(300, getOption("timeout")))
options(download.file.method.GEOquery = "wget")
BMDM= getGEO('GSE157878',GSEMatrix = T) 
BMDM= pData(BMDM[[1]])
colnames(BMDM)[colnames(BMDM)=='title']= 'sampleID'
BMDM$timepoint= sub('.*time point (.*)','\\1',BMDM$sampleID) %>% as.numeric()
BMDM$timepoint= BMDM$timepoint -12 #according to the paper, 16h post serum shock equals to ZT4 of Mac 
# BMDM_raw= getGEOSuppFiles('GSE157878')
tmp= list()
for(i in list.files('GSE157878/',pattern = '.gz')){
  tmp[[i]]= read.table(paste0('GSE157878/',i),row.names = 1)
}
BMDM_raw= do.call(cbind,tmp)
colnames(BMDM_raw)= list.files('GSE157878/',pattern = '.gz') %>% str_split(.,'_',simplify = T) %>% .[,1]
BMDM_raw= BMDM_raw[,BMDM$geo_accession]
geneIDs= read.csv(file = 'mouse_ensembl_length.csv',row.names = 1)
count2tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}
BMDM_raw= BMDM_raw[rownames(BMDM_raw) %in% geneIDs$symbol,]
gene_length= geneIDs[match(rownames(BMDM_raw),geneIDs$symbol),'length']
BMDM_tpm= apply(BMDM_raw,2,count2tpm,gene_length) %>% t()
pseudo_cell_others= BMDM_tpm
timeTest= BMDM$timepoint%%24/24

x= pseudo_cell[,genes_to_use_for_prediction %>% .[.%in%rownames(BMDM_raw)]]
xTest= pseudo_cell_others[,genes_to_use_for_prediction %>% .[.%in%rownames(BMDM_raw)]]
nObs_others= nrow(xTest)
nObs= nrow(x)

#correct batch effect
library(sva)
batch= c(rep('mMel_in1',8),
         rep('mMel_in4',8),
         rep('mMel_in2',8),
         rep(c('mSkin_in21','mSkin_in22'),4),
         rep('mMel_in3',8),
         rep(c('mSkin_in31','mSkin_in32'),4),
         rep('mSkin_pub3',18),
         
         rep('Skin_GSE115104_tpm',length(bulk_sample1)),
         rep('Skin_GSE83855_tpm',length(bulk_sample2)),
         
         rep('BMDM',nObs_others/3),
         rep('BMDM2',nObs_others/3),
         rep('BMDM3',nObs_others/3)
)
combat_edata = ComBat(dat= cbind(t(x),t(xTest)), batch=batch, mean.only = F) #mean.only=F work better than mean.only=T
x= t(combat_edata[,1:94])
xTest= t(combat_edata[,95:ncol(combat_edata)])

#Train a model on the full test dataset
fitResultFinal = zeitzeigerFit(x, time)
spcResultFinal = zeitzeigerSpc(
  fitResultFinal$xFitMean, fitResultFinal$xFitResid, sumabsv = 3) #but the optimal sumabsv_optimal is 3

#test
predResult = zeitzeigerPredict(x, time, xTest, spcResultFinal, nSpc = 4) #but the optimal SPCs_optimal is 4
#alculate the difference between predicted time and observed time
dfTest = data.frame(
  timeObs = timeTest, timePred = predResult$timePred,
  timeError = getCircDiff(predResult$timePred, timeTest))

dfTest= 24*dfTest
avg_error= (sum(abs(dfTest$timeError))/nrow(dfTest)) %>% round(.,digits = 2)
p4=ggplot(dfTest,aes(x = as.factor(timeObs), y = timeError)) +
  geom_boxplot(outlier.color=NA) + #to avoid plotting outlier as dot
  geom_jitter()+
  labs(x = 'Observed time (h)', y = 'Error (h)',caption = paste0('Average error: ',avg_error,'h')) +
  theme_bw()+theme(text = element_text(size=18))
p4 #all big errors are from mMel_in
ggsave(paste0(output_dir,'_BMDM_prediction_batch_corrected_only_coEfficient_genes.pdf'),p4,width = 10,height = 10)
dfTest_list[['BMDM']]=dfTest




#another: mSkin microarray data:GSE174155: the data fetched by GEOquery used MAS5 to normalized data, which is not log2 transformed. RMA normalization seems to generate better results.
library(GEOquery)
Skin_microarray_GSE174155= getGEO('GSE174155',GSEMatrix = T) 
if (length(Skin_microarray_GSE174155) > 1) idx <- grep("GPL6246", attr(Skin_microarray_GSE174155, "names")) else idx <- 1
Skin_microarray_GSE174155 <- Skin_microarray_GSE174155[[idx]]
# Skin_microarray_GSE174155_raw= exprs(Skin_microarray_GSE174155) %>% as.data.frame()
Skin_microarray_GSE174155= pData(Skin_microarray_GSE174155)
library(affy)
Skin_microarray_GSE174155_raw= read.affybatch(paste0('GSE174155_RAW/',list.celfiles('GSE174155_RAW/')))
Skin_microarray_GSE174155_raw= rma(Skin_microarray_GSE174155_raw) %>% exprs() %>% as.data.frame()#normalization
#convert ID to symbol
Skin_microarray_GSE17415_soft= read.delim('GSE174155_RAW/GSE174155_family.soft') %>% splitstackshape::cSplit('Gene.Symbol','//')
Skin_microarray_GSE174155_raw$ID= Skin_microarray_GSE17415_soft[match(rownames(Skin_microarray_GSE174155_raw),Skin_microarray_GSE17415_soft$ID),]$Gene.Symbol_001

pseudo_cell_others <- Skin_microarray_GSE174155_raw %>% 
  dplyr::filter(!is.na(ID)) %>%
  dplyr::group_by(ID) %>%
  dplyr::summarise_all(mean)  %>% #use mean value as the expression for duplicated genes
  as.data.frame()%>% column_to_rownames('ID') %>% t()
timeTest=  Skin_microarray_GSE174155$characteristics_ch1.2 %>% sub('.*CT(.*)','\\1',.) %>% as.numeric() %%24/24

x= pseudo_cell[,genes_to_use_for_prediction %>% .[.%in%Skin_microarray_GSE174155_raw$ID]]
xTest= pseudo_cell_others[,genes_to_use_for_prediction %>% .[.%in%Skin_microarray_GSE174155_raw$ID]]
nObs_others= nrow(xTest)
nObs= nrow(x)

#correct batch effect
library(sva)
batch= c(rep('mMel_in1',8),
         rep('mMel_in4',8),
         rep('mMel_in2',8),
         rep(c('mSkin_in21','mSkin_in22'),4),
         rep('mMel_in3',8),
         rep(c('mSkin_in31','mSkin_in32'),4),
         rep('mSkin_pub3',18),
         
         rep('Skin_GSE115104_tpm',length(bulk_sample1)),
         rep('Skin_GSE83855_tpm',length(bulk_sample2)),
         
         rep('KO',nObs_others/2),
         rep('WT',nObs_others/2)
)
combat_edata = ComBat(dat= cbind(t(x),t(xTest)), batch=batch, mean.only = F) #mean.only=F work better than mean.only=T
x= t(combat_edata[,1:94])
xTest= t(combat_edata[,95:ncol(combat_edata)])

#Train a model on the full test dataset
fitResultFinal = zeitzeigerFit(x, time)
spcResultFinal = zeitzeigerSpc(
  fitResultFinal$xFitMean, fitResultFinal$xFitResid, sumabsv = 3) #but the optimal sumabsv_optimal is 3

#test
predResult = zeitzeigerPredict(x, time, xTest, spcResultFinal, nSpc = 4) #but the optimal SPCs_optimal is 4
#alculate the difference between predicted time and observed time
dfTest = data.frame(
  timeObs = timeTest, timePred = predResult$timePred,
  timeError = getCircDiff(predResult$timePred, timeTest))

dfTest= 24*dfTest
dfTest$genotype= Skin_microarray_GSE174155$`genotype:ch1`
avg_error= (sum(abs(dfTest$timeError))/nrow(dfTest)) %>% round(.,digits = 2)
p4=ggplot(dfTest,aes(x = as.factor(timeObs), y = timeError)) +
  geom_boxplot(outlier.color=NA) + #to avoid plotting outlier as dot
  geom_jitter(aes(color=genotype))+
  labs(x = 'Observed time (h)', y = 'Error (h)',caption = paste0('Average error: ',avg_error,'h')) +
  theme_bw()+theme(text = element_text(size=18))
p4 #all big errors are from mMel_in
ggsave(paste0(output_dir,'_Skin_microarray_GSE174155_prediction_batch_corrected_only_coEfficient_genes.pdf'),p4,width = 10,height = 10)
dfTest_list[['Skin_microarray_GSE174155']]=dfTest




#mSkin microarray data: the data fetched by GEOquery works similarly to rma normalized raw data
library(GEOquery)
Skin_microarray= getGEO('GSE38622',GSEMatrix = T) 
if (length(Skin_microarray) > 1) idx <- grep("GPL6246", attr(Skin_microarray, "names")) else idx <- 1
Skin_microarray <- Skin_microarray[[idx]]
# Skin_microarray_raw= exprs(Skin_microarray) %>% as.data.frame()
Skin_microarray= pData(Skin_microarray)
library(affy)
Skin_microarray_raw= read.affybatch(paste0('GSE38622_RAW/',list.celfiles('GSE38622_RAW/')))
Skin_microarray_raw= rma(Skin_microarray_raw) %>% exprs() %>% as.data.frame() #normalization
#convert ID to symbol
require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("mmusculus_gene_ensembl", mart)
annotLookup <- getBM(
  mart=mart,
  attributes=c(
    "affy_mogene_1_0_st_v1",
    "ensembl_gene_id",
    "gene_biotype",
    "external_gene_name"
  ))
Skin_microarray_raw= Skin_microarray_raw[rownames(Skin_microarray_raw) %in% annotLookup$affy_mogene_1_0_st_v1,]
genes_bulk= annotLookup[match(rownames(Skin_microarray_raw),annotLookup$affy_mogene_1_0_st_v1),]$external_gene_name
Skin_microarray_raw$ID= genes_bulk

pseudo_cell_others <- Skin_microarray_raw %>% 
  .[stats::complete.cases(.[,-1]), ] %>%
  dplyr::group_by(ID) %>%
  dplyr::summarise_all(mean)  %>% #use mean value as the expression for duplicated genes
  as.data.frame()%>% column_to_rownames('ID') %>% t()
timeTest=  Skin_microarray$`sample collection:ch1` %>% sub('hr','',.) %>% as.numeric() %%24/24

x= pseudo_cell[,genes_to_use_for_prediction %>% .[.%in%Skin_microarray_raw$ID]]
xTest= pseudo_cell_others[,genes_to_use_for_prediction %>% .[.%in%Skin_microarray_raw$ID]]
nObs_others= nrow(xTest)
nObs= nrow(x)

#correct batch effect
library(sva)
batch= c(rep('mMel_in1',8),
         rep('mMel_in4',8),
         rep('mMel_in2',8),
         rep(c('mSkin_in21','mSkin_in22'),4),
         rep('mMel_in3',8),
         rep(c('mSkin_in31','mSkin_in32'),4),
         rep('mSkin_pub3',18),
         
         rep('Skin_GSE115104_tpm',length(bulk_sample1)),
         rep('Skin_GSE83855_tpm',length(bulk_sample2)),
         
         rep('BMDM',nObs_others)
)
combat_edata = ComBat(dat= cbind(t(x),t(xTest)), batch=batch, mean.only = F) #mean.only=F work better than mean.only=T
x= t(combat_edata[,1:94])
xTest= t(combat_edata[,95:ncol(combat_edata)])

#Train a model on the full test dataset
fitResultFinal = zeitzeigerFit(x, time)
spcResultFinal = zeitzeigerSpc(
  fitResultFinal$xFitMean, fitResultFinal$xFitResid, sumabsv = 3) #but the optimal sumabsv_optimal is 3

#test
predResult = zeitzeigerPredict(x, time, xTest, spcResultFinal, nSpc = 4) #but the optimal SPCs_optimal is 4
#alculate the difference between predicted time and observed time
dfTest = data.frame(
  timeObs = timeTest, timePred = predResult$timePred,
  timeError = getCircDiff(predResult$timePred, timeTest))

dfTest= 24*dfTest
avg_error= (sum(abs(dfTest$timeError))/nrow(dfTest)) %>% round(.,digits = 2)
p4=ggplot(dfTest,aes(x = as.factor(timeObs), y = timeError)) +
  geom_boxplot(outlier.color=NA) + #to avoid plotting outlier as dot
  geom_jitter()+
  labs(x = 'Observed time (h)', y = 'Error (h)',caption = paste0('Average error: ',avg_error,'h')) +
  theme_bw()+theme(text = element_text(size=18))
p4 #all big errors are from mMel_in
ggsave(paste0(output_dir,'_Skin_microarray_prediction_batch_corrected_only_coEfficient_genes.pdf'),p4,width = 10,height = 10)
dfTest_list[['Skin_microarray']]=dfTest



#mSkin microarray data: GSE38623 --> skin at different hair growth fate as above
library(GEOquery)
Skin_microarray_GSE38623= getGEO('GSE38623',GSEMatrix = T) 
if (length(Skin_microarray_GSE38623) > 1) idx <- grep("GPL6246", attr(Skin_microarray_GSE38623, "names")) else idx <- 1
Skin_microarray_GSE38623 <- Skin_microarray_GSE38623[[idx]]
# Skin_microarray_GSE38623_raw= exprs(Skin_microarray_GSE38623) %>% as.data.frame()
Skin_microarray_GSE38623= pData(Skin_microarray_GSE38623)
library(affy)
Skin_microarray_GSE38623_raw= read.affybatch(paste0('GSE38623_RAW/',list.celfiles('GSE38623_RAW/')))
Skin_microarray_GSE38623_raw= rma(Skin_microarray_GSE38623_raw) %>% exprs() %>% as.data.frame() #normalization
#convert ID to symbol
require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("mmusculus_gene_ensembl", mart)
annotLookup <- getBM(
  mart=mart,
  attributes=c(
    "affy_mogene_1_0_st_v1",
    "ensembl_gene_id",
    "gene_biotype",
    "external_gene_name"
  ))
Skin_microarray_GSE38623_raw= Skin_microarray_GSE38623_raw[rownames(Skin_microarray_GSE38623_raw) %in% annotLookup$affy_mogene_1_0_st_v1,]
genes_bulk= annotLookup[match(rownames(Skin_microarray_GSE38623_raw),annotLookup$affy_mogene_1_0_st_v1),]$external_gene_name
Skin_microarray_GSE38623_raw$ID= genes_bulk

pseudo_cell_others <- Skin_microarray_GSE38623_raw %>% 
  .[stats::complete.cases(.[,-1]), ] %>%
  dplyr::group_by(ID) %>%
  dplyr::summarise_all(mean)  %>% #use mean value as the expression for duplicated genes
  as.data.frame()%>% column_to_rownames('ID') %>% t()
timeTest=  Skin_microarray_GSE38623$`sample collection (zt):ch1` %>% sub('hr','',.) %>% as.numeric() %%24/24

x= pseudo_cell[,genes_to_use_for_prediction %>% .[.%in%Skin_microarray_GSE38623_raw$ID]]
xTest= pseudo_cell_others[,genes_to_use_for_prediction %>% .[.%in%Skin_microarray_GSE38623_raw$ID]]
nObs_others= nrow(xTest)
nObs= nrow(x)

#correct batch effect
library(sva)
batch= c(rep('mMel_in1',8),
         rep('mMel_in4',8),
         rep('mMel_in2',8),
         rep(c('mSkin_in21','mSkin_in22'),4),
         rep('mMel_in3',8),
         rep(c('mSkin_in31','mSkin_in32'),4),
         rep('mSkin_pub3',18),
         
         rep('Skin_GSE115104_tpm',length(bulk_sample1)),
         rep('Skin_GSE83855_tpm',length(bulk_sample2)),
         
         rep('BMDM',nObs_others)
)
combat_edata = ComBat(dat= cbind(t(x),t(xTest)), batch=batch, mean.only = F) #mean.only=F work better than mean.only=T
x= t(combat_edata[,1:94])
xTest= t(combat_edata[,95:ncol(combat_edata)])

#Train a model on the full test dataset
fitResultFinal = zeitzeigerFit(x, time)
spcResultFinal = zeitzeigerSpc(
  fitResultFinal$xFitMean, fitResultFinal$xFitResid, sumabsv = 3) #but the optimal sumabsv_optimal is 3

#test
predResult = zeitzeigerPredict(x, time, xTest, spcResultFinal, nSpc = 4) #but the optimal SPCs_optimal is 4
#alculate the difference between predicted time and observed time
dfTest = data.frame(
  timeObs = timeTest, timePred = predResult$timePred,
  timeError = getCircDiff(predResult$timePred, timeTest))

dfTest= 24*dfTest
avg_error= (sum(abs(dfTest$timeError))/nrow(dfTest)) %>% round(.,digits = 2)
p4=ggplot(dfTest,aes(x = as.factor(timeObs), y = timeError)) +
  geom_boxplot(outlier.color=NA) + #to avoid plotting outlier as dot
  geom_jitter()+
  labs(x = 'Observed time (h)', y = 'Error (h)',caption = paste0('Average error: ',avg_error,'h')) +
  theme_bw()+theme(text = element_text(size=18))
p4 #all big errors are from mMel_in
ggsave(paste0(output_dir,'_Skin_microarray_GSE38623_prediction_batch_corrected_only_coEfficient_genes.pdf'),p4,width = 10,height = 10)
dfTest_list[['Skin_microarray_GSE38623']]=dfTest

####################predict the time of human samples: sumabsv=1.5 and nSpc=2 works better for human############
# human skin microarray: GSE35635
library(GEOquery)
hSkin_microarray_GSE35635= getGEO('GSE35635',GSEMatrix = T) 
hSkin_microarray_GSE35635 <- hSkin_microarray_GSE35635[[1]]
hSkin_microarray_GSE35635_raw= exprs(hSkin_microarray_GSE35635) %>% as.data.frame()
hSkin_microarray_GSE35635= pData(hSkin_microarray_GSE35635)
colnames(hSkin_microarray_GSE35635)[colnames(hSkin_microarray_GSE35635)=='title']= 'sampleID'
hSkin_microarray_GSE35635$timepoint= hSkin_microarray_GSE35635$characteristics_ch1.1 %>% sub('.*point (.*) \\(.*','\\1',.) %>% 
  case_match(.,'0'~2.5,'1'~7.5,'2'~12.5) 
hSkin_microarray_GSE35635$patientID= hSkin_microarray_GSE35635$sampleID %>% sub('.*subject (.*)','\\1',.)

hSkin_microarray_GSE35635_soft= read.delim('GSE35635_RAW/GSE35635_family.soft') 
hSkin_microarray_GSE35635_raw$ID= hSkin_microarray_GSE35635_soft[match(rownames(hSkin_microarray_GSE35635_raw),hSkin_microarray_GSE35635_soft$ID),]$GENE_SYMBOL
hSkin_microarray_GSE35635_raw$ID= hSkin_microarray_GSE35635_raw$ID %>% capitalize_1st(.,original = 'upper') #transform to mouse gene format
pseudo_cell_others <- hSkin_microarray_GSE35635_raw  %>% 
  dplyr::filter(!is.na(ID)) %>%
  dplyr::filter(ID!='') %>%
  dplyr::group_by(ID) %>%
  dplyr::summarise_all(mean)  %>% #use mean value as the expression for duplicated genes
  as.data.frame()%>% column_to_rownames('ID') %>% t()
timeTest=  hSkin_microarray_GSE35635$timepoint %%24/24

x= pseudo_cell[,genes_to_use_for_prediction %>% .[.%in%hSkin_microarray_GSE35635_raw$ID]]
xTest= pseudo_cell_others[,genes_to_use_for_prediction %>% .[.%in%hSkin_microarray_GSE35635_raw$ID]]
nObs_others= nrow(xTest)
nObs= nrow(x)

#correct batch effect
library(sva)
batch= c(rep('mMel_in1',8),
         rep('mMel_in4',8),
         rep('mMel_in2',8),
         rep(c('mSkin_in21','mSkin_in22'),4),
         rep('mMel_in3',8),
         rep(c('mSkin_in31','mSkin_in32'),4),
         rep('mSkin_pub3',18),
         
         rep('Skin_GSE115104_tpm',length(bulk_sample1)),
         rep('Skin_GSE83855_tpm',length(bulk_sample2)),
         
         rep('hSKIN',nObs_others)
)
combat_edata = ComBat(dat= cbind(t(x),t(xTest)), batch=batch, mean.only = F) #mean.only=F work better than mean.only=T
x= t(combat_edata[,1:94])
xTest= t(combat_edata[,95:ncol(combat_edata)])

#Train a model on the full test dataset
fitResultFinal = zeitzeigerFit(x, time)
spcResultFinal = zeitzeigerSpc(
  fitResultFinal$xFitMean, fitResultFinal$xFitResid, sumabsv = 1.5) #but the optimal sumabsv_optimal is 1.5

#test
predResult = zeitzeigerPredict(x, time, xTest, spcResultFinal, nSpc = 2) #but the optimal SPCs_optimal is 2
#alculate the difference between predicted time and observed time

#to correct the human and mouse behavior
dfTest = data.frame(
  timeObs = timeTest, timePred = (predResult$timePred+0.5)%%1,
  timeError = getCircDiff((predResult$timePred+0.5)%%1, timeTest))
dfTest= 24*dfTest
avg_error= (sum(abs(dfTest$timeError))/nrow(dfTest)) %>% round(.,digits = 2)
p4=ggplot(dfTest,aes(x = as.factor(timeObs), y = timeError)) +
  geom_boxplot(outlier.color=NA) + #to avoid plotting outlier as dot
  geom_jitter()+
  labs(x = 'Observed time (h)', y = 'Error (h)',caption = paste0('Average error: ',avg_error,'h')) +
  theme_bw()+theme(text = element_text(size=18))
p4 #all big errors are from mMel_in
ggsave(paste0(output_dir,'_hSkin_microarray_GSE35635_prediction_batch_corrected_only_coEfficient_genes.pdf'),p4,width = 10,height = 10)
dfTest_list[['hSkin_microarray_GSE35635']]=dfTest





# human skin microarray: GSE205155
library(GEOquery)
hSkin_microarray_GSE205155= getGEO('GSE205155',GSEMatrix = T) 
hSkin_microarray_GSE205155 <- hSkin_microarray_GSE205155[[1]]
hSkin_microarray_GSE205155_raw= exprs(hSkin_microarray_GSE205155) %>% as.data.frame()
hSkin_microarray_GSE205155= pData(hSkin_microarray_GSE205155)
colnames(hSkin_microarray_GSE205155)[colnames(hSkin_microarray_GSE205155)=='title']= 'sampleID'
hSkin_microarray_GSE205155$timepoint= hSkin_microarray_GSE205155$`time of collection:ch1` %>% 
  case_match(.,'8:00h'~1,'12:00h'~5,'16:00h'~9,'20:00h'~13,'00:00h day after'~17,'4:00h day after'~21,'8:00h day after'~1) 
hSkin_microarray_GSE205155_soft= read.delim('GSE205155_RAW/GSE205155_family.soft') 
hSkin_microarray_GSE205155_raw$ID= hSkin_microarray_GSE205155_soft[match(rownames(hSkin_microarray_GSE205155_raw),hSkin_microarray_GSE205155_soft$ID),]$GENE_SYMBOL
hSkin_microarray_GSE205155_raw$ID= hSkin_microarray_GSE205155_raw$ID %>% capitalize_1st(.,original = 'upper') #transform to mouse gene format

pseudo_cell_others <- hSkin_microarray_GSE205155_raw  %>% 
  dplyr::filter(!is.na(ID)) %>%
  dplyr::filter(ID!='') %>%
  dplyr::group_by(ID) %>%
  dplyr::summarise_all(mean)  %>% #use mean value as the expression for duplicated genes
  as.data.frame()%>% column_to_rownames('ID') %>% t()
timeTest=  hSkin_microarray_GSE205155$timepoint %%24/24

x= pseudo_cell[,genes_to_use_for_prediction %>% .[.%in%hSkin_microarray_GSE205155_raw$ID]]
xTest= pseudo_cell_others[,genes_to_use_for_prediction %>% .[.%in%hSkin_microarray_GSE205155_raw$ID]]
nObs_others= nrow(xTest)
nObs= nrow(x)

#correct batch effect
library(sva)
batch= c(rep('mMel_in1',8),
         rep('mMel_in4',8),
         rep('mMel_in2',8),
         rep(c('mSkin_in21','mSkin_in22'),4),
         rep('mMel_in3',8),
         rep(c('mSkin_in31','mSkin_in32'),4),
         rep('mSkin_pub3',18),
         
         rep('Skin_GSE115104_tpm',length(bulk_sample1)),
         rep('Skin_GSE83855_tpm',length(bulk_sample2)),
         
         hSkin_microarray_GSE205155$`subject:ch1`
)
combat_edata = ComBat(dat= cbind(t(x),t(xTest)), batch=batch, mean.only = F) #mean.only=F work better than mean.only=T
x= t(combat_edata[,1:94])
xTest= t(combat_edata[,95:ncol(combat_edata)])

#Train a model on the full test dataset
fitResultFinal = zeitzeigerFit(x, time)
spcResultFinal = zeitzeigerSpc(
  fitResultFinal$xFitMean, fitResultFinal$xFitResid, sumabsv = 1.5) #but the optimal sumabsv_optimal is 1.5

#test
predResult = zeitzeigerPredict(x, time, xTest, spcResultFinal, nSpc = 2) #but the optimal SPCs_optimal is 2
#alculate the difference between predicted time and observed time

#to correct the human and mouse behavior
dfTest = data.frame(
  timeObs = timeTest, timePred = (predResult$timePred+0.5)%%1,
  timeError = getCircDiff((predResult$timePred+0.5)%%1, timeTest))
dfTest= 24*dfTest
dfTest$day= '1'
dfTest$day[hSkin_microarray_GSE205155$`time of collection:ch1` %>% grep('after',.)]='2'
dfTest$tissue= hSkin_microarray_GSE205155$`tissue:ch1`
dfTest$gender= hSkin_microarray_GSE205155$`Sex:ch1`
avg_error= (sum(abs(dfTest$timeError))/nrow(dfTest)) %>% round(.,digits = 2)
p4=ggplot(dfTest,aes(x = as.factor(timeObs), y = timeError)) +
  geom_boxplot(outlier.color=NA) + #to avoid plotting outlier as dot
  geom_jitter(aes(color=tissue))+
  labs(x = 'Observed time (h)', y = 'Error (h)',caption = paste0('Average error: ',avg_error,'h')) +
  theme_bw()+theme(text = element_text(size=18))
p4 #all big errors are from mMel_in
ggsave(paste0(output_dir,'_hSkin_microarray_GSE205155_prediction_batch_corrected_only_coEfficient_genes.pdf'),p4,width = 10,height = 10)
dfTest_list[['hSkin_microarray_GSE205155']]=dfTest




# human skin microarray: GSE112660
library(GEOquery)
hSkin_microarray_GSE112660= getGEO('GSE112660',GSEMatrix = T) 
hSkin_microarray_GSE112660 <- hSkin_microarray_GSE112660[[1]]
hSkin_microarray_GSE112660_raw= exprs(hSkin_microarray_GSE112660) %>% as.data.frame()
hSkin_microarray_GSE112660= pData(hSkin_microarray_GSE112660)
colnames(hSkin_microarray_GSE112660)[colnames(hSkin_microarray_GSE112660)=='title']= 'sampleID'
hSkin_microarray_GSE112660$timepoint= hSkin_microarray_GSE112660$`collection time:ch1` %>% 
  case_match(.,'0'~17,'1200'~5,'1800'~11,'600'~23,'NA'~NA) 
hSkin_microarray_GSE112660$patientID= hSkin_microarray_GSE112660$sampleID %>% sub('(.*)_M.*','\\1',.)
#filter NA samples
hSkin_microarray_GSE112660_raw= hSkin_microarray_GSE112660_raw[,!is.na(hSkin_microarray_GSE112660$timepoint)]
hSkin_microarray_GSE112660= hSkin_microarray_GSE112660[!is.na(hSkin_microarray_GSE112660$timepoint),]

hSkin_microarray_GSE112660_soft= getGEO('GPL13667')
hSkin_microarray_GSE112660_soft= hSkin_microarray_GSE112660_soft@dataTable@table %>% as.data.frame() %>% splitstackshape::cSplit('Gene Symbol','//')
hSkin_microarray_GSE112660_raw$ID= hSkin_microarray_GSE112660_soft[match(rownames(hSkin_microarray_GSE112660_raw),hSkin_microarray_GSE112660_soft$ID),]$`Gene Symbol_01`
hSkin_microarray_GSE112660_raw$ID= hSkin_microarray_GSE112660_raw$ID %>% capitalize_1st(.,original = 'upper') #transform to mouse gene format

pseudo_cell_others <- hSkin_microarray_GSE112660_raw  %>% 
  dplyr::filter(!is.na(ID)) %>%
  dplyr::filter(ID!='') %>%
  dplyr::group_by(ID) %>%
  dplyr::summarise_all(mean)  %>% #use mean value as the expression for duplicated genes
  as.data.frame()%>% column_to_rownames('ID') %>% t()
timeTest=  hSkin_microarray_GSE112660$timepoint %%24/24

x= pseudo_cell[,genes_to_use_for_prediction %>% .[.%in%hSkin_microarray_GSE112660_raw$ID]]
xTest= pseudo_cell_others[,genes_to_use_for_prediction %>% .[.%in%hSkin_microarray_GSE112660_raw$ID]]
nObs_others= nrow(xTest)
nObs= nrow(x)

#correct batch effect
library(sva)
batch= c(rep('mMel_in1',8),
         rep('mMel_in4',8),
         rep('mMel_in2',8),
         rep(c('mSkin_in21','mSkin_in22'),4),
         rep('mMel_in3',8),
         rep(c('mSkin_in31','mSkin_in32'),4),
         rep('mSkin_pub3',18),
         
         rep('Skin_GSE115104_tpm',length(bulk_sample1)),
         rep('Skin_GSE83855_tpm',length(bulk_sample2)),
         
         hSkin_microarray_GSE112660$patientID
)
combat_edata = ComBat(dat= cbind(t(x),t(xTest)), batch=batch, mean.only = F) #mean.only=F work better than mean.only=T
x= t(combat_edata[,1:94])
xTest= t(combat_edata[,95:ncol(combat_edata)])

#Train a model on the full test dataset
fitResultFinal = zeitzeigerFit(x, time)
spcResultFinal = zeitzeigerSpc(
  fitResultFinal$xFitMean, fitResultFinal$xFitResid, sumabsv = 1.5) #but the optimal sumabsv_optimal is 1.5

#test
predResult = zeitzeigerPredict(x, time, xTest, spcResultFinal, nSpc = 2) #but the optimal SPCs_optimal is 2
#alculate the difference between predicted time and observed time

#to correct the human and mouse behavior
dfTest = data.frame(
  timeObs = timeTest, timePred = (predResult$timePred+0.5)%%1,
  timeError = getCircDiff((predResult$timePred+0.5)%%1, timeTest))
dfTest= 24*dfTest
avg_error= (sum(abs(dfTest$timeError))/nrow(dfTest)) %>% round(.,digits = 2)
p4=ggplot(dfTest,aes(x = as.factor(timeObs), y = timeError)) +
  geom_boxplot(outlier.color=NA) + #to avoid plotting outlier as dot
  geom_jitter()+
  labs(x = 'Observed time (h)', y = 'Error (h)',caption = paste0('Average error: ',avg_error,'h')) +
  theme_bw()+theme(text = element_text(size=18))
p4 #
ggsave(paste0(output_dir,'_hSkin_microarray_GSE112660_prediction_batch_corrected_only_coEfficient_genes.pdf'),p4,width = 10,height = 10)
dfTest_list[['hSkin_microarray_GSE112660']]=dfTest



# human skin microarray: GSE139305 --> this dataset includes hSkin_microarray_GSE112660
library(GEOquery)
hSkin_microarray_GSE139305= getGEO('GSE139305',GSEMatrix = T) 
hSkin_microarray_GSE139305 <- hSkin_microarray_GSE139305[[1]]
hSkin_microarray_GSE139305_raw= exprs(hSkin_microarray_GSE139305) %>% as.data.frame()
hSkin_microarray_GSE139305= pData(hSkin_microarray_GSE139305)
colnames(hSkin_microarray_GSE139305)[colnames(hSkin_microarray_GSE139305)=='title']= 'sampleID'
hSkin_microarray_GSE139305$timepoint= hSkin_microarray_GSE139305$characteristics_ch1.1 %>% 
  case_match(.,'collection time: 0'~17,'collection time: 1200'~5,'collection time: 1800'~11,'collection time: 600'~23,'NA'~NA) 
hSkin_microarray_GSE139305$patientID= hSkin_microarray_GSE139305$sampleID %>% sub('(.*) M.*','\\1',.)
#filter NA samples
hSkin_microarray_GSE139305_raw= hSkin_microarray_GSE139305_raw[,!is.na(hSkin_microarray_GSE139305$timepoint)]
hSkin_microarray_GSE139305= hSkin_microarray_GSE139305[!is.na(hSkin_microarray_GSE139305$timepoint),]

hSkin_microarray_GSE139305_soft= getGEO('GPL13667')
hSkin_microarray_GSE139305_soft= hSkin_microarray_GSE139305_soft@dataTable@table %>% as.data.frame() %>% splitstackshape::cSplit('Gene Symbol','//')
hSkin_microarray_GSE139305_raw$ID= hSkin_microarray_GSE139305_soft[match(rownames(hSkin_microarray_GSE139305_raw),hSkin_microarray_GSE139305_soft$ID),]$`Gene Symbol_01`
hSkin_microarray_GSE139305_raw$ID= hSkin_microarray_GSE139305_raw$ID %>% capitalize_1st(.,original = 'upper') #transform to mouse gene format

pseudo_cell_others <- hSkin_microarray_GSE139305_raw  %>% 
  dplyr::filter(!is.na(ID)) %>%
  dplyr::filter(ID!='') %>%
  dplyr::group_by(ID) %>%
  dplyr::summarise_all(mean)  %>% #use mean value as the expression for duplicated genes
  as.data.frame()%>% column_to_rownames('ID') %>% t()
timeTest=  hSkin_microarray_GSE139305$timepoint %%24/24

x= pseudo_cell[,genes_to_use_for_prediction %>% .[.%in%hSkin_microarray_GSE139305_raw$ID]]
xTest= pseudo_cell_others[,genes_to_use_for_prediction %>% .[.%in%hSkin_microarray_GSE139305_raw$ID]]
nObs_others= nrow(xTest)
nObs= nrow(x)

#correct batch effect
library(sva)
batch= c(rep('mMel_in1',8),
         rep('mMel_in4',8),
         rep('mMel_in2',8),
         rep(c('mSkin_in21','mSkin_in22'),4),
         rep('mMel_in3',8),
         rep(c('mSkin_in31','mSkin_in32'),4),
         rep('mSkin_pub3',18),
         
         rep('Skin_GSE115104_tpm',length(bulk_sample1)),
         rep('Skin_GSE83855_tpm',length(bulk_sample2)),
         
         hSkin_microarray_GSE139305$patientID
)
combat_edata = ComBat(dat= cbind(t(x),t(xTest)), batch=batch, mean.only = F) #mean.only=F work better than mean.only=T
x= t(combat_edata[,1:94])
xTest= t(combat_edata[,95:ncol(combat_edata)])

#Train a model on the full test dataset
fitResultFinal = zeitzeigerFit(x, time)
spcResultFinal = zeitzeigerSpc(
  fitResultFinal$xFitMean, fitResultFinal$xFitResid, sumabsv = 1.5) #but the optimal sumabsv_optimal is 1.5

#test
predResult = zeitzeigerPredict(x, time, xTest, spcResultFinal, nSpc = 2) #but the optimal SPCs_optimal is 2
#alculate the difference between predicted time and observed time

#to correct the human and mouse behavior
dfTest = data.frame(
  timeObs = timeTest, timePred = (predResult$timePred+0.5)%%1,
  timeError = getCircDiff((predResult$timePred+0.5)%%1, timeTest))
dfTest= 24*dfTest
dfTest$tissue= hSkin_microarray_GSE139305$`tissue:ch1`
avg_error= (sum(abs(dfTest$timeError))/nrow(dfTest)) %>% round(.,digits = 2)
p4=ggplot(dfTest,aes(x = as.factor(timeObs), y = timeError)) +
  geom_boxplot(outlier.color=NA) + #to avoid plotting outlier as dot
  geom_jitter()+
  labs(x = 'Observed time (h)', y = 'Error (h)',caption = paste0('Average error: ',avg_error,'h')) +
  theme_bw()+theme(text = element_text(size=18))
p4 #
ggsave(paste0(output_dir,'_hSkin_microarray_GSE139305_prediction_batch_corrected_only_coEfficient_genes.pdf'),p4,width = 10,height = 10)
dfTest_list[['hSkin_microarray_GSE139305']]=dfTest


####################predict the time of human melanoma TCGA##############
TCGA_dir<- "./TCGA_data/"
if(!dir.exists(TCGA_dir)){dir.create(TCGA_dir,recursive = T)} 
library(TCGAbiolinks)
library(survival)
library(survminer)
library(tidyverse)
library(DT)
library(SummarizedExperiment)
SKCM.clinical= read.csv(file = paste0(TCGA_dir,"TCGA_SKCM_clinical_primary_tumor.csv"),row.names = 1)
SKCM.tpm= read.csv(row.names = 1,file = paste0(TCGA_dir,"TCGA_SKCM_tpm_primary_tumor.csv"))
SKCM_TCGA= SKCM.tpm %>% dplyr::filter(SYMBOL!='')
colnames(SKCM_TCGA)[colnames(SKCM_TCGA)=='SYMBOL']='ID'
SKCM_TCGA$ID= capitalize_1st(SKCM_TCGA$ID,original = 'upper')
SKCM_TCGA$ID[SKCM_TCGA$ID=='Bmal1']='Arntl' #In TCGA,there is no ARNTL but BMAL1
pseudo_cell_others <- SKCM_TCGA  %>% 
  dplyr::filter(!is.na(ID)) %>%
  dplyr::filter(ID!='') %>%
  dplyr::group_by(ID) %>%
  dplyr::summarise_all(mean)  %>% #use mean value as the expression for duplicated genes
  as.data.frame()%>% column_to_rownames('ID') %>% t()

x= pseudo_cell[,genes_to_use_for_prediction %>% .[.%in%SKCM_TCGA$ID]]
xTest= pseudo_cell_others[,genes_to_use_for_prediction %>% .[.%in%SKCM_TCGA$ID]]
nObs_others= nrow(xTest)
nObs= nrow(x)
#correct batch effect
library(sva)
batch= c(rep('mMel_in1',8),
         rep('mMel_in4',8),
         rep('mMel_in2',8),
         rep(c('mSkin_in21','mSkin_in22'),4),
         rep('mMel_in3',8),
         rep(c('mSkin_in31','mSkin_in32'),4),
         rep('mSkin_pub3',18),
         
         rep('Skin_GSE115104_tpm',length(bulk_sample1)),
         rep('Skin_GSE83855_tpm',length(bulk_sample2)),
         
         rep('hSKIN',nObs_others)
)
combat_edata = ComBat(dat= cbind(t(x),t(xTest)), batch=batch, mean.only = F) #mean.only=F work better than mean.only=T
x= t(combat_edata[,1:94])
xTest= t(combat_edata[,95:ncol(combat_edata)])

#Train a model on the full test dataset
fitResultFinal = zeitzeigerFit(x, time)
spcResultFinal = zeitzeigerSpc(
  fitResultFinal$xFitMean, fitResultFinal$xFitResid, sumabsv = 1.5) #but the optimal sumabsv_optimal is 3
#test
predResult = zeitzeigerPredict(x, time, xTest, spcResultFinal, nSpc = 2) #but the optimal SPCs_optimal is 4

#distribution
tmp= as.data.frame((predResult$timePred *24+19)%%24)
colnames(tmp)='timePred'
p4= ggplot(tmp,aes(timePred))+
  geom_density(color='blue',adjust=1/5)+ ylab('Probability')+
  geom_vline(xintercept = 13, linetype='dashed')+
  theme(text=element_text(size = 21))+
  scale_x_continuous(expand = expansion(),breaks=seq(0,24,4),limits = c(0,24)) + 
  scale_y_continuous(expand = expansion(),breaks=seq(0,0.2,0.04),limits = c(0,0.2)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p4
ggsave(filename = paste0(output_dir,project,'_TCGA_SKCM_primary_time_distribution_precise.pdf') ,p4,width = 10,height = 10)

#Plot the behavior of the SPCs over time
z = x %*% spcResultFinal$v[, 1:3]
colnames(z) = c('SPC 1', 'SPC 2', 'SPC 3')
zMelt = melt(
  data.table(z, obs = 1:nObs, Time = time, check.names = FALSE),
  id.vars = c('obs', 'Time'), variable.name = 'SPC', value.name = 'Abundance')
ggplot(zMelt) +
  facet_grid(vars(SPC), scales = 'free_y') +
  geom_point(aes(x = Time, y = Abundance), size = 2, shape = 1) + theme_bw()

#Plot the coefficients of the features for the SPCs
v = data.frame(spcResultFinal$v[, 1:3])
colnames(v) = c('SPC 1', 'SPC 2', 'SPC 3')
rownames(v)= colnames(x)
v = v[apply(v, 1, function(r) any(r != 0)), ]
v[v == 0] = NA
v = v[do.call(order, v), ]
v$feature = rownames(v)
vMelt = melt(setDT(v), id.vars = 'feature', variable.name = 'spc',
             value.name = 'Coefficient')
vMelt[, feature := factor(feature, rev(v$feature))]

ggplot(vMelt) +
  facet_wrap(vars(spc), nrow = 1) +
  geom_bar(aes(x = feature, y = Coefficient), stat = 'identity') +
  labs(x = 'Feature') + coord_flip() +
  theme_bw() + theme(panel.spacing = unit(1.2, 'lines'))

#test predicted time and survival
SKCM.clinical$predicted_time= (predResult$timePred[match(SKCM.clinical$submitter_id,colnames(SKCM_TCGA))] *24 +19) %%24
SKCM.clinical$time_group= 'unsure'
SKCM.clinical$time_group[SKCM.clinical$predicted_time>=0 & SKCM.clinical$predicted_time<=13]='morning'
SKCM.clinical$time_group[SKCM.clinical$predicted_time<=24 & SKCM.clinical$predicted_time>13]='afternoon'
SKCM.clinical$time_group= factor(SKCM.clinical$time_group,levels = c('morning','afternoon','unsure'))
data= SKCM.clinical
notDead <- is.na(data$days_to_death)
if (any(notDead == TRUE)) {
  data[notDead, "days_to_death"] <- data[notDead, "days_to_last_follow_up"]
}
data$s <- grepl("dead|deceased", data$vital_status, ignore.case = TRUE) # change alive/dead to FALSE/TRUE
data$days_to_death=as.numeric(data$days_to_death)
write.csv(SKCM.clinical,file = paste0(TCGA_dir,"TCGA_SKCM_clinical_primary_tumor_predicted_time.csv"))
# SKCM.clinical= read.csv(file = paste0(TCGA_dir,"TCGA_SKCM_clinical_primary_tumor_predicted_time.csv"),row.names = 1)


group= 'time_group'
TCGAanalyze_survival(data,clusterCol = group,conf.int = F,risk.table = F,legend = group,
                     filename = paste0(output_dir,group,'_time_group_SKCM_survival_primary_tumor_only.png'))

#test predicted time and immune infiltration
TCGA_infiltration= read.csv('infiltration_estimation_for_tcga.csv')
SKCM_infiltration= TCGA_infiltration %>% mutate(ID= sub('(.*)-[0-9+].*','\\1',cell_type)) %>% 
  mutate(ID= gsub('-','.',ID)) %>% dplyr::filter(ID %in% colnames(SKCM.tpm))
table(SKCM_infiltration$ID %>% duplicated()) #has duplicated samples from tumor --> use mean
SKCM_infiltration= SKCM_infiltration %>% .[,as.numeric(substr(.$cell_type,14,15)) < 6] %>% #6 means metastatic
  column_to_rownames('cell_type') %>%
  dplyr::group_by(ID) %>%
  dplyr::summarise_all(mean)  %>%
  as.data.frame()
SKCM_infiltration= SKCM_infiltration[match(SKCM.clinical$submitter_id,SKCM_infiltration$ID),]
SKCM_infiltration$predicted_time= (predResult$timePred[match(SKCM_infiltration$ID,colnames(SKCM_TCGA))] *24 +19) %%24 #convert to time in life
write.csv(SKCM_infiltration,file = paste0(TCGA_dir,"infiltration_estimation_for_tcga_predicted_time.csv"))
# SKCM_infiltration= read.csv(file = paste0(TCGA_dir,"infiltration_estimation_for_tcga_predicted_time.csv"),row.names = 1)


SKCM_infiltration$time_group= 'unsure'
SKCM_infiltration$time_group[SKCM_infiltration$predicted_time>0 | SKCM_infiltration$predicted_time<=13]='morning'
SKCM_infiltration$time_group[SKCM_infiltration$predicted_time<=24 & SKCM_infiltration$predicted_time>13]='afternoon'
SKCM_infiltration$time_group= factor(SKCM_infiltration$time_group,levels = c('morning','afternoon'))
suffix= 'loose_group'

SKCM_infiltration$time_group= 'unsure'
SKCM_infiltration$time_group[SKCM_infiltration$predicted_time>=3.5 & SKCM_infiltration$predicted_time<=9.5]='morning'
SKCM_infiltration$time_group[SKCM_infiltration$predicted_time<=20.5 & SKCM_infiltration$predicted_time>=15.5]='afternoon'
SKCM_infiltration$time_group= factor(SKCM_infiltration$time_group,levels = c('morning','afternoon','unsure'))
suffix= 'strict_group'


#use self-calculated CD8_ex and CD8_em for analysis
SKCM_infiltration_self= read.csv(paste0(TCGA_dir,'CIBERSORTx_Job5_Results_more_markers_for_Tcm_500_permutation.csv')) #calculated from SKCM_GSE120575_aPD1aCTLA4_prior_treatment_cell_signature
SKCM_infiltration_self= SKCM_infiltration_self[match(SKCM.clinical$submitter_id,SKCM_infiltration_self$Mixture),]
SKCM_infiltration_self$predicted_time= (predResult$timePred[match(SKCM_infiltration_self$Mixture,colnames(SKCM_TCGA))] *24 +19) %%24 #convert to time in life

SKCM_infiltration_self$time_group= 'unsure'
SKCM_infiltration_self$time_group[SKCM_infiltration_self$predicted_time>0 | SKCM_infiltration_self$predicted_time<=13]='morning'
SKCM_infiltration_self$time_group[SKCM_infiltration_self$predicted_time<=24 & SKCM_infiltration_self$predicted_time>13]='afternoon'
SKCM_infiltration_self$time_group= factor(SKCM_infiltration_self$time_group,levels = c('morning','afternoon'))
suffix= 'loose_group'


#CD8 ratio
cibersort_ratio= (SKCM_infiltration_self[,c("CD8Tex"),drop=F] %>% rowSums())/
  (SKCM_infiltration_self[,c("CD8Tcm" ,"CD8Tem" 
  ),drop=F]%>% rowSums()) #according to a review, Tfh in solid tumor is likely to be anti-tumor

data= data.frame(row.names = rownames(SKCM_infiltration_self),
                 cibersort_ratio= cibersort_ratio,
                 timePred= SKCM_infiltration_self$predicted_time,
                 time_group= SKCM_infiltration_self$time_group) 
p2=list()
for(j in c('cibersort_ratio')){
  data2= data %>% dplyr::filter(!is.nan(!!sym(j))) %>% dplyr::filter(time_group!='unsure') %>% dplyr::filter(!is.infinite(!!sym(j)))#%>% dplyr::filter(!!sym(j)<20*median(.[,j])) 
  p1=ggplot(data2,aes(time_group,!!sym(j)))+ 
    geom_boxplot(outlier.color=NA) + #to avoid plotting outlier as dot
    geom_jitter()+ #theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15))+
    stat_compare_means(method = "t.test",size=6)+
    theme(panel.background = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(legend.position= "none",text = element_text(size=18))+
    ylab('CD8_ex/CD8_non_ex')
  p2[[j]] = p1
}
p3= cowplot::plot_grid(plotlist = p2,ncol=1)
ggsave(filename = paste0(output_dir,project,'_CD8ex_CD8non_ratio_SKCM_time_',suffix,'.pdf') ,p3,width = 5,height = 5)



#calculate anti/pro lymphocyte ratio
cibersort_ratio= (SKCM_infiltration_self[,c("CD8Tex","Treg")] %>% rowSums())/
  (SKCM_infiltration_self[,c("CD8Tcm" ,"CD8Tem" ,"CD4Tconv","gdT","NK")]%>% rowSums()) #according to a review, Tfh in solid tumor is likely to be anti-tumor
data= data.frame(row.names = rownames(SKCM_infiltration_self),
                 cibersort_ratio= cibersort_ratio,
                 timePred= SKCM_infiltration_self$predicted_time,
                 time_group= SKCM_infiltration_self$time_group) 
p2=list()
for(j in c('cibersort_ratio')){
  data2= data %>% dplyr::filter(!is.nan(!!sym(j))) %>% dplyr::filter(time_group!='unsure') %>%  dplyr::filter(!is.infinite(!!sym(j)))%>% dplyr::filter(!!sym(j)<20*median(.[,j])) 
  p1=ggplot(data2,aes(time_group,!!sym(j)))+ 
    geom_boxplot(outlier.color=NA) + #to avoid plotting outlier as dot
    geom_jitter()+ #theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15))+
    stat_compare_means(method = "t.test",size=6)+
    theme(panel.background = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(legend.position= "none",text = element_text(size=18))+
    ylab('Pro/Anti-tumor lymphocyte')
  p2[[j]] = p1
}
p3= cowplot::plot_grid(plotlist = p2,ncol=1)
ggsave(filename = paste0(output_dir,project,'_anti_pro_TNK_ratio_self_SKCM_time_',suffix,'.pdf') ,p3,width = 5,height = 5)

#############summarize###########
save(dfTest_list,file = paste0(output_dir,project,'_dfTest_list.RData'))
dfTest_list$BMDM$genotype='BMDM'
dfTest_list$Mac_age$genotype='Mac'
res= dfTest_list[c(1,4:7)] %>% do.call(bind_rows,.) #did not include Mac_age or mMZ (MC38 tumro)
res= res %>% mutate(CellType= case_match(genotype, c(NA,'WT','K14-Cre (WT)') ~ 'WT',
                                         'BMDM' ~ 'BMDM',
                                         c('Bmal1-stopFL (KO)','Bmal1-stopFL/K14-Cre (RE)','Cry1(-/-) : Cry2(-/-)')~'KO'))
res$CellType= factor(res$CellType,levels=c('WT','KO','BMDM'))
avg_error= (sum(abs(res$timeError))/nrow(res)) %>% round(.,digits = 2)
tmp= res %>% dplyr::filter(CellType=='WT')
avg_error_wt= (sum(abs(tmp$timeError))/nrow(tmp)) %>% round(.,digits = 2)
tmp= res %>% dplyr::filter(CellType=='KO')
avg_error_ko= (sum(abs(tmp$timeError))/nrow(tmp)) %>% round(.,digits = 2)
tmp= res %>% dplyr::filter(CellType=='BMDM')
avg_error_bmdm= (sum(abs(tmp$timeError))/nrow(tmp)) %>% round(.,digits = 2)

p4=ggplot(res,aes(x = as.factor(timeObs), y = timeError)) +
  geom_boxplot(outlier.color=NA) + #to avoid plotting outlier as dot
  geom_jitter(aes(color=CellType))+
  labs(x = 'Observed time (h)', y = 'Error (h)',caption = paste0('Average error: ',avg_error,'h')) +
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text = element_text(size=32),
                   panel.border = element_blank(),axis.line = element_line(colour = "black"))
p4 #all big errors are from mMel_in
ggsave(paste0(output_dir,'_mouse_summary_prediction_batch_corrected_only_coEfficient_genes1.pdf'),p4,width = 10,height = 10)
p4=ggplot(res,aes(x = as.factor(CellType), y = timeError)) +
  geom_boxplot(outlier.color=NA) + #to avoid plotting outlier as dot
  geom_jitter()+
  labs(x = 'Observed time (h)', y = 'Error (h)',caption = paste0('Average error (h):\n',
                                                                 'WT     ','KO  ','BMDM\n',
                                                                 avg_error_wt,'    ',avg_error_ko,'      ',avg_error_bmdm)) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text = element_text(size=32),
        panel.border = element_blank(),axis.line = element_line(colour = "black"))
p4 #all big errors are from mMel_in
ggsave(paste0(output_dir,'_mouse_summary_prediction_batch_corrected_only_coEfficient_genes.pdf'),p4,width = 10,height = 10)


dfTest_list$hSkin_microarray_GSE35635$tissue= 'Epidermis'
res= dfTest_list[c(8:10,12)] %>% do.call(bind_rows,.) %>% dplyr::filter(tissue %notin% c('others','skin'))
res$tissue[grep('epider',res$tissue,ignore.case = T)]='Epidermis'
res$tissue[grep(' der',res$tissue,ignore.case = T)]='Dermis'
res$tissue[grep('skin_w/o',res$tissue,ignore.case = T)]='Melanoma'
res$tissue= factor(res$tissue,levels=c('Epidermis','Dermis','Melanoma'))
tmp= res %>% dplyr::filter(tissue=='Epidermis')
avg_error_epi= (sum(abs(tmp$timeError))/nrow(tmp)) %>% round(.,digits = 2)
tmp= res %>% dplyr::filter(tissue=='Dermis')
avg_error_der= (sum(abs(tmp$timeError))/nrow(tmp)) %>% round(.,digits = 2)
tmp= res %>% dplyr::filter(tissue=='Melanoma')
avg_error_mel= (sum(abs(tmp$timeError))/nrow(tmp)) %>% round(.,digits = 2)
p4=ggplot(res,aes(x = as.factor(timeObs), y = timeError)) +
  geom_boxplot(outlier.color=NA) + #to avoid plotting outlier as dot
  geom_jitter(aes(color=tissue))+
  labs(x = 'Observed time (h)', y = 'Error (h)',caption = paste0('Average error (h):\n',
                                                                 'Epiermis   Dermis  Melanoma\n',
                                                                 avg_error_epi,'          ',avg_error_der,'            ',avg_error_mel)) +
  theme_bw()+theme(text = element_text(size=18))
p4 #all big errors are from mMel_in
ggsave(paste0(output_dir,'_human_summary_prediction_batch_corrected_only_coEfficient_genes.pdf'),p4,width = 10,height = 10)
p4=ggplot(res,aes(x = as.factor(tissue), y = timeError)) +
  geom_boxplot(outlier.color=NA) + #to avoid plotting outlier as dot
  geom_jitter()+
  labs(x = 'Observed time (h)', y = 'Error (h)',caption = paste0('Average error (h):\n',
                                                                 'Epiermis   Dermis  Melanoma\n',
                                                                 avg_error_epi,'          ',avg_error_der,'            ',avg_error_mel)) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text = element_text(size=32),
        panel.border = element_blank(),axis.line = element_line(colour = "black"))
p4 #all big errors are from mMel_in
ggsave(paste0(output_dir,'_human_summary2_prediction_batch_corrected_only_coEfficient_genes.pdf'),p4,width = 10,height = 10)


dfTest_list$hSkin_microarray_GSE35635$tissue= 'Epidermis'
res= dfTest_list[c(8:10,12)] %>% do.call(bind_rows,.) %>% dplyr::filter(tissue %notin% c('others','skin'))
res$tissue[grep('der',res$tissue,ignore.case = T)]='Healthy'
res$tissue[grep('skin_w/o',res$tissue,ignore.case = T)]='Melanoma'
res$tissue= factor(res$tissue,levels=c('Healthy','Melanoma'))
tmp= res %>% dplyr::filter(tissue=='Healthy')
avg_error_epi= (sum(abs(tmp$timeError))/nrow(tmp)) %>% round(.,digits = 2)
tmp= res %>% dplyr::filter(tissue=='Melanoma')
avg_error_mel= (sum(abs(tmp$timeError))/nrow(tmp)) %>% round(.,digits = 2)
p4=ggplot(res[res$tissue=='Healthy',],aes(x = as.factor(tissue), y = timeError)) +
  geom_violin(color='#BCBBDD',fill='#A7C6E8',draw_quantiles = c(0.25, 0.5, 0.75)) + #
  geom_violin_quantiles(quantiles = c(0.25, 0.5, 0.75), linesize = c(1, 2, 1),linetype='dashed')+
  labs(x = 'Observed time (h)', y = 'Error (h)',caption = paste0('Average error (h): ',avg_error_epi)) +
  theme_bw()+theme(plot.caption = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text = element_text(size=32),
        panel.border = element_blank(),axis.line = element_line(colour = "black"))
p4 #all big errors are from mMel_in
ggsave(paste0(output_dir,'_human_summary5_prediction_batch_corrected_only_coEfficient_genes.pdf'),p4,width = 5,height = 10)
p4=ggplot(res[res$tissue=='Melanoma',],aes(x = as.factor(tissue), y = timeError)) +
  geom_boxplot(outlier.color=NA) + #to avoid plotting outlier as dot
  geom_jitter(size=6)+
  labs(x = 'Observed time (h)', y = 'Error (h)',caption = paste0('Average error (h): ',avg_error_mel)) +
  ylim(c(-12,12))+
  theme_bw()+theme(plot.caption = element_text(hjust = 0.5),)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text = element_text(size=32),
        panel.border = element_blank(),axis.line = element_line(colour = "black"))
p4 #all big errors are from mMel_in
ggsave(paste0(output_dir,'_human_summary4_prediction_batch_corrected_only_coEfficient_genes.pdf'),p4,width = 6,height = 10)


res$timeObs= (res$timeObs+6.99)%%24 #convert to o'clock
res$timePred= (res$timePred+6.99)%%24
precision= res %>% dplyr::filter((timePred>=13 & timePred<=24 & timeObs<=13) | 
                                   (timePred<=13 & timeObs>13& timeObs<23.99) )
p4=ggplot(res) +
  geom_point(aes(x = timeObs, y = timePred), size = 2, shape = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
  geom_abline(slope = 0,intercept = 13,linetype='dotted')+
  geom_vline(xintercept = 13,linetype='dotted',)+
  scale_x_continuous(expand = expansion(),breaks=seq(0,24,6)) + scale_y_continuous(expand = expansion(),breaks=seq(0,24,6)) +
  labs(x = 'Observed time (oclock)', y = 'Predicted time') + theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text = element_text(size=15))+
  annotate('rect',xmin = 0,xmax=13.01,ymin=0,ymax = 13,alpha=0.8,fill='#BCBBDD')+
  annotate('rect',xmin = 13,xmax=24,ymin=13,ymax = 24,alpha=0.8,fill='#F19B9B')+
  labs(caption = paste0('Morning/Afternoon pediction precision: ',round(1-nrow(precision)/nrow(res),digits=3)))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text = element_text(size=32))+
  theme(plot.margin = margin(1,1,0,0, "cm"))
p4
ggsave(paste0(output_dir,'_human_summary3_prediction_batch_corrected_only_coEfficient_genes.pdf'),p4,width = 10,height = 10)




########################compute the conserved Coefficients#############
coefficient_list_final_more= list()
error_list_more= list()
for(k in 1998:2018){
  coefficient_list=list()
  
  other_cells= list(mMel_in,mSkin_in,mSkin_pub)
  names(other_cells)= c('mMel_in','mSkin_in','mSkin_pub')
  genes_shared= intersect(rownames(mMel_in),rownames(mSkin_in)) %>% intersect(.,rownames(mSkin_pub))
  other_cells_sample= list()
  time_sample= list()
  for(seed in 1:5){
    other_cells_sample[[seed]]=list()
    time_sample[[seed]]=list()
    addTaskCallback(function(...) {set.seed(seed+k);TRUE}) #this will keep the seed 2023 used for the whole session! 
    for(j in names(other_cells)){
      cells_to_sample= sample(rownames(other_cells[[j]]@meta.data),round(nrow(other_cells[[j]]@meta.data)*0.6))
      other_cells_sample[[seed]][[j]]= other_cells[[j]] %>% subset(cells=cells_to_sample,features= genes_shared) %>% 
        AverageExpression(assays = 'RNA',return.seurat = F,group.by = 'mouseID') %>% 
        .$RNA %>% t()
      meta= table(other_cells[[j]]$timepoint,other_cells[[j]]$mouseID) %>% reshape2::melt() %>% dplyr::filter(value!=0) %>%.[,1:2] %>% 
        .[match(rownames(other_cells_sample[[seed]][[j]]),.$Var2),]
      colnames(meta)=c('timepoint','Cell')
      if(j!='mSkin_pub'){meta$timepoint=  case_match(meta$timepoint,'ZT01'~1,'ZT07'~7,'ZT13'~13,'ZT19'~19)}else{meta$timepoint=as.numeric(meta$timepoint)}
      time_sample[[seed]][[j]]= meta$timepoint
    }
  }
  pseudo_cell= other_cells_sample[1:4] %>% unlist(recursive = F) 
  pseudo_cell[c(2,3,5,6,9)]=NULL #mSKin_pub has 18 samples while mSKin_in only 8. To balance the samples, mix at 1:2 ratio. melanoma is underpresented. Therefore, mix with total healthy samples at 1:2
  pseudo_cell= pseudo_cell%>% do.call(rbind,.) *100 #tmp like
  time= time_sample[1:4]%>% unlist(recursive = F) 
  time[c(2,3,5,6,9)]= NULL
  time= time %>% unlist()/24
  
  pseudo_cell_others= other_cells_sample[[5]] %>% do.call(rbind,.)*100 #tmp like
  timeTest= time_sample[[5]]%>% unlist() /24
  
  
  #add in bulk data
  genes_shared= intersect(rownames(mMel_in),rownames(mSkin_in)) %>% intersect(.,rownames(mSkin_pub))%>% 
    intersect(.,rownames(Skin_GSE115104_tpm))%>% intersect(.,rownames(Skin_GSE83855_tpm))
  bulk_sample1= sample(1:ncol(Skin_GSE115104_tpm),round(ncol(Skin_GSE115104_tpm)*0.6))
  bulk_sample2= sample(1:ncol(Skin_GSE83855_tpm),round(ncol(Skin_GSE83855_tpm)*0.6))
  
  pseudo_cell= rbind(pseudo_cell[,genes_shared],t(Skin_GSE115104_tpm[genes_shared,bulk_sample1]),t(Skin_GSE83855_tpm[genes_shared,bulk_sample2]))
  time= c(time,Skin_GSE115104$timepoint[bulk_sample1]/24,Skin_GSE83855$timepoint[bulk_sample2]%%24/24)
  
  pseudo_cell_others= rbind(pseudo_cell_others[,genes_shared],t(Skin_GSE115104_tpm[genes_shared,-bulk_sample1]),t(Skin_GSE83855_tpm[genes_shared,-bulk_sample2]))
  timeTest= c(timeTest,Skin_GSE115104$timepoint[-bulk_sample1]/24,Skin_GSE83855$timepoint[-bulk_sample2]%%24/24)
  
  x=pseudo_cell
  xTest= pseudo_cell_others
  nObs_others= nrow(xTest)
  nObs= nrow(x)
  
  #correct batch effect
  library(sva)
  batch= c(rep('mMel_in1',8),
           rep('mMel_in4',8),
           rep('mMel_in2',8),
           rep(c('mSkin_in21','mSkin_in22'),4),
           rep('mMel_in3',8),
           rep(c('mSkin_in31','mSkin_in32'),4),
           rep('mSkin_pub3',18),
           
           rep('Skin_GSE115104_tpm',length(bulk_sample1)),
           rep('Skin_GSE83855_tpm',length(bulk_sample2)),
           
           rep('mMel_in5',8),
           rep(c('mSkin_in51','mSkin_in52'),4),
           rep('mSkin_pub5',18),
           rep('Skin_GSE115104_tpm',ncol(Skin_GSE115104_tpm)-length(bulk_sample1)),
           rep('Skin_GSE83855_tpm',ncol(Skin_GSE83855_tpm)-length(bulk_sample2)))
  combat_edata = ComBat(dat= cbind(t(pseudo_cell),t(pseudo_cell_others)), batch=batch, mean.only = F) #mean.only=F work better than mean.only=T
  x= t(combat_edata[,1:94])
  xTest= t(combat_edata[,95:ncol(combat_edata)])
  
  
  sumabsv = c(1, 1.5, 2, 3)
  nSpc = 1:4
  nFolds = 10
  timeTrain = time
  xTrain= x
  foldid = sample(rep(1:nFolds, length.out = nObs))
  fitResultList = zeitzeigerFitCv(x, time, foldid)
  spcResultList = list()
  for (ii in seq_len(length(sumabsv))) {
    spcResultList[[ii]] = zeitzeigerSpcCv(fitResultList, sumabsv = sumabsv[ii])}
  
  predResultList = list()
  for (ii in seq_len(length(sumabsv))) {
    predResultList[[ii]] = zeitzeigerPredictCv(
      x, time, foldid, spcResultList[[ii]], nSpc = nSpc)}
  
  #cross validation results
  timePredList = lapply(predResultList, function(a) a$timePred)
  
  cvResult = data.table(
    do.call(rbind, timePredList),
    timeObs = rep(timeTrain, length(sumabsv)),
    sumabsv = rep(sumabsv, each = nObs),
    obs = rep(1:nObs, length(sumabsv)))
  
  cvResultMelt = melt(
    cvResult, id.vars = c('obs', 'timeObs', 'sumabsv'), variable.name = 'nSpc',
    value.name = 'timePred', variable.factor = FALSE)
  cvResultMelt[, nSpc := as.integer(substr(nSpc, 2, 2))]
  cvResultMelt[, sumabsv := factor(sumabsv)]
  cvResultMelt[, timeError := getCircDiff(timePred, timeObs)]
  
  cvResultMeltGroup =
    cvResultMelt[, .(medae = median(abs(timeError))), by = .(sumabsv, nSpc)]
  ggplot(cvResultMeltGroup) +
    geom_point(aes(x = nSpc, y = medae, shape = sumabsv, color = sumabsv), size = 2) +
    labs(x = 'Number of SPCs', y = 'Median absolute error') +
    theme_bw() + theme(legend.position = c(0.7, 0.7))
  
  #get the best sumabsv and SPCs
  SPCs_optimal= cvResultMeltGroup %>% arrange(medae) %>% pull(nSpc) %>% .[1]
  sumabsv_optimal= cvResultMeltGroup %>% arrange(medae) %>% pull(sumabsv) %>% .[1] %>% as.character() %>% as.numeric()
  
  #Train a model on the full test dataset
  fitResultFinal = zeitzeigerFit(x, time)
  spcResultFinal = zeitzeigerSpc(
    fitResultFinal$xFitMean, fitResultFinal$xFitResid, sumabsv = 3) #but the optimal sumabsv_optimal is 3
  
  #test
  predResult = zeitzeigerPredict(x, time, xTest, spcResultFinal, nSpc = 4) #but the optimal SPCs_optimal is 4
  #alculate the difference between predicted time and observed time
  dfTest = data.frame(
    timeObs = timeTest, timePred = predResult$timePred,
    timeError = getCircDiff(predResult$timePred, timeTest))
  
  dfTest= 24*dfTest
  avg_error= (sum(abs(dfTest$timeError))/nrow(dfTest)) %>% round(.,digits = 2)
  error_list_more[[k]]= avg_error
  
  
  #Plot the coefficients of the features for the SPCs
  for(i in 1:length(spcResultList)){
    coefficient_list[[i]]= list()
    for(j in 1:length(spcResultList[[i]])){
      v = data.frame(spcResultList[[i]][[j]]$v[, 1:3])
      colnames(v) = c('SPC 1', 'SPC 2', 'SPC 3')
      rownames(v)= colnames(x)
      v = v[apply(v, 1, function(r) any(r != 0)), ]
      v[v == 0] = NA
      v = v[do.call(order, v), ]
      v$feature = rownames(v)
      vMelt = melt(setDT(v), id.vars = 'feature', variable.name = 'spc',
                   value.name = 'Coefficient')
      vMelt[, feature := factor(feature, rev(v$feature))]
      
      coefficient_list[[i]][[j]]= vMelt
    }
  }
  
  coefficient_list_final_more[[k]]= coefficient_list
}


coefficient_list_final_less= list()
error_list_less= list()
for(k in 1998:2018){
  coefficient_list=list()
  
  other_cells= list(mMel_in,mSkin_in,mSkin_pub)
  names(other_cells)= c('mMel_in','mSkin_in','mSkin_pub')
  genes_shared= intersect(rownames(mMel_in),rownames(mSkin_in)) %>% intersect(.,rownames(mSkin_pub))
  other_cells_sample= list()
  time_sample= list()
  for(seed in 1:3){
    other_cells_sample[[seed]]=list()
    time_sample[[seed]]=list()
    addTaskCallback(function(...) {set.seed(seed+k);TRUE}) #this will keep the seed 2023 used for the whole session! 
    for(j in names(other_cells)){
      cells_to_sample= sample(rownames(other_cells[[j]]@meta.data),round(nrow(other_cells[[j]]@meta.data)*0.6))
      other_cells_sample[[seed]][[j]]= other_cells[[j]] %>% subset(cells=cells_to_sample,features= genes_shared) %>% 
        AverageExpression(assays = 'RNA',return.seurat = F,group.by = 'mouseID') %>% 
        .$RNA %>% t()
      meta= table(other_cells[[j]]$timepoint,other_cells[[j]]$mouseID) %>% reshape2::melt() %>% dplyr::filter(value!=0) %>%.[,1:2] %>% 
        .[match(rownames(other_cells_sample[[seed]][[j]]),.$Var2),]
      colnames(meta)=c('timepoint','Cell')
      if(j!='mSkin_pub'){meta$timepoint=  case_match(meta$timepoint,'ZT01'~1,'ZT07'~7,'ZT13'~13,'ZT19'~19)}else{meta$timepoint=as.numeric(meta$timepoint)}
      time_sample[[seed]][[j]]= meta$timepoint
    }
  }
  pseudo_cell= other_cells_sample[1:2] %>% unlist(recursive = F) 
  pseudo_cell[c(6)]=NULL #mSKin_pub has 18 samples while mSKin_in only 8. To balance the samples, mix at 1:2 ratio. melanoma is underpresented. Therefore, mix with total healthy samples at 1:2
  pseudo_cell= pseudo_cell%>% do.call(rbind,.) *100 #tmp like
  time= time_sample[1:2]%>% unlist(recursive = F) 
  time[c(6)]= NULL
  time= time %>% unlist()/24
  
  pseudo_cell_others= other_cells_sample[[3]] %>% do.call(rbind,.)*100 #tmp like
  timeTest= time_sample[[3]]%>% unlist() /24
  
  
  #add in bulk data
  genes_shared= intersect(rownames(mMel_in),rownames(mSkin_in)) %>% intersect(.,rownames(mSkin_pub))%>% 
    intersect(.,rownames(Skin_GSE115104_tpm))%>% intersect(.,rownames(Skin_GSE83855_tpm))
  bulk_sample1= sample(1:ncol(Skin_GSE115104_tpm),round(ncol(Skin_GSE115104_tpm)*0.6))
  bulk_sample2= sample(1:ncol(Skin_GSE83855_tpm),round(ncol(Skin_GSE83855_tpm)*0.6))
  
  pseudo_cell= rbind(pseudo_cell[,genes_shared],t(Skin_GSE115104_tpm[genes_shared,bulk_sample1]),t(Skin_GSE83855_tpm[genes_shared,bulk_sample2]))
  time= c(time,Skin_GSE115104$timepoint[bulk_sample1]/24,Skin_GSE83855$timepoint[bulk_sample2]%%24/24)
  
  pseudo_cell_others= rbind(pseudo_cell_others[,genes_shared],t(Skin_GSE115104_tpm[genes_shared,-bulk_sample1]),t(Skin_GSE83855_tpm[genes_shared,-bulk_sample2]))
  timeTest= c(timeTest,Skin_GSE115104$timepoint[-bulk_sample1]/24,Skin_GSE83855$timepoint[-bulk_sample2]%%24/24)
  
  x=pseudo_cell
  xTest= pseudo_cell_others
  nObs_others= nrow(xTest)
  nObs= nrow(x)
  
  #correct batch effect
  library(sva)
  batch= c(rep('mMel_in2',8),
           rep(c('mSkin_in21','mSkin_in22'),4),
           rep('mMel_in3',8),
           rep(c('mSkin_in31','mSkin_in32'),4),
           rep('mSkin_pub3',18),
           
           rep('Skin_GSE115104_tpm',length(bulk_sample1)),
           rep('Skin_GSE83855_tpm',length(bulk_sample2)),
           
           rep('mMel_in5',8),
           rep(c('mSkin_in51','mSkin_in52'),4),
           rep('mSkin_pub5',18),
           rep('Skin_GSE115104_tpm',ncol(Skin_GSE115104_tpm)-length(bulk_sample1)),
           rep('Skin_GSE83855_tpm',ncol(Skin_GSE83855_tpm)-length(bulk_sample2)))
  combat_edata = ComBat(dat= cbind(t(pseudo_cell),t(pseudo_cell_others)), batch=batch, mean.only = F) #mean.only=F work better than mean.only=T
  x= t(combat_edata[,1:78])
  xTest= t(combat_edata[,79:ncol(combat_edata)])
  
  
  sumabsv = c(1, 1.5, 2, 3)
  nSpc = 1:4
  nFolds = 10
  timeTrain = time
  xTrain= x
  foldid = sample(rep(1:nFolds, length.out = nObs))
  fitResultList = zeitzeigerFitCv(x, time, foldid)
  spcResultList = list()
  for (ii in seq_len(length(sumabsv))) {
    spcResultList[[ii]] = zeitzeigerSpcCv(fitResultList, sumabsv = sumabsv[ii])}
  
  predResultList = list()
  for (ii in seq_len(length(sumabsv))) {
    predResultList[[ii]] = zeitzeigerPredictCv(
      x, time, foldid, spcResultList[[ii]], nSpc = nSpc)}
  
  #cross validation results
  timePredList = lapply(predResultList, function(a) a$timePred)
  
  cvResult = data.table(
    do.call(rbind, timePredList),
    timeObs = rep(timeTrain, length(sumabsv)),
    sumabsv = rep(sumabsv, each = nObs),
    obs = rep(1:nObs, length(sumabsv)))
  
  cvResultMelt = melt(
    cvResult, id.vars = c('obs', 'timeObs', 'sumabsv'), variable.name = 'nSpc',
    value.name = 'timePred', variable.factor = FALSE)
  cvResultMelt[, nSpc := as.integer(substr(nSpc, 2, 2))]
  cvResultMelt[, sumabsv := factor(sumabsv)]
  cvResultMelt[, timeError := getCircDiff(timePred, timeObs)]
  
  cvResultMeltGroup =
    cvResultMelt[, .(medae = median(abs(timeError))), by = .(sumabsv, nSpc)]
  ggplot(cvResultMeltGroup) +
    geom_point(aes(x = nSpc, y = medae, shape = sumabsv, color = sumabsv), size = 2) +
    labs(x = 'Number of SPCs', y = 'Median absolute error') +
    theme_bw() + theme(legend.position = c(0.7, 0.7))
  
  #get the best sumabsv and SPCs
  SPCs_optimal= cvResultMeltGroup %>% arrange(medae) %>% pull(nSpc) %>% .[1]
  sumabsv_optimal= cvResultMeltGroup %>% arrange(medae) %>% pull(sumabsv) %>% .[1] %>% as.character() %>% as.numeric()
  
  #Train a model on the full test dataset
  fitResultFinal = zeitzeigerFit(x, time)
  spcResultFinal = zeitzeigerSpc(
    fitResultFinal$xFitMean, fitResultFinal$xFitResid, sumabsv = 3) #but the optimal sumabsv_optimal is 3
  
  #test
  predResult = zeitzeigerPredict(x, time, xTest, spcResultFinal, nSpc = 4) #but the optimal SPCs_optimal is 4
  #alculate the difference between predicted time and observed time
  dfTest = data.frame(
    timeObs = timeTest, timePred = predResult$timePred,
    timeError = getCircDiff(predResult$timePred, timeTest))
  
  dfTest= 24*dfTest
  avg_error= (sum(abs(dfTest$timeError))/nrow(dfTest)) %>% round(.,digits = 2)
  error_list_less[[k]]= avg_error
  
  
  #Plot the coefficients of the features for the SPCs
  for(i in 1:length(spcResultList)){
    coefficient_list[[i]]= list()
    for(j in 1:length(spcResultList[[i]])){
      v = data.frame(spcResultList[[i]][[j]]$v[, 1:3])
      colnames(v) = c('SPC 1', 'SPC 2', 'SPC 3')
      rownames(v)= colnames(x)
      v = v[apply(v, 1, function(r) any(r != 0)), ]
      v[v == 0] = NA
      v = v[do.call(order, v), ]
      v$feature = rownames(v)
      vMelt = melt(setDT(v), id.vars = 'feature', variable.name = 'spc',
                   value.name = 'Coefficient')
      vMelt[, feature := factor(feature, rev(v$feature))]
      
      coefficient_list[[i]][[j]]= vMelt
    }
  }
  
  coefficient_list_final_less[[k]]= coefficient_list
}

error_final= unlist(error_list_more) #no difference
error_final= unlist(error_list_less) #no difference
save(coefficient_list_final_less,coefficient_list_final_more,file = paste0(output_dir,project,'_coefficient_list.RData'))

res= unlist(coefficient_list_final_less,recursive = F) %>% lapply(.,rbindlist)%>% rbindlist %>%
  dplyr::filter(!is.na(Coefficient)) 
res_shared_features= table(res$feature)%>% sort(decreasing = T)  #shared coefficient
genes_to_use_for_prediction_final_less= names(res_shared_features)[res_shared_features>100]
genes_to_use_for_prediction_final_less= c('Nr1d2','Ucp2','Cry1','Dbp','Arntl','Tef','Polq','Rorc','Arfgap3','Per3','Htr2b','Pdgfrl',
                                          'Lrig1','Per2','Fanci','Gja1','Fam76a','Polb','Rev1','Tppp3','Pafah1b3','Npas2','Lonp2','Scaper',
                                          'Mthfd1l','Rsad1','Guk1','Leo1','Has2','Nr1d1','Kdelr3', #herer >150
                                          'Lama3','Ift80','Col6a1','Prr13','Bscl2','Mrpl45','D330050G23Rik','Zfp493','Cstf1','Col5a2','Col1a2',
                                          'Nr1i3','Secisbp2','Rad54b','Wrap53','Ikbip','Pcolce','Steap4','Tubd1','Ceacam19','Hlf','Ino80','Ahsa1')

res= unlist(coefficient_list_final_more,recursive = F) %>% lapply(.,rbindlist)%>% rbindlist %>%
  dplyr::filter(!is.na(Coefficient)) 
res_shared_features= table(res$feature)%>% sort(decreasing = T)  #shared coefficient
genes_to_use_for_prediction_final_more= names(res_shared_features)[res_shared_features>100]
genes_to_use_for_prediction_final_more= c('Nr1d2','D330050G23Rik','Tomm6','Dbp','Tef','Per3','Prr13','Cstf1','Ttc21b','Arid4b','Arntl',
                                          'Rev1','Cry1','Gpsm3','Guk1','Ucp2','Rbm3','Sh2d1b1','Brca1','Snrpd3','Thop1','Enho','Orai3',
                                          'Vrk3','Polb','Tppp3','Dnajc30','Fyb','Fam76a','Rrp7a','Cacybp','Pi4k2b','Hnrnpa1','Enox1', #here >150
                                          'Per2','Ofd1','Incenp','Malat1','Pdcl3','Prrc1','Slc25a13','Fzr1','Htr2b','Rad23a','Zfp511',
                                          'Gnai2','Cul4b','Ndufaf2','Tial1')


