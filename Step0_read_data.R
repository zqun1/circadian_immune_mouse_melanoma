rm(list = ls())
setwd('/media/zq/DATA/Zeng_Qun_UNIGE/Collaboration_Chen_melanoma')
library(Seurat)
library(data.table)
library(plyr)
library(tidyverse)
library(ggplot2)
library(dittoSeq)
input_dir <- "./input/"
if(!dir.exists(input_dir)){dir.create(input_dir,recursive = T)} 
output_dir <- "./output/mMel_in/"
if(!dir.exists(output_dir)){dir.create(output_dir,recursive = T)} 
marker_dir <- "./marker/mMel_in/"
if(!dir.exists(marker_dir)){dir.create(marker_dir,recursive = T)} 

project <- 'mMel_in'
source('Modified function.R')

#######read in data##########
mMel_in= readRDS(paste0(input_dir,'seu_obj_all_merged_filtered_processed_pop_def_man_1_alt_obj.rds'))
mMel_in #50053 features across 33601 samples within 2 assays
colnames(mMel_in@meta.data)
mMel_in@meta.data$Cell= rownames(mMel_in@meta.data)
FeatureScatter(mMel_in,feature1 = 'nFeature_RNA',feature2 = 'percent_mt') #already filtered.
#########dowstream#######
min.cells=20 
num.cells <- Matrix::rowSums(x = mMel_in@assays$RNA@counts > 0)

Mel_in <- NormalizeData(mMel_in)
# Score cells for cell cycle
ccgene_mouse= readxl::read_xlsx('mouse cell cycle genes.xlsx')
mMel_in <- CellCycleScoring(mMel_in, s.features = ccgene_mouse[ccgene_mouse$category=='s',]$symbol, g2m.features =ccgene_mouse[ccgene_mouse$category=='g2m',]$symbol)
# Turn mitoRatio into categorical factor vector based on quartile values
fn_mt= fivenum(mMel_in@meta.data$percent_mt)
mMel_in@meta.data$mitoFr <- cut(mMel_in@meta.data$percent_mt, 
                               breaks=c(-Inf, fn_mt[2], fn_mt[3], fn_mt[4], Inf), #adjust accordingly
                               labels=c("Low","Medium","Medium high", "High"))
# regress out nUMI and G2M/S differences.
mMel_in$CC.Difference <- mMel_in$S.Score - mMel_in$G2M.Score
mMel_in <- ScaleData(mMel_in, vars.to.regress = c("percent_mt","CC.Difference"))
mMel_in <- RunPCA(mMel_in, npcs = 100, verbose = FALSE) 
#basic visualization & clustering
#determine the appropriate number of PCs to use for further analysis
ElbowPlot(mMel_in, ndims = 50)
DimHeatmap(mMel_in, dims = c(1:12), cells = 500, balanced = TRUE)
DimHeatmap(mMel_in, dims = c(13:24), cells = 500, balanced = TRUE)
#for the EC atlas paper, we picked 8 dimensions
n.pcs <- 20
n.neib= 30
mMel_in <- FindNeighbors(mMel_in, reduction = "pca", dims = 1:n.pcs, k.param = n.neib) #k.param may need to be adjusted

res.used <- seq(0.1,1.5,by=0.2)
# Loop over and perform clustering of different resolutions
for(i in res.used){
  mMel_in <- FindClusters(object = mMel_in, resolution = i)
}
mMel_in <- FindClusters(mMel_in, resolution = 1) #resolution may need to be adjusted - lower resolution = less clusters
# Make plot
library(clustree)
clus.tree.out <- clustree(mMel_in) +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Set1") +
  scale_edge_color_continuous(low = "grey80", high = "red")
# Plot
clus.tree.out  #15,24,33,41,42,(46),46,48
ggsave(paste(output_dir,project,'_clustree.png', sep=""), 
       clus.tree.out, width = 10, height = 10, limitsize = F)
#Use resolution=1
mMel_in <- RunUMAP(mMel_in, dims = 1:n.pcs, min.dist = 0.3, n.neighbors = n.neib) #n.neighbors may need to be adjusted
mMel_in <- FindClusters(mMel_in, resolution = 0.5) #resolution may need to be adjusted - lower resolution = less clusters
p1=DimPlot_modified(mMel_in, label = T,prefix='C',legend_col=6,label.box=F,label.size=6)
p1
ggsave(paste(output_dir,project,'_UMAP.png', sep=""),
       p1, width = 10, height = 10, limitsize = F)
library(dittoSeq)
p1= dittoDimPlot(mMel_in,'seurat_clusters', split.by = "mouseID")+theme(strip.text.x = element_text(size = 12)) 
p1
ggsave(paste(output_dir,project,'_mouseID_UMAP1.png', sep=""), 
       p1, width = 10, height = 10, limitsize = F)


p1=DimPlot_modified(mMel_in, group.by = "mouseID",prefix= 'M',legend_col=4,label.box=F,label=F)
p1
ggsave(paste(output_dir,project,'_mouseID_UMAP.png', sep=""), 
       p1, width = 30, height = 10, limitsize = F)


#Very important. Set the default assay back to RNA for DE analysis and value plotting.
#Find markers
#find differentially expressed markers for every cluster compared to all remaining cells, and report only the positive ones
Idents(mMel_in)=mMel_in$seurat_clusters
library(future)
plan("multiprocess", workers = 4)
options(future.globals.maxSize= 50*1024^3)
DefaultAssay(mMel_in)<- "RNA"
mMel_in.markers <- FindAllMarkers(mMel_in, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.5,test.use = 'roc')
tmp=mMel_in.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
#save marker genes as .csv file
write.csv(mMel_in.markers, paste(marker_dir,project,".markers.csv", sep = ""), row.names = FALSE)
#mMel_in.markers= read.csv(paste(marker_dir,"mMel_in.markers.csv", sep = ""),header = T)

mMel_in.markers %>% filter(myAUC > 0.6) %>% count(cluster, name = 'number') #all clusters have >0 markers of AUC>0.6. Increasing the resolution is possible.


#annotation
Idents(mMel_in)= mMel_in$seurat_clusters
markers = c('Hba-a1','Hbb-b1', #Erythrocytes 
            "Cldn5", "Cdh5", "Pecam1","Vwf",  #EC "RAMP2","FLT1", 
            "Epcam","Krt19", "Cdh1", "Krt18",  
            'Mlana', 'Mitf','Dct', #Melanocyte
            "Col1a1", "Bgn","Col1a2",'S100a4', #Fibroblast ACTA2, TAGLN, MHY11, PDGFB "DCN","C1R","THY1",'FAP',
            "Rgs5","Cspg4","Abcc9","Kcnj8",'Pdgfrb', #Pericytes "TAGLN","ACTA2","MCAM","STEAP4", Ito cell/Hepatic stellate cell
            "Tagln","Acta2","Myh11","Synpo2", #SMC "CNN1","DES",MYLK,MYL9
            "S100b", "Plp1",'Gfap', #Glia 'AIF1': IBA1, 
            "Cd79a", "Igkc","Ighm", 'Mzb1', #B "IGLC3","IGHG3","IGHA2","CD79B",
            "Ncam1", "Nkg7","Klrd1", #NK "GNLY",'KLRF1',
            "Kit","Ms4a2", #Mast cell "GATA2", "TPSAB1", "CPA3",
            "G0s2","Csf3r","Cd33", #Neutrophils "FCGR3B",
            "Siglech","Gzmb","Bst2",'Irf8', #pDC 
            "H2-Aa","H2-Ab1","Itgax", "Clec9a",  #DC "LILRA4", "CCL17",
            "Cd207",'Cd209a',  #Langerhans cell
            "Lyz1","Marco","Cd68","Fcgr4", "Aif1",#Myeloid
            'Cd163','Vsig4','Mafb', #Kuffer
            "Cd3d","Trac","Trdc", #T "CD3G","TRBC2","CD3E","TRBC1",
            "Ptprc",
            "Birc5","Top2a") #Mitotic "MKI67","CKS1B",
p1=DotPlot(mMel_in, features = markers, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90,hjust=0.5,vjust = 0.5),plot.background = element_rect(fill = "white")) 
ggsave(paste(output_dir,project,'_Dotplot_marker.png', sep=""), 
       p1, width = 10, height = 10, limitsize = F)
p1
mMel_in.markers[mMel_in.markers$cluster==40,]%>%slice_max(n=50,order_by = avg_log2FC)

new.cluster.ids <- c("Myeloid",
                     "NK", "T", "B", "Myeloid", "Myeloid", #
                     "cDC2", "T_mitotic", "T","T","Myeloid",  #8 Treg. 9 CD8T
                     "Mac", "Neutrophil","T","cDC1","pDC", #
                     "B", "Melanocyte") #
names(new.cluster.ids) <- levels(mMel_in)
mMel_in <- RenameIdents(mMel_in, new.cluster.ids)
order= c("B","Myeloid",'Mac','cDC1','cDC2','pDC','Neutrophil',"T",'T_mitotic','NK','Melanocyte')
table(mMel_in@active.ident)
Idents(mMel_in)=factor(Idents(mMel_in), levels = order)
table(mMel_in@active.ident)
mMel_in <- AddMetaData(mMel_in, metadata = mMel_in@active.ident, col.name = 'subclustering_round1')
p1=DimPlot_modified(mMel_in, prefix = 'C',group.by = 'subclustering_round1',label.size=8,legend_col = 5,label.box = F)
ggsave(paste(output_dir,project,'_annotated_round1.png', sep=""),
       p1, width = 10, height = 10, limitsize = F)
p1
markers= c('Lyz2',"Cd79a",'G0s2','Siglech',"Clec9a",'Itgax',"Cd3d",'Nkg7','Birc5')
p1=FeaturePlot(mMel_in, features = markers, label = F)
p1
ggsave(paste(output_dir,project,'_annotated_round1_Immune.png', sep=""),
       p1, width = 10, height = 10, limitsize = F)



#Find markers
#find differentially expressed markers for every cluster compared to all remaining cells, and report only the positive ones
mMel_in.annotated.markers <- FindAllMarkers(mMel_in, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
mMel_in.annotated.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
#save marker genes as .csv file
write.csv(mMel_in.annotated.markers, paste(marker_dir,"mMel_in.markers.csv", sep = ""), row.names = FALSE)
#analyse the expression of marker genes in different ways:

save(mMel_in, mMel_in.markers, mMel_in.annotated.markers,file = paste(input_dir, project,".RData", sep = ""))
#load(paste(input_dir,project, ".RData", sep = ""))




#######plot Cell distribution between different time points########

mMel_in@meta.data$timepoint= factor(mMel_in@meta.data$timepoint,levels= c("ZT01",'ZT07',"ZT13",'ZT19'))
mMel_in@meta.data$Sample_origin=mMel_in@meta.data$timepoint
mMel_in@meta.data$patientID=mMel_in@meta.data$mouseID
plotFractionPvalue(metadata = mMel_in@meta.data,name = 'total',parameter = 'subclustering_round1',
                   origin_levels = levels(mMel_in@meta.data$Sample_origin),paired = T,width = 20)

library(dittoSeq)
p1= dittoBarPlot(mMel_in, 'subclustering_round1',group.by = 'timepoint',x.labels.rotate = 0,retain.factor.levels = T)
p1
ggsave(paste(output_dir,project,'_cell_fraction_stacked.png', sep=""), 
       p1, width = 10, height = 10, limitsize = F)
