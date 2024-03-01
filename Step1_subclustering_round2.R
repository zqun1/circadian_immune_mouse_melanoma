rm(list = ls())
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
load(paste(input_dir,project, ".RData", sep = ""))

##########UMAP and dotplot######
mMel_in= subset(mMel_in,subset= subclustering_round2!='doublets'&subclustering_round2!='Low_quality') 
mMel_in #50053 features across 30603 samples within 2 assays

cat(paste(shQuote(unique(mMel_in$subclustering_round1.5))),sep=',')
order= c('B','cDC1','cDC2','cDC_activated','DC_mitotic','pDC','Mac','Mo','Neutrophil','NK','NKT','DETC','gdT','CD4','Treg','CD8','T','T_mitotic','Melanocyte')
mMel_in@meta.data$subclustering_round1.5= factor(mMel_in@meta.data$subclustering_round1.5,levels = order)
p1=DimPlot_modified(mMel_in,group.by='subclustering_round1.5',prefix='C',label.box=F,label.size = 6,legend_col=6)
p1
ggsave(paste(output_dir,project,'_UMAP_all_annotated_round1.5.pdf', sep=""), 
       p1, width = 10, height = 10, limitsize = F)

cat(paste(shQuote(unique(mMel_in$subclustering_round2))),sep=',')
order= c('B','cDC1','cDC1_mitotic','cDC2','cDC2_mitotic','cDC_activated','pDC',
         'Mac_MHC','Mac_C1qc','Mac_Ifn','Mac_Spp1','Mac_Mrc1','Mac_mitotic','cMo','ncMo','Neutrophil',
         'NK','NK_Ifn','NK_mitotic','NKT','DETC','gdT','CD4_n','CD4_ef','Treg','Treg_mitotic',
         'CD8_ef','CD8_em','CD8_ex','CD8_mitotic','T_Ifn','T_ribo','Melanocyte')
mMel_in@meta.data$subclustering_round2= factor(mMel_in@meta.data$subclustering_round2,levels = order)

markers= c("Cd79a",'G0s2','Siglech','Lyz2','H2-Aa',"Clec9a",'Ccr7','Cd68','Cd207',"Cd3d",'Cd8a','Foxp3','Trdc','Nkg7','Birc5','Mlana')
p1=FeaturePlot(mMel_in ,features = markers, label = F)
p1
ggsave(paste(output_dir,project,'_annotated_round1.5_Immune.png', sep=""),
       p1, width = 15, height = 15, limitsize = F)


p1=DimPlot_modified(mMel_in,group.by='subclustering_round2',prefix='C',label.box=F,label.size = 6,legend_col=6)
p1
ggsave(paste(output_dir,project,'_UMAP_all_annotated_round2.pdf', sep=""), 
       p1, width = 10, height = 10, limitsize = F)


markers = c("Cd79a", "Igkc", #B "IGLC3","IGHG3","IGHA2","CD79B",
            "H2-Aa","H2-Ab1",'Itgax',
            "Clec9a","Xcr1",
            "Irf4",'Cd209a',"Clec10a",
            "Ccr7","Fscn1", #LAMP3+
            "Siglech","Bst2", #pDC
            "Itgam",'Adgre1',"Cd68", #'Ctsb',
            'Spp1',"C1qc",'Isg15','Mrc1', 
            'Fcn1','Csf1r', # monocyte 'Adgre1',"Itgam",
            'Ly6c2','Sell','Ccr2',#Ly6C+/- monocyte Ccr2- Sell- 'Ly6c1',
            "G0s2","Csf3r",
            "Nkg7","Klrd1", #NK "GNLY",'KLRF1',
            'Klrb1c','Klra1',#NKT KLRB1=NK1.1 'Cd7','Klra7',
            "Cd3d",'Cd3e',
            'Trdv4', #DETC 'Col27a1',
            "Tcrg-C1","Trdc", #gdT cell 
            "Trac","Trbc1", #T "CD3G","TRBC2","CD3E","TRBC1",
            'Cd4',
            "Foxp3","Ctla4",
            "Cd8a","Cd8b1", #CD8 T 
            "Birc5","Top2a", #Mitotic "MKI67","CKS1B",
            'Mlana', 'Dct') #Melanocyte
p1=DotPlot(mMel_in, features = markers,group.by = 'subclustering_round2', cols = "RdBu") + theme(axis.text.x = element_text(angle = 90,hjust=0.5,vjust = 0.5),plot.background = element_rect(fill = "white")) 
ggsave(paste(output_dir,project,'_Dotplot_marker_annotated_round2.pdf', sep=""),
       p1, width = 15, height = 10, limitsize = F)
p1
p1=DotPlot(mMel_in, features = markers,group.by = 'subclustering_round1.5', cols = "RdBu") + theme(axis.text.x = element_text(angle = 90,hjust=0.5,vjust = 0.5),plot.background = element_rect(fill = "white")) 
ggsave(paste(output_dir,project,'_Dotplot_marker_annotated_round1.5.pdf', sep=""),
       p1, width = 15, height = 10, limitsize = F)
p1

save(mMel_in, mMel_in.markers, mMel_in.annotated.markers,file = paste(input_dir, project,".RData", sep = ""))

#######plot heatmap according to subclustering_round1.5#######
plan("multiprocess", workers = 4)
options(future.globals.maxSize= 50*1024^3)

Idents(mMel_in)= mMel_in@meta.data$subclustering_round1.5
mMel_in.annotated1.5.markers <- FindAllMarkers(mMel_in, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
write.csv(mMel_in.annotated1.5.markers,file= paste0(marker_dir,project,"_all_annotated_1.5_markers.csv"))
#mMel_in.annotated1.5.markers= read.csv(row.names=1,file= paste0(marker_dir,project,"_all_annotated_1.5_markers.csv"))

features<- mMel_in.annotated1.5.markers %>% group_by(cluster) %>% slice_max(n=5,order_by =avg_log2FC) 
p1= plotmarkers_heatmap(mMel_in,features,group.by = 'subclustering_round1.5',annotation_colors = 'dittoseq', #default is hue_pal colors
                        annotation_legend = F)
pdf(file = paste(output_dir,project,'_top5_all_round1.5_ave_heatmap.pdf', sep=""), width = 10, height = 15)
p1
dev.off()


genes_to_label= c('Clec9a',
                     'Cd209a','Tnfsf9',
                     'Ccr7','Il12b','Cd83',
                     'Siglech',
                     'C1qc','Spp1','Lgals3','Il1b',
                     'Csf1r',
                     'Csf3r',
                     'Prf1','Fasl',
                     'Cd4','Icos',
                     'Foxp3',
                     'Ifng','Pdcd1','Lag3','Tox','Tnfrsf9',
                     'Trdc')
p1= plotmarkers_heatmap(mMel_in,mMel_in.annotated1.5.markers,group.by = 'subclustering_round1.5',annotation_colors = 'dittoseq', #default is hue_pal colors
                        label_genes=genes_to_label)
pdf(file = paste(output_dir,project,'_all_round1.5_ave_heatmap.pdf', sep=""), width = 10, height = 15)
p1
dev.off()


