rm(list = ls())
library(Seurat)
library(data.table)
library(plyr)
library(tidyverse)
library(ggplot2)
library(dittoSeq)
input_dir <- "./input/"
if(!dir.exists(input_dir)){dir.create(input_dir,recursive = T)} 
output_dir= './output/mMel_in/revision/' #
if(!dir.exists(output_dir)){dir.create(output_dir,recursive = T)} 
marker_dir <- "./marker/mMel_in/"
if(!dir.exists(marker_dir)){dir.create(marker_dir,recursive = T)} 

project <- 'mMel_in'
source('Modified function.R')
load(paste(input_dir,project, ".RData", sep = ""))
color_use= c('#F19B9B','#F1C3C1','#BCBBDD','#A7C6E8')

#############Further subcluster CD8#############
seu= subset(mMel_in,subclustering_round1.5=='CD8')
seu= seu %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() 
ElbowPlot(seu, ndims = 50)
n.pcs <- 15
n.neib=20 
seu= seu%>% FindNeighbors(dims = 1:n.pcs,k.param = n.neib) %>% FindClusters(resolution = 1.0) %>% RunUMAP(dims = 1:n.pcs, n.neighbors = n.neib)
p1= DimPlot_modified(seu)
p1
ggsave(paste(output_dir,project,'_Dimplot_CD8T_round2.pdf', sep=""), 
       p1, width = 10, height = 10, limitsize = F)

markers=c('Tbx21', #ILC1 'CD69','IFNG','ITGA1','IL7low'
          'Gata3','Areg','Bcl11b',#ILC2
          "Ncam1","Fcgr4","B3gat1","Klrc1","Ncr1", #NK
          "Cd3d","Cd3e", #T
          "Trac","Trbc1","Trbc2",
          "Klrb1","Gzmh",'Cd1d1','Cdc42se1','Ckb',#NKT KLRB1=NK1.1 CD3E, FCGR3A/NCAM1
          #'RACK1','GAS5','SARAF','UQCRB', #NKT markers calculated 'TMSB4XP4','MTATP6P1',
          "Slc4a10","Ccr5","Ccr6","Ccl20","Rora", #MAIT, also express KLRB1 "CXCR6",
          "Cd8a","Cd8b1", #CD8 T 
          "Cd4", 
          "Cd40lg","Cxcr5","Tnfrsf4", #Tfh TNFRSF4=OX40. PDCD1, IL7R, B3GAT1=CD57
          "Sell","Ccr7","Lef1", #Naive/TCM SELL=CD62L
          "Gzma","Gzmb","Gnly","Prf1","Nkg7","Il2","Ifng","Tnf","Il4ra", #Effector IL7r-
          'Ikzf2', #terminal effector
          "Cd44","Gzmk",'Klrg1', #Effector memory IL7r+
          'Il2ra','Cd27', #Central memory IL7r+
          "Cxcr6","Itga1","Itgae","Cd69", #TRM  "Ccr10",'Cxcr3','Ccr4','Ccr8','Cla','Il15ra','Fabp4',
          "Lag3","Pdcd1","Layn",'Entpd1',"Havcr2","Eomes",#Exhausted HAVCR2=TIM3 'CD244',
          # "Layn","Eomes",#'CD244',
          'Xcl1','Tox','Tnfsf4','Cx3cr1','Tcf7', #pre-exhausted Tcf1,'Il7r' high. Gzma, Gzmb, Prf1 low
          "Foxp3","Ctla4", #Treg ,"IL10","KLRG1","FOLR4","Il2ra",
          "Mki67","Top2a" #Mitotic
)

p1= DotPlot(seu,features = markers, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90,hjust=0.5,vjust = 0.5), plot.background = element_rect(fill = "white")) 
p1
ggsave(paste(output_dir,project,'_Dotplot_CD8T_round2.pdf', sep=""), 
       p1, width = 15, height = 10, limitsize = F)


Idents(seu)=seu@meta.data$seurat_clusters
new.cluster.ids <- c("CD8_ex_term", #0 Klr*
                     "CD8_ex_int", "CD8_em", "CD8_ex_int", "CD8_em", "CD8_ex_prog", #
                     "CD8_ef","CD8_ef",'CD8_ex_prog')  
names(new.cluster.ids) <- levels(seu)
seu <- RenameIdents(seu, new.cluster.ids)
order= c('CD8_ef',"CD8_em","CD8_ex_prog","CD8_ex_int","CD8_ex_term")
Idents(seu)=factor(Idents(seu), levels = order)
seu <- AddMetaData(seu, metadata = seu@active.ident, col.name = 'subclustering_round3')
p1=DimPlot_modified(seu, label = T,label.box = F,prefix = 'CD8_')
ggsave(paste(output_dir,project,'_UMAP_T_annotated_round3.pdf', sep=""), 
       p1, width = 10, height = 10, limitsize = F)
p1


p1=FeaturePlot(seu,features = c('Gzmk','Ifng','Gzma',"Pdcd1",'Tox',"Havcr2","Entpd1","Lag3","Tcf7",'Slamf6','Ctla4','Tnfrsf9'),order = T)
p1
ggsave(paste(output_dir,project,'_CD8T_features.pdf', sep=""), 
       p1, width = 15, height = 10, limitsize = F)

p1=StackedVlnPlot_modified(seu,features = 'Pdcd1',split.by = 'timepoint',x_axis_text_size = 20,angle.x = 0,hjust.x = 0.5,vjust.x = 1,
                        pt.size = 0.1,color.use = color_use)& 
  theme(legend.position = 'none',text = element_text(size=20), axis.title.y = element_text(face='bold',angle=90, vjust = 0.5))&
  stat_summary(fun = mean, geom='point', size = 8, colour = "black", shape = 95,position = position_dodge(width = 0.9))
p1
ggsave(filename = paste0(output_dir,project,'_CD8T_PD1_ZT.pdf'),p1,width = 10,height = 4)

plotFractionPvalue(seu@meta.data,name = 'CD8T_sub_annotated',parameter = 'subclustering_round3',compare = F,paired = F,height = 6,
                   origin_levels = c('ZT01','ZT07','ZT13','ZT19'))


seu.annotated.markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
#VlnPlot of top marker genes
library(CellChat)
for(i in levels(seu$subclustering_round3)){
  genes= seu.annotated.markers %>% filter(cluster==i) %>% slice_max(avg_log2FC,n=20) %>% pull(gene)
  p1=StackedVlnPlot(seu,features = genes,color.use = dittoSeq::dittoColors())
  ggsave(filename = paste0(output_dir,project,'_',i,'_marker_gene_Vln.pdf'),p1,width = 10,height = 15)
}


#replot total T cell UMAP
load(file = paste(input_dir, project,"_T.RData", sep = ""))
T_mMel_in@meta.data$subclustering_round3= as.character(T_mMel_in@meta.data$subclustering_round2)
T_mMel_in@meta.data[match(seu@meta.data$Cell,T_mMel_in@meta.data$Cell),'subclustering_round3']= seu@meta.data$subclustering_round3 %>% as.character()
new_levels= c(levels(T_mMel_in$subclustering_round2)[1:3],levels(seu$subclustering_round3),levels(T_mMel_in$subclustering_round2)[7:13])
T_mMel_in@meta.data$subclustering_round3= factor(T_mMel_in@meta.data$subclustering_round3,levels = new_levels)
p1=DimPlot_modified(T_mMel_in, label = T,label.size = 6,prefix = 'T',legend_col=4,label.box=F,group.by = 'subclustering_round3')
p1
ggsave(paste(output_dir,project,'_UMAP_T_annotated_round2.pdf', sep=""), 
       p1, width = 10, height = 10, limitsize = F)
#fraction
plotFractionPvalue(T_mMel_in@meta.data,name='T',parameter = 'subclustering_round3',
                   origin_levels=c('ZT01','ZT07','ZT13','ZT19'),
                   paired = T,psize = 5,fontsize = 15,width = 23)

#Stacked fraction
library(dittoSeq)
p1=dittoBarPlot(T_mMel_in,"subclustering_round3", group.by = "Sample_origin",x.labels.rotate = 0,retain.factor.levels = T) + 
  theme(text= element_text(size=18), axis.text.x = element_text(size=18) ,plot.title = element_blank(),axis.title.x = element_blank())
p1
ggsave(paste(output_dir,project,'_round3_T_fraction_stacked.pdf', sep=""), p1, width = 10, height = 10)


############circadian clock genes no scatter: expm ed#############
mMel_in= subset(mMel_in,subset= subclustering_round1!= 'Melanocyte') #remove melanocyte
cells= mMel_in$subclustering_round1.5 %>% droplevels() %>% levels()
genes= c('Arntl','Per1','Per2','Per3','Clock','Cry1','Cry2','Nr1d1','Nr1d2','Bhlhe40','Bhlhe41','Dbp','Nfil3','Rora','Rorb') #
Idents(mMel_in)= mMel_in@meta.data$subclustering_round1.5

library(ggpubr)
source('Circadian_fit.R')
library(gggap)
exp_table= mMel_in@assays$RNA@data[genes,] %>% as.matrix()%>% t() %>% as.data.frame() %>% cbind(.,mMel_in@meta.data)  #%>% expm1() 
exp_table$timepoint= exp_table$timepoint %>% stringr::str_extract(.,"[[:digit:]]+") %>% as.integer() 

p1= list()
i='CD45'
p1[[i]]= list()
tmp= exp_table
for(j in genes){
  output= fit_cosinor_cosinor_no_jitter(data= tmp,variable = 'timepoint',observation = j,replicated_obs = F,group = NULL,color_use = color_use,
                                        two_axis = F,bar = F,stat_summary = F,errorbar = T,CI=F)
  p1[[i]][[j]]= output$p+ ggtitle(paste(i,j))
  #p1[[j]]=p1[[j]] %>% gggap(segments= c(0.1999999,0.2),ylim=c(0,4),tick_width = c(0.1,1))
}

for(i in cells){
  p1[[i]]= list()
  tmp= exp_table %>% dplyr::filter(subclustering_round1.5==i) 
  for(j in genes){
    output= fit_cosinor_cosinor_no_jitter(data= tmp,variable = 'timepoint',observation = j,replicated_obs = F,group = NULL,color_use = color_use,
                                          two_axis = F,bar = F,stat_summary = F,errorbar = T,CI=F)
    p1[[i]][[j]]= output$p+ ggtitle(paste(i,j))
    #p1[[j]]=p1[[j]] %>% gggap(segments= c(0.1999999,0.2),ylim=c(0,4),tick_width = c(0.1,1))
  }
}

combined_plot= cowplot::plot_grid(plotlist= list(p1$CD45$Per1),nrow= 1,align = "v")
ggsave(paste(output_dir,project,'_annotated_round1.5_CD45_Per1_errorbar.pdf', sep=""),combined_plot,width = 4.5, height = 4)


p2= list(p1[['Mo']][['Arntl']],p1[['Mo']][['Per1']],p1[['Mo']][['Per2']],p1[['Mo']][['Nr1d2']],p1[['Mo']][['Nr1d1']],p1[['Mo']][['Nfil3']],p1[['Mo']][['Cry2']],p1[['Mo']][['Cry1']],p1[['Mo']][['Bhlhe40']],
         p1[['Mac']][['Arntl']],p1[['Mac']][['Per1']],p1[['Mac']][['Per2']],p1[['Mac']][['Nr1d2']],p1[['Mac']][['Nr1d1']],p1[['Mac']][['Nfil3']],p1[['Mac']][['Cry2']],p1[['Mac']][['Cry1']],p1[['Mac']][['Bhlhe40']],
         p1[['cDC2']][['Arntl']],p1[['cDC2']][['Per1']],p1[['cDC2']][['Per2']],p1[['cDC2']][['Nr1d2']],p1[['cDC2']][['Nr1d1']],p1[['cDC2']][['Nfil3']],p1[['cDC2']][['Cry1']],p1[['cDC2']][['Cry2']],p1[['cDC2']][['Bhlhe40']],
         p1[['cDC1']][['Arntl']],p1[['cDC1']][['Per1']],p1[['cDC1']][['Per2']],p1[['cDC1']][['Nr1d2']],p1[['cDC1']][['Nr1d1']],p1[['cDC1']][['Nfil3']],p1[['cDC1']][['Cry2']],p1[['cDC1']][['Cry1']],p1[['cDC1']][['Bhlhe40']],
         p1[['CD8']][['Arntl']],p1[['CD8']][['Per1']],p1[['CD8']][['Per2']],p1[['CD8']][['Nr1d2']],p1[['CD8']][['Nr1d1']],p1[['CD8']][['Nfil3']],p1[['CD8']][['Cry2']],p1[['CD8']][['Cry1']],p1[['CD8']][['Bhlhe40']],
         p1[['CD4']][['Arntl']],p1[['CD4']][['Per1']],p1[['CD4']][['Per2']],p1[['CD4']][['Nr1d2']],p1[['CD4']][['Nr1d1']],p1[['CD4']][['Nfil3']],p1[['CD4']][['Cry1']],p1[['CD4']][['Cry2']],p1[['CD4']][['Bhlhe40']],
         p1[['Treg']][['Arntl']],p1[['Treg']][['Per1']],p1[['Treg']][['Per2']],p1[['Treg']][['Nr1d2']],p1[['Treg']][['Nr1d1']],p1[['Treg']][['Nfil3']],p1[['Treg']][['Cry2']],p1[['Treg']][['Cry1']],p1[['Treg']][['Bhlhe40']],
         p1[['NK']][['Arntl']],p1[['NK']][['Per2']],p1[['NK']][['Per1']],p1[['NK']][['Nr1d2']],p1[['NK']][['Nr1d1']],p1[['NK']][['Nfil3']],p1[['NK']][['Cry1']],p1[['NK']][['Cry2']],p1[['NK']][['Bhlhe41']])
combined_plot= cowplot::plot_grid(plotlist= p2,ncol=9,nrow= 8,align = "v")
ggsave(paste(output_dir,project,'_annotated_round1.5_circadian_genes_selected_pub_errorbar.pdf', sep=""),combined_plot,width = 40.5, height = 32)
p2= list(p1[['Mo']][['Arntl']],p1[['Mo']][['Per1']],p1[['Mo']][['Cry2']],p1[['Mo']][['Nfil3']],p1[['Mo']][['Bhlhe40']],
         p1[['Mac']][['Arntl']],p1[['Mac']][['Per1']],p1[['Mac']][['Cry2']],p1[['Mac']][['Nfil3']],p1[['Mac']][['Bhlhe40']],
         p1[['cDC2']][['Arntl']],p1[['cDC2']][['Per1']],p1[['cDC2']][['Cry1']],p1[['cDC2']][['Nfil3']],p1[['cDC2']][['Bhlhe40']],
         p1[['cDC1']][['Arntl']],p1[['cDC1']][['Per1']],p1[['cDC1']][['Cry2']],p1[['cDC1']][['Nfil3']],p1[['cDC1']][['Bhlhe40']],
         p1[['CD8']][['Arntl']],p1[['CD8']][['Per1']],p1[['CD8']][['Cry2']],p1[['CD8']][['Nfil3']],p1[['CD8']][['Bhlhe40']],
         p1[['CD4']][['Arntl']],p1[['CD4']][['Per1']],p1[['CD4']][['Cry1']],p1[['CD4']][['Nfil3']],p1[['CD4']][['Bhlhe40']]
         #p1[['Treg']][['Arntl']],p1[['Treg']][['Per1']],p1[['Treg']][['Cry2']],p1[['Treg']][['Nfil3']],p1[['Treg']][['Bhlhe40']],
         #p1[['NK']][['Arntl']],p1[['NK']][['Per2']],p1[['NK']][['Cry1']],p1[['NK']][['Nfil3']],p1[['NK']][['Bhlhe41']]
         )
combined_plot= cowplot::plot_grid(plotlist= p2,ncol=5,nrow= 6,align = "v")
ggsave(paste(output_dir,project,'_annotated_round1.5_circadian_genes_selected_pub_part_errorbar.pdf', sep=""),combined_plot,width = 22.5, height = 27)


p2= list(p1[['Mo']][['Per1']],p1[['Mac']][['Per1']],p1[['cDC2']][['Per1']],p1[['cDC1']][['Per1']],p1[['B']][['Per1']],p1[['CD8']][['Per1']],
         p1[['CD4']][['Per1']],
         p1[['Mo']][['Cry2']],p1[['Mac']][['Cry2']],p1[['cDC2']][['Cry2']],p1[['cDC1']][['Cry2']],p1[['B']][['Cry2']],p1[['CD8']][['Cry2']],
         p1[['CD4']][['Cry2']],
         p1[['Mo']][['Nr1d2']],p1[['Mac']][['Nr1d2']],p1[['cDC2']][['Nr1d2']],p1[['cDC1']][['Nr1d2']],p1[['B']][['Nr1d2']],p1[['CD8']][['Nr1d2']],
         p1[['CD4']][['Nr1d2']])
combined_plot= cowplot::plot_grid(plotlist= p2,ncol=7,nrow= 3,align = "v")
ggsave(paste(output_dir,project,'_annotated_round1.5_circadian_genes_Per1_Cry2_Nr1d2_errorbar.pdf', sep=""),combined_plot,width = 31.5, height = 13.5)


p2= list(p1[['Mo']][['Arntl']],p1[['Mac']][['Arntl']],p1[['cDC2']][['Arntl']],p1[['cDC1']][['Arntl']],p1[['B']][['Arntl']],p1[['CD8']][['Arntl']],
         p1[['CD4']][['Arntl']],
         p1[['Mo']][['Per2']],p1[['Mac']][['Per2']],p1[['cDC2']][['Per2']],p1[['cDC1']][['Per2']],p1[['B']][['Per2']],p1[['CD8']][['Per2']],
         p1[['CD4']][['Per2']],
         p1[['Mo']][['Nr1d1']],p1[['Mac']][['Nr1d1']],p1[['cDC2']][['Nr1d1']],p1[['cDC1']][['Nr1d1']],p1[['B']][['Nr1d1']],p1[['CD8']][['Nr1d1']],
         p1[['CD4']][['Nr1d1']])
combined_plot= cowplot::plot_grid(plotlist= p2,ncol=7,nrow= 3,align = "v")
ggsave(paste(output_dir,project,'_annotated_round1.5_circadian_genes_Per2_Arntl_Nr1d1_errorbar.pdf', sep=""),combined_plot,width = 31.5, height = 13.5)



############random sample to generate pseudocells############
seu= subset(mMel_in,subclustering_round1.5=='CD8')
other_cells= list(seu)
names(other_cells)= c('CD8')
genes_shared= rownames(seu)
other_cells_sample= list()
source('Circadian_fit.R')
genes= c('Pdcd1','Arntl','Per1','Per2','Per3','Clock','Cry1','Cry2','Nr1d1','Nr1d2','Bhlhe40','Bhlhe41','Dbp','Nfil3','Rora','Rorb') #

for(seed in 2009:2018){
  other_cells_sample[[seed]]=list()
  addTaskCallback(function(...) {set.seed(seed);TRUE}) #this will keep the seed 2023 used for the whole session! 
  for(j in names(other_cells)){
    cells_to_sample= sample(rownames(other_cells[[j]]@meta.data),round(nrow(other_cells[[j]]@meta.data)*0.2))
    other_cells_sample[[seed]]= other_cells[[j]] %>% subset(cells=cells_to_sample,features= genes_shared) %>% 
      AverageExpression(assays = 'RNA',return.seurat = F,group.by = 'mouseID') %>% 
      .$RNA %>% t()%>% as.data.frame() %>% rownames_to_column('mouseID')
    other_cells_sample[[seed]]$timepoint= str_split(other_cells_sample[[seed]]$mouseID,'_',simplify = T)[,1]
  }
}

exp_table= other_cells_sample %>% do.call(rbind,.)
exp_table$timepoint= exp_table$timepoint %>% stringr::str_extract(.,"[[:digit:]]+") %>% as.integer() 

p1= list()
i='CD8'
p1[[i]]= list()
tmp= exp_table
for(j in genes){
  #output= fit_cosinor_cosinor(data= tmp,jitter_width = 2,variable = 'timepoint',observation = j,replicated_obs = F,group = NULL,color_use = color_use,two_axis = F,stat_summary = T)
  output= fit_cosinor_cosinor_no_jitter(data= tmp,variable = 'timepoint',observation = j,replicated_obs = F,group = NULL,color_use = color_use,two_axis = F,
                                        bar = F,stat_summary = F,errorbar = T,CI = F)
  p1[[i]][[j]]= output$p+ ggtitle(paste(i,j))+ylim(0,NA)
  #p1[[j]]=p1[[j]] %>% gggap(segments= c(0.1999999,0.2),ylim=c(0,4),tick_width = c(0.1,1))
}
combined_plot= cowplot::plot_grid(plotlist= p1$CD8,nrow= 4,align = "v")
ggsave(paste(output_dir,project,'_annotated_round1.5_CD8_sampled_pseudo_circadian_errobar.pdf', sep=""),combined_plot,width = 16, height = 16)

ggplot(tmp,aes(timepoint, Pdcd1,colors=timepoint))+geom_jitter()


###########Vlnplot and percent data##########
#CD8 percent data
seu= subset(mMel_in,subclustering_round1.5=='CD8')
seu@meta.data$subclustering_round2= droplevels(seu@meta.data$subclustering_round2)
plotFractionPvalue(seu@meta.data,name = 'CD8T_round2',parameter = 'subclustering_round2',compare = F,paired = F,height = 4,width = 4,
                   origin_levels = c('ZT01','ZT07','ZT13','ZT19'))






