rm(list = ls())
library(Seurat)
library(data.table)
library(plyr)
library(tidyverse)
library(ggplot2)
library(dittoSeq)
set.seed(2023)
input_dir <- "./input/"
if(!dir.exists(input_dir)){dir.create(input_dir,recursive = T)} 
output_dir <- "./output/mMel_in/DiscoRhythm_singlecell/"
if(!dir.exists(output_dir)){dir.create(output_dir,recursive = T)} 
marker_dir <- "./marker/mMel_in/"
if(!dir.exists(marker_dir)){dir.create(marker_dir,recursive = T)} 

project <- 'mMel_in'
source('Modified function.R')
load(paste(input_dir,project, ".RData", sep = ""))

#########single cell: subclustering_1.5######
output_dir <- "./output/mMel_in/DiscoRhythm_singlecell/round1.5/"
if(!dir.exists(output_dir)){dir.create(output_dir,recursive = T)} 
library(DiscoRhythm)
library(SummarizedExperiment)
library(openxlsx)
library(rmarkdown)
cells=levels(mMel_in$subclustering_round1.5) 

output= list()
for(i in cells){
  tmp= subset(mMel_in,subset= subclustering_round1.5==i) %>% GetAssayData(assay = "RNA", slot = "data")
  tmp= tmp %>% as.data.frame() %>% rownames_to_column(var='ID') #input requires an ID column at first
  se <- discoDFtoSE(tmp) #make summarizedExperiment object
  selectDataSE <- discoCheckInput(se) #row-wise checks for missing values and constant values
  #Since a lot of zeros and no technical replicates, it is not reasonable to remove outliers or combine replicates
  FinalSE= selectDataSE
  colData(FinalSE)= colData(FinalSE)[,-1]
  discoODAres <- discoODAs(FinalSE,period=24,method=c("CS",'JTK','LS'),
                       ncores=8,circular_t=FALSE)

  genes= lapply(discoODAres,rownames_to_column,var='gene') %>% lapply(.,filter,pvalue<0.05) %>% lapply(.,pull,gene)
  tmp3= list(cells=i,CS_JTK= intersect(genes$CS,genes$JTK),CS_LS=intersect(genes$CS,genes$LS),JTK_LS=intersect(genes$LS,genes$JTK))
  discoODAres$intersection= plyr::ldply(tmp3, rbind) %>% t() %>% as.data.frame() #very useful

  #write.xlsx(discoODAres,file = paste0(output_dir,project,'_DiscoRhythm_result_',i,'.xlsx'),rowNames=T)
  output[[i]]= discoODAres
}
save(output,file = paste0(output_dir,project,'_DiscoRhythm_output.RData'))
#load(file = paste0(output_dir,project,'_DiscoRhythm_output.RData'))

genes_count= lapply(output,FUN = function(x){lapply(x[1:2],rownames_to_column,var='gene') %>% 
    lapply(.,filter,pvalue<0.05) %>% lapply(.,pull,gene) %>% unlist %>% unique %>% length}) %>% 
  do.call(rbind,.) %>% as.data.frame()
genes_count$Cell= rownames(genes_count)
library(ggpubr)
p1= ggdotchart(genes_count,
  x = "Cell", y = "V1",
  add = "segments", sorting = "descending",
  ylab = "#Circadian genes", xlab='',title = "CS/JTK"
)+theme(text = element_text(size=24))
p1
ggsave(filename = paste0(output_dir,project,'_#circadian_genes_union_round1.5.pdf'),p1,width = 10,height = 10)

genes_count= lapply(output,FUN = function(x){x[[4]][,2] %>% .[.!='CS_JTK' &!is.na(.)] %>% length}) %>% 
  do.call(rbind,.) %>% as.data.frame()
genes_count$Cell= rownames(genes_count)
library(ggpubr)
p1= ggdotchart(genes_count,
               x = "Cell", y = "V1",
               add = "segments", sorting = "descending",
               ylab = "#Circadian genes", xlab='',title = "CS&JTK",size = 6
)+theme(text = element_text(size=36,face = 'bold'))
p1
ggsave(filename = paste0(output_dir,project,'_#circadian_genes_intersection_round1.5.pdf'),p1,width = 10,height = 10)

library(ComplexHeatmap)
library(circlize)
#therefore, individual heatmap and then show side by side
genes= lapply(output,FUN = function(x){lapply(x[1:2],rownames_to_column,var='gene') %>% 
    lapply(.,dplyr::filter,pvalue<0.05) %>% lapply(.,pull,gene) %>% unlist %>% unique}) 
p1= list()
genes_all= genes%>% do.call(c,.) %>% unique()
Cell= levels(mMel_in$subclustering_round1.5)
mMel_in@meta.data$CellType= interaction(mMel_in@meta.data$subclustering_round1.5,mMel_in@meta.data$timepoint,sep = '-',lex.order = TRUE) #interaction keep levels
exp_data= AverageExpression(mMel_in,features = genes_all,
                            assays = 'RNA',group.by = 'CellType') 
exp_data= exp_data$RNA
circadian_gene= c('Arntl','Per1','Clock','Cry1','Cry2','Nr1d2','Bhlhe40','Dbp','Nfil3','Rora','Per2','Per3','Nr1d1','Bhlhe41','Rorb')
#no common genes. used threshold 0.3
genes_highlight= genes%>% do.call(c,.) %>% table %>% sort(decreasing = T) %>% .[.> 0.9*length(Cell)] %>% names()
label_genes= c(circadian_gene) #genes_highlight, 
for(i in 1:length(Cell)){
  tmp= genes[[i]]
  data= exp_data[tmp,(4*i-3):(4*i)]
  colnames(data)= c('ZT01','ZT07','ZT13','ZT19')
  data= t(scale(t(data))) 
  top_annotation= columnAnnotation(Cell=anno_block(gp=gpar(fill=dittoColors()[i]),
                                                   labels =Cell[i], labels_gp = gpar(fontsize=12)))
  right_label_annotation= rowAnnotation(link = anno_mark(at = which(rownames(data)%in%label_genes), #custome label genes
                                                         labels = rownames(data)[which(rownames(data)%in%label_genes)],
                                                         labels_gp = gpar(fontsize=10)))
  tmp= Heatmap(data,name='Avg_exp',show_column_names=T,
              cluster_rows = T,cluster_columns = F,show_row_names = F,na_col = 'grey',
              column_names_gp = grid::gpar(fontsize = 12),
              show_row_dend = F,show_heatmap_legend = F,
              top_annotation = top_annotation)+right_label_annotation
  p1[[i]]= grid.grabExpr(draw(tmp)) 
}
p1= p1[genes_count$V1 %>% order(decreasing = T)] #arrange according to the number of genes oscillating
p1[[length(Cell)+1]]= grid.grabExpr(draw(Legend(at=c(-2,-1,0,1,2),title = 'Avg_exp',col_fun = circlize::colorRamp2(c(-2, 0, 1), c("blue", "white", "red")))))
p2= cowplot::plot_grid(plotlist = p1,nrow = 2)
ggsave(filename = paste0(output_dir,project,'_circadian_genes_heatmap_round1.5.pdf'),p2,width = 20,height = 12)


#ORA: only shared genes
library(clusterProfiler)
library(GSEABase)
library(msigdbr)
library(stringr)
geneset=msigdbr(species = "mouse")
geneset.use= geneset[geneset$gs_cat=="H",] #Hallmark, C2=curated gene sets, C5=Ontology.
geneset.use= geneset.use %>% dplyr::select(gs_name, gene_symbol,entrez_gene)
for(i in names(output)){
  whatever=5
  tmp= data.frame(gene=output[[i]][[4]][,2],avg_log2FC= whatever)
  tmp= tmp[-1,]
  n_gene= nrow(tmp)
  if(n_gene>=15){ #only do ORA if >=15 genes are significant
    PlotORA_avg(tmp,cell= i,avg_log2FC = 1,species = 'mouse',only.pos = T,GO_only = T)
  }
}



