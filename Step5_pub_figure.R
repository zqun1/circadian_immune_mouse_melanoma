rm(list = ls())
library(Seurat)
library(data.table)
library(plyr)
library(tidyverse)
library(ggplot2)
library(dittoSeq)
input_dir <- "./input/"
if(!dir.exists(input_dir)){dir.create(input_dir,recursive = T)} 
output_dir= './output/mMel_in/pub/' #
if(!dir.exists(output_dir)){dir.create(output_dir,recursive = T)} 
marker_dir <- "./marker/mMel_in/"
if(!dir.exists(marker_dir)){dir.create(marker_dir,recursive = T)} 

project <- 'mMel_in'
source('Modified function.R')
load(paste(input_dir,project, ".RData", sep = ""))
color_use= c('#F19B9B','#F1C3C1','#BCBBDD','#A7C6E8')

#######circadian genes in jitter plot round1.5#####
mMel_in= subset(mMel_in,subset= subclustering_round1!= 'Melanocyte') #remove melanocyte
cells= mMel_in$subclustering_round1.5 %>% droplevels() %>% levels()
genes= c('Arntl','Per1','Per2','Per3','Clock','Cry1','Cry2','Nr1d1','Nr1d2','Bhlhe40','Bhlhe41','Dbp','Nfil3','Rora','Rorb') #
Idents(mMel_in)= mMel_in@meta.data$subclustering_round1.5

library(ggpubr)
source('Circadian_fit.R')
library(gggap)
exp_table= mMel_in@assays$RNA@data[genes,] %>% as.matrix() %>% t() %>% as.data.frame() %>% cbind(.,mMel_in@meta.data) 
exp_table$timepoint= exp_table$timepoint %>% stringr::str_extract(.,"[[:digit:]]+") %>% as.integer() 

p1= list()
i='CD45'
p1[[i]]= list()
tmp= exp_table
for(j in genes){
  output= fit_cosinor_cosinor(data= tmp,variable = 'timepoint',observation = j,replicated_obs = F,group = NULL,color_use = color_use,jitter_width = 0.4,two_axis = T)
  p1[[i]][[j]]= output$p+ ggtitle(j)
  #p1[[j]]=p1[[j]] %>% gggap(segments= c(0.1999999,0.2),ylim=c(0,4),tick_width = c(0.1,1))
}
p1[[i]][['cells']]= ggplot()+annotate('text',label=i,x = 4, y = 25)+theme_void()

for(i in cells){
  p1[[i]]= list()
  tmp= exp_table %>% dplyr::filter(subclustering_round1.5==i) 
  for(j in genes){
    output= fit_cosinor_cosinor(data= tmp,variable = 'timepoint',observation = j,replicated_obs = F,group = NULL,color_use = color_use,jitter_width = 0.4,two_axis = T)
    p1[[i]][[j]]= output$p+ ggtitle(j)
    #p1[[j]]=p1[[j]] %>% gggap(segments= c(0.1999999,0.2),ylim=c(0,4),tick_width = c(0.1,1))
  }
  p1[[i]][['cells']]= ggplot()+annotate('text',label=i,x = 4, y = 25)+theme_void()
}

p2= list(p1[['Mo']][['cells']],p1[['Mo']][['Arntl']],p1[['Mo']][['Per1']],p1[['Mo']][['Nr1d2']],p1[['Mo']][['Nfil3']],p1[['Mo']][['Cry2']],p1[['Mo']][['Bhlhe40']],
         p1[['Mac']][['cells']],p1[['Mac']][['Arntl']],p1[['Mac']][['Per1']],p1[['Mac']][['Nr1d2']],p1[['Mac']][['Nfil3']],p1[['Mac']][['Cry2']],p1[['Mac']][['Bhlhe40']],
         p1[['cDC2']][['cells']],p1[['cDC2']][['Arntl']],p1[['cDC2']][['Per1']],p1[['cDC2']][['Nr1d2']],p1[['cDC2']][['Nfil3']],p1[['cDC2']][['Cry1']],p1[['cDC2']][['Bhlhe40']],
         p1[['cDC1']][['cells']],p1[['cDC1']][['Arntl']],p1[['cDC1']][['Per1']],p1[['cDC1']][['Nr1d2']],p1[['cDC1']][['Nfil3']],p1[['cDC1']][['Cry2']],p1[['cDC1']][['Bhlhe40']],
         p1[['CD8']][['cells']],p1[['CD8']][['Arntl']],p1[['CD8']][['Per1']],p1[['CD8']][['Nr1d2']],p1[['CD8']][['Nfil3']],p1[['CD8']][['Cry2']],p1[['CD8']][['Bhlhe40']],
         p1[['CD4']][['cells']],p1[['CD4']][['Arntl']],p1[['CD4']][['Per1']],p1[['CD4']][['Nr1d2']],p1[['CD4']][['Nfil3']],p1[['CD4']][['Cry1']],p1[['CD4']][['Bhlhe40']],
         p1[['Treg']][['cells']],p1[['Treg']][['Arntl']],p1[['Treg']][['Per1']],p1[['Treg']][['Nr1d2']],p1[['Treg']][['Nfil3']],p1[['Treg']][['Cry2']],p1[['Treg']][['Bhlhe40']],
         p1[['NK']][['cells']],p1[['NK']][['Arntl']],p1[['NK']][['Per2']],p1[['NK']][['Nr1d2']],p1[['NK']][['Nfil3']],p1[['NK']][['Cry2']],p1[['NK']][['Bhlhe41']])
combined_plot= cowplot::plot_grid(plotlist= p2,ncol=7,nrow= 8,align = "v")
ggsave(paste(output_dir,project,'_annotated_round1.5_circadian_genes_selected_pub.pdf', sep=""),combined_plot,width = 21, height = 20)


p2= list(p1[['Mo']][['cells']],p1[['cDC2']][['cells']],p1[['cDC1']][['cells']],p1[['B']][['cells']],p1[['CD8']][['cells']],p1[['CD4']][['cells']],p1[['NK']][['cells']],p1[['Mac']][['cells']],p1[['Treg']][['cells']],
         p1[['Mo']][['Per1']],p1[['cDC2']][['Per1']],p1[['cDC1']][['Per1']],p1[['B']][['Per1']],p1[['CD8']][['Per1']],p1[['CD4']][['Per1']],p1[['NK']][['Per2']],p1[['Mac']][['Per1']],p1[['Treg']][['Per1']],
         p1[['Mo']][['Cry2']],p1[['cDC2']][['Cry1']],p1[['cDC1']][['Cry2']],p1[['B']][['Cry2']],p1[['CD8']][['Cry2']],p1[['CD4']][['Cry2']],p1[['NK']][['Cry2']],p1[['Mac']][['Cry2']],p1[['Treg']][['Cry2']])
combined_plot= cowplot::plot_grid(plotlist= p2,ncol=9,nrow= 3,align = "v")
ggsave(paste(output_dir,project,'_annotated_round1.5_circadian_genes_selected_pub2.pdf', sep=""),combined_plot,width = 27, height = 7.5)

###############single-cell: upset plot of enriched GO of oscillatory genes###########
output_dir <- "./output/mMel_in/DiscoRhythm_singlecell/round1.5/"
load(file = paste0(output_dir,project,'_DiscoRhythm_output.RData'))
output_dir= './output/mMel_in/pub/' #

#ORA: only shared genes
library(clusterProfiler)
library(GSEABase)
library(msigdbr)
library(ComplexHeatmap)
p1= list()
for(i in names(output)){
  tmp= output[[i]][[4]][,2]
  tmp= tmp[-1]
  n_gene= length(tmp)
  if(n_gene>=15){ #only do ORA if >=15 genes are significant
    EC_out.GO= enrichGO(tmp,OrgDb = 'org.Mm.eg.db',keyType = 'SYMBOL',ont = 'BP')
    EC_out.GO.sim= clusterProfiler::simplify(EC_out.GO,cutoff=0.7)  
    if(nrow(fortify(EC_out.GO.sim))){ ##May have result but no enrichedResult. fortify to convert the enrichedResult to a dataframe. 
      p1[[i]] =clusterProfiler::dotplot(EC_out.GO.sim,showCategory=20)+
        scale_y_discrete(labels = function(x) str_wrap(x, width = 81)) #forcce break line if y label length is >81
    }
  }
}
save(p1,file = paste0(input_dir,project,'_DiscoRhythm_single_cell_round1.5_GO.RData'))

data= lapply(p1,function(x){x$data$Description %>% as.character()})
tmp= do.call(c,data)
tmp2= tmp %>% str_split(pattern = ' ') %>% unlist() %>% table() %>% sort(decreasing = T) 
grep('activation',tmp,value = T)
data=lapply(data,gsub,pattern='.*(cycle|mitotic|proliferation|DNA replication).*',replacement='Proliferation',ignore.case = T)
data=lapply(data,gsub,pattern='.*(T cell activation).*',replacement='Regulation of T cell activation',ignore.case = T)
data=lapply(data,gsub,pattern='.*(adhesion).*',replacement='Leukocyte adhesion',ignore.case = T)
data=lapply(data,gsub,pattern='.*(migration).*',replacement='Migration',ignore.case = T)
data=lapply(data,gsub,pattern='.*(apopto).*',replacement='Regulation of apoptosis',ignore.case = T)
data=lapply(data,gsub,pattern='.*(macromolecule biosynthetic).*',replacement='Macromolecule biosynthetic',ignore.case = T)
data=lapply(data,gsub,pattern='.*(mRNA).*',replacement='mRNA',ignore.case = T)
data=lapply(data,gsub,pattern='.*(amide).*',replacement='Amide metabolism',ignore.case = T)
data=lapply(data,gsub,pattern='.*(translation).*',replacement='Translation',ignore.case = T)
data=lapply(data,gsub,pattern='.*(ribonucleoprotein).*',replacement='Ribonucleoprotein complex',ignore.case = T)
data[['Melanocyte']]=NULL
names(data)[names(data)=='T']='T_unclear'
tmp= list_to_matrix(data) %>% make_comb_mat(mode = 'distinct')
tmp2= do.call(paste0,as.data.frame(list_to_matrix(data))) #concatnate the code
tmp3= do.call(paste0,as.data.frame(t(as.data.frame(tmp)))) #how to get the code?

col_name= match(tmp3,tmp2) %>% rownames(list_to_matrix(data))[.] %>% .[1:15]
p2=UpSet(tmp[comb_degree(tmp)>=7],#comb_col = c("black", "black","darkgreen","darkgreen","blue","red")[comb_degree(tmp)], #comb_degree: how many sets share a term
         column_labels =col_name,top_annotation = NULL) 
pdf(file = paste0(output_dir,project,'_round1.5_disorhythm_singlecell_GO_upset.pdf'),width = 10,height = 10)
p2
dev.off()

tmp4= lapply(data,unique) %>% do.call(c,.) %>% table %>% sort(decreasing = T) %>% as.data.frame() %>% dplyr::filter(Freq>6)
colnames(tmp4)[1]='Term'
tmp4$Term= as.character(tmp4$Term)
tmp4$Term[grep('organelle',tmp4$Term)]= 'Organelle assembly'
tmp4$Term[tmp4$Term=='viral process']= 'Viral process'
tmp4$Term[tmp4$Term=='myeloid cell differentiation']= 'Myeloid cell differentiation'
tmp4$Term[tmp4$Term=='regulation of hemopoiesis']= 'Regulation of hemopoiesis'
tmp4$Term[tmp4$Term=='homeostasis of number of cells']= 'Homeostasis of number of cells'
tmp4$Term=factor(tmp4$Term,levels = tmp4$Term)
tmp4$highlight= 'N'
tmp4[grep('myeloid|activation|number|adhesion',tmp4$Term,ignore.case = T),]$highlight= 'Y'
p1=ggplot(aes(y = Freq, x = Term, fill = highlight), data = tmp4) + geom_bar(
  position= "stack",stat="identity") +theme(legend.position="top")+ scale_fill_hue(direction = -1)+
  coord_flip()+xlab('')+
  theme(panel.background = element_rect(fill = "transparent", colour = NA),  plot.background = element_rect(fill = "transparent", colour = NA))+
  theme(text = element_text(size=32),axis.text.y = element_text(size=24),legend.position = 'none')
p1
ggsave(filename = paste0(output_dir,project,'_round1.5_disorhythm_singlecell_GO_bar.pdf'),p1,width = 10,height = 10)

###########heatmap of CD8 genes################
genes= c('Pdcd1','Havcr2','Tigit','Lag3','Ctla4','Btla','Cd101','Il10','Cd160', #'Cd244',
         'Cd28','Tnfrsf9','Cd27','Icos','Cd226','Tnfrsf4','Ccl7',
         "Gzma","Gzmb",'Lamp1',"Prf1","Nkg7","Il2","Ifng","Tnf","Il4ra",'Il21', #"Gnly",
         'Il2ra','Il2rb','Il2rg','Ifngr1','Il7r','Il12rb2','Il18r1','Il18rap','Il1rl1','Il6st','Il10ra','Il12rb1','Il21r',
         'Klf2','Lef1','Tcf7','Tbx21','Irf4','Gata3','Rorc','Batf','Eomes','Foxp1','Foxp3','Tox','Tox2','Maf',
         'Ccr7','Ccr2','Cxcr3','Cxcr4','Cx3cr1','S1pr1','Itga1','Itga4','Itgae','Itgb1','Itgb7','Cd44','Ly6c2','Cxcr5')
mMel_in@meta.data$tmp= mMel_in@meta.data$tmp= interaction(mMel_in@meta.data$timepoint ,mMel_in@meta.data$subclustering_round2)
Avg_exp= AverageExpression(mMel_in,assays = 'RNA',features = genes,group.by = 'tmp') %>%.[[1]]
pdf(file = paste0(output_dir,project,'_CD8_functional_genes_heatmap.pdf'),width = 20,height = 6)
p1=pheatmap::pheatmap(t(Avg_exp),scale = 'column',cluster_rows = F,cluster_cols = F,
                      cellwidth = 15,cellheight = 15,gaps_row = c(4,8),gaps_col = c(9,16,26,39,53),
                      fontsize_col = 15,
                      color = c(colorRampPalette(colors = c("blue","white"))(10),colorRampPalette(colors = c("white","red"))(10)))
dev.off()

p1= dittoBarPlot(mMel_in,'subclustering_round2',group.by = 'timepoint')+
  theme(text = element_text(size=20))
p1
ggsave(filename = paste0(output_dir,project,'_CD8_composition_barplot.pdf'),p1,width = 10,height = 10)




