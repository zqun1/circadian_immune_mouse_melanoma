#########general########
#plot heatmap: give a seu and its DEGs, draw a heatmap of avg expression according to the group.by.
#DEGs should be calculated by findallmarkers
plotmarkers_heatmap<- function(seu,data.frame,group.by,name='log2FC', cluster_rows = F,cluster_columns = F,show_column_names=T,
                               show_colnames = T,annotation_col = T,annotation_row = T,show_row_names = F,annotation_legend=F,
                               annotation_colors = c('dittoSeq','scales'),label_genes=NULL,
                               cell_width = NULL,cell_height = NULL,column_names_gp = gpar(fontsize = 12),
                               row_names_gp = gpar(fontsize = 12),...){
  library(scales)
  library(dittoSeq)
  library(ComplexHeatmap)
  
  if(is.null(levels(seu@meta.data[,group.by]))){
    seu@meta.data[,group.by]= factor(seu@meta.data[,group.by],levels= unique(seu@meta.data[,group.by]))
  }
  cell_levels= levels(seu@meta.data[,group.by])
  data.frame$cluster= factor(data.frame$cluster,levels=cell_levels)
  features= data.frame %>% group_by(cluster) %>% dplyr::arrange(-avg_log2FC,.by_group = T) %>% .$gene
  table1=AverageExpression(seu,features = features,group.by = group.by)
  table2=t(scale(t(table1$RNA)))
  table2=table2[features,]

  if(tolower(annotation_colors)=='dittoseq'){annotation_colors= dittoSeq::dittoColors()[1:length(cell_levels)]}
  else{annotation_colors= hue_pal()(length(cell_levels))}
  if(annotation_col==T){
    col_use_top= annotation_colors
    names(col_use_top)= cell_levels
    top_annotation= HeatmapAnnotation(show_legend = annotation_legend,Cells= cell_levels,col = list(Cells=col_use_top),show_annotation_name = F)
  }else{top_annotation = NULL}
  if(annotation_row==T){
    col_use_left= as.factor(data.frame$cluster)
    levels(col_use_left)= annotation_colors
    names(col_use_left)= data.frame$cluster
    right_annotation=rowAnnotation(show_legend = annotation_legend,Cells= data.frame$cluster, #side bar
                                   col=list(Cells=col_use_left) #color
                                   )
  }else{right_annotation = NULL}
  
  library(circlize)
  if(!is.null(label_genes)){
    right_label_annotation= rowAnnotation(link = anno_mark(at = which(rownames(table2)%in%label_genes), #custome label genes
                                                                       labels = rownames(table2)[which(rownames(table2)%in%label_genes)]))
  }else{right_label_annotation = NULL}
  
  if(!is.null(cell_width)){width = ncol(table2)*unit(cell_width, "mm")}else{width = NULL}
  if(!is.null(cell_height)){height = nrow(table2)*unit(cell_height, "mm")}else{height= NULL}

  p1= Heatmap(table2,name=name,col = colorRamp2(c(min(table2),0,max(table2)), c("blue", "white", "red")),show_column_names=show_column_names,
              cluster_rows = cluster_rows,cluster_columns = cluster_columns,show_row_names = show_row_names,
              top_annotation = top_annotation,right_annotation = right_annotation,
              width = width,height = height,row_names_gp = gpar(fontsize = 12),
              column_names_gp = gpar(fontsize = 12),...) + right_label_annotation
  
  return(p1)
}

#for bulk
plotmarkers_heatmap_bulk<- function(data.frame,metadata,group.by,name='log2FC', cluster_rows = F,cluster_columns = F,show_column_names=T,
                               show_colnames = T,annotation_col = T,annotation_row = F,show_row_names = F,annotation_legend=F,
                               annotation_colors = NULL,label_genes=NULL,
                               cell_width = NULL,cell_height = NULL,column_names_gp = gpar(fontsize = 12),
                               row_names_gp = gpar(fontsize = 12),...){
  library(scales)
  library(dittoSeq)
  library(ComplexHeatmap)
  
  if(is.null(levels(metadata[,group.by]))){
    metadata[,group.by]= factor(metadata[,group.by],levels= unique(metadata[,group.by]))
  }
  cell_levels= levels(metadata[,group.by])
  table2=t(scale(t(data.frame)))
  
  if(!is.null(annotation_colors)){annotation_colors= annotation_colors[1:length(cell_levels)]}
  else{annotation_colors= dittoSeq::dittoColors()[1:length(cell_levels)]}
  if(annotation_col==T){
    col_use_top= annotation_colors
    names(col_use_top)= cell_levels
    top_annotation= HeatmapAnnotation(show_legend = annotation_legend,Cells= metadata[,group.by],
                                      col = list(Cells=col_use_top[metadata[,group.by]]),show_annotation_name = F)
  }else{top_annotation = NULL}
  if(annotation_row==T){
    message('not developed yet')
  }else{right_annotation = NULL}
  
  library(circlize)
  if(!is.null(label_genes)){
    right_label_annotation= rowAnnotation(link = anno_mark(at = which(rownames(table2)%in%label_genes), #custome label genes
                                                           labels = rownames(table2)[which(rownames(table2)%in%label_genes)]))
  }else{right_label_annotation = NULL}
  
  if(!is.null(cell_width)){width = ncol(table2)*unit(cell_width, "mm")}else{width = NULL}
  if(!is.null(cell_height)){height = nrow(table2)*unit(cell_height, "mm")}else{height= NULL}
  
  p1= Heatmap(table2,name=name,col = colorRamp2(c(min(table2),0,max(table2)), c("blue", "white", "red")),
              show_column_names=show_column_names,show_row_names = show_row_names,
              cluster_rows = cluster_rows,cluster_columns = cluster_columns,
              top_annotation = top_annotation,right_annotation = right_annotation,
              width = width,height = height,row_names_gp = gpar(fontsize = 12),
              column_names_gp = gpar(fontsize = 12),...) + right_label_annotation
  
  return(p1)
}

`%notin%` <- Negate(`%in%`) #very slow

# convert Ensembl to Symbols and calculate the average as the expression for duplicated Symbols
# library(AnnotationHub)
# library(ensembldb)
# hub <- AnnotationHub()
# dm <- query(hub, c("EnsDb", "mmusculus")) #AH78783 AH79689 AH83216 AH89180 AH89426 | Ensembl 103 EnsDb for Homo sapiens
# geneIDs2= dm[['AH109655']] #
# save(geneIDs2,file = '../geneIDs.RData')
Mouse_Ensembl_to_Symbol <- function(expression_matrix){
  expression_matrix= as.data.frame(expression_matrix)
  ensembl.genes= rownames(expression_matrix)
  geneIDs= read.csv('../Collaboration_Chen_melanoma/mouse_ensembl_length.csv',row.names = 1)
  table(duplicated(geneIDs$symbol))
  table(geneIDs$symbol=='') #has duplicated IDs and blank IDs
  expression_matrix$ID= geneIDs[match(rownames(expression_matrix),geneIDs$gene_id),]$symbol
  #use mean as the duplicated gene expression levels
  expression_matrix <- expression_matrix%>%
    dplyr::filter(ID!='') %>% 
    dplyr::group_by(ID) %>%
    dplyr::summarise_all(mean)  %>%
    as.data.frame() %>% 
    column_to_rownames('ID')
  return(expression_matrix)
}


# Basic function to convert mouse to human gene names
convertMouseGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  genesV2 = genesV2[match(x,genesV2$MGI.symbol),]
  return(genesV2)
}
# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  genesV2 = genesV2[match(x,genesV2$hgnc_symbol),]
  return(genesV2)
}

#find duplicated rows
Find_duplicated<- function(dataframe,colnames=1){
  # duplicated_values= list()
  duplicated_rows= list()
  for(i in colnames){
    duplicated_values = dataframe[duplicated(dataframe[,i]),i]
    duplicated_rows[[i]] = which(dataframe[,i] %in% duplicated_values)
  }
  duplicated_rows= unlist(duplicated_rows) %>% unique()
  return(dataframe[duplicated_rows,])
}

plotFractionPvalue<- function(metadata,name,parameter='subclustering_round1', compare= T, #whether do statistical analysis
                              psize=3.5,fontsize=13,width=10,height=10,origin_levels= NULL, paired= F){ 
  #compare the percentage of cells from less than 3 different origins.
  library(ggpubr)
  table1= as.data.frame(table(metadata$Sample_origin,metadata$patientID,metadata[,parameter]))
  colnames(table1)=c("origin","ID","type","freq")
  table2= table(metadata$Sample_origin,metadata$patientID)
  for (i in 1:nrow(table1)) {table1$percent[i]= 100*table1$freq[i]/table2[table1$origin[i],table1$ID[i]]}
  if(is.null(origin_levels))  origin_levels=unique(metadata$Sample_origin)
  table1$origin=factor(table1$origin,levels = origin_levels)
  library(scales)
  color= hue_pal(direction = -1)(length(origin_levels))
  combinations <- combn(origin_levels,m=2)
  my_comparisons= list()
  for(i in 1:ncol(combinations)){my_comparisons[[i]]= combinations[,i]}
  p1 <- ggpaired(table1, x = "origin", y = "percent",color = "origin", id = "ID",palette = color,
                 line.color = "gray", line.size = 0.4,facet.by = "type", short.panel.labs = FALSE)+
    facet_grid(cols =  vars(type)) +  theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5,size = fontsize), legend.position= "none") + xlab("") + ylab("Percentage")+
    theme(text = element_text(size=fontsize))
  if(compare){
    p1=p1+stat_compare_means(comparisons = my_comparisons,label = "p.format", paired = paired, method = "wilcox.test",size=psize)
  }
  ggsave(p1,filename = paste0(output_dir,project,"_",name,"_",parameter,"_fraction_NT_pval_wilcox.pdf"),width = width, height = height)
  
  p1 <- ggpaired(table1, x = "origin", y = "percent",color = "origin", id = "ID",palette = color,
                 line.color = "gray", line.size = 0.4,facet.by = "type", short.panel.labs = FALSE)+
    facet_grid(cols =  vars(type)) +  theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5,size = fontsize), legend.position= "none") + xlab("")+ ylab("Percentage") +
    theme(text = element_text(size=fontsize))
  if(compare){
    p1=p1+stat_compare_means(comparisons = my_comparisons,label = "p.format", paired = paired, method = "wilcox.test",size=psize)
  }
  ggsave(p1,filename = paste0(output_dir,project,"_",name,"_",parameter,"_fraction_NT_pval_t.pdf"),width = width, height = height)
}



##############seurat#################
#DimPlot_modified modifies the labeling of DimPlot to make it read-friendly. Use it as DimPlot
DimPlot_modified<- function(seu, prefix='cluster',group.by=NULL,label=T,label.box=T,
                            cols=dittoSeq::dittoColors(),legend_col=3,legend_fontsize=6,...){
  #if group.by=NULL, then use the active idents
  seu[["ident"]] <- Idents(object = seu)
  orig.groups <- group.by
  group.by <- group.by %||% "ident"
  
  lengths= length(unique(seu@meta.data[,group.by]))
  tmp<- paste0(prefix, seq(from=0,length.out=lengths))
  seu@meta.data[,group.by]= as.factor(seu@meta.data[,group.by])
  seu@meta.data$tmp <- factor(tmp[seu@meta.data[,group.by]], levels = tmp)
  index_cluster <- paste0(prefix, seq(from=0,length.out=lengths)," ",levels(seu@meta.data[,group.by]))
  seu@meta.data$cluster_index <- 
    factor(index_cluster[seu@meta.data[,group.by]], levels = index_cluster)
  p <-
    DimPlot(seu, group.by = "tmp", label = label, label.box = label.box, 
            cols = cols,...) +
    scale_color_discrete(
      type = cols,
      labels = levels(seu@meta.data$cluster_index)) +
    theme(plot.title = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(size = legend_fontsize),
          ) +
    guides(color = guide_legend(
      ncol = legend_col, #how many columns you want for legend
      override.aes=list(size = 3))) +
    scale_fill_manual(values = rep("white",lengths))
  return(p)
}

#Read10X_path reads vectors of matrix.dir, features.dir, and barcodes.dir that specify the location of each individual file.
#Compared to Read10X, it does not use the mother dir and does not require the names to be exactly matrix.tsv, features.tsv, barcodes.tsv
Read10X_path <- function (matrix.dir, features.dir, barcodes.dir, gene.column = 2, cell.column = 1, unique.features = TRUE, 
                              strip.suffix = FALSE) #genes.tsv.gz and features.tsv.gz problem
{
  library(Matrix)
  full.data <- list()
  for (i in seq(length(matrix.dir))) {
    barcode.loc <- barcodes.dir[i]
    features.loc <- features.dir[i]
    matrix.loc <- matrix.dir[i]
    data <- readMM(file = matrix.loc)
    cell.barcodes <- read.table(file = barcode.loc, header = FALSE, 
                                sep = "\t", row.names = NULL)
    if (ncol(x = cell.barcodes) > 1) {
      cell.names <- cell.barcodes[, cell.column]
    }
    else {
      cell.names <- readLines(con = barcode.loc)
    }
    if (all(grepl(pattern = "\\-1$", x = cell.names)) & 
        strip.suffix) {
      cell.names <- as.vector(x = as.character(x = sapply(X = cell.names, 
                                                          FUN = ExtractField, field = 1, delim = "-")))
    }
    if (is.null(x = names(x = matrix.dir))) {
      if (length(x = matrix.dir) < 2) {
        colnames(x = data) <- cell.names
      }
      else {
        colnames(x = data) <- paste0(i, "_", cell.names)
      }
    }
    else {
      colnames(x = data) <- paste0(names(x = matrix.dir)[i], 
                                   "_", cell.names)
    }
    feature.names <- read.delim(file = features.loc, header = FALSE, 
                                stringsAsFactors = FALSE)
    if (any(is.na(x = feature.names[, gene.column]))) {
      warning("Some features names are NA. Replacing NA names with ID from the opposite column requested", 
              call. = FALSE, immediate. = TRUE)
      na.features <- which(x = is.na(x = feature.names[, 
                                                       gene.column]))
      replacement.column <- ifelse(test = gene.column == 
                                     2, yes = 1, no = 2)
      feature.names[na.features, gene.column] <- feature.names[na.features, 
                                                               replacement.column]
    }
    if (unique.features) {
      fcols = ncol(x = feature.names)
      if (fcols < gene.column) {
        stop(paste0("gene.column was set to ", gene.column, 
                    " but feature.tsv.gz (or genes.tsv) only has ", 
                    fcols, " columns.", " Try setting the gene.column argument to a value <= to ", 
                    fcols, "."))
      }
      rownames(x = data) <- make.unique(names = feature.names[, 
                                                              gene.column])
    }
    if (ncol(x = feature.names) > 2) {
      data_types <- factor(x = feature.names$V3)
      lvls <- levels(x = data_types)
      if (length(x = lvls) > 1 && length(x = full.data) == 
          0) {
        message("10X data contains more than one type and is being returned as a list containing matrices of each type.")
      }
      expr_name <- "Gene Expression"
      if (expr_name %in% lvls) {
        lvls <- c(expr_name, lvls[-which(x = lvls == 
                                           expr_name)])
      }
      data <- lapply(X = lvls, FUN = function(l) {
        return(data[data_types == l, , drop = FALSE])
      })
      names(x = data) <- lvls
    }
    else {
      data <- list(data)
    }
    full.data[[length(x = full.data) + 1]] <- data
  }
  return(full.data)
}

Read10X_modified <- function (data.dir, gene.column = 2, cell.column = 1, unique.features = TRUE, 
          strip.suffix = FALSE) #genes.tsv.gz and features.tsv.gz problem
{
  library(Matrix)
  full.data <- list()
  for (i in seq_along(along.with = data.dir)) {
    run <- data.dir[i]
    if (!dir.exists(paths = run)) {
      stop("Directory provided does not exist")
    }
    barcode.loc <- file.path(run, "barcodes.tsv")
    gene.loc <- file.path(run, "genes.tsv")
    features.loc <- file.path(run, "features.tsv")
    matrix.loc <- file.path(run, "matrix.mtx")
    pre_ver_3 <- file.exists(gene.loc)
    if (!pre_ver_3) {
      addgz <- function(s) {
        return(paste0(s, ".gz"))
      }
      barcode.loc <- addgz(s = barcode.loc)
      matrix.loc <- addgz(s = matrix.loc)
      gene.loc<- addgz(s = gene.loc)
      features.loc= gene.loc #if not pre_ver3, genes.tsv.gz instead of features.tsv.gz is present. The same???
    }
    if (!file.exists(barcode.loc)) {
      stop("Barcode file missing. Expecting ", basename(path = barcode.loc))
    }
    if (!pre_ver_3 && !file.exists(features.loc)) {
      stop("Gene name or features file missing. Expecting ", 
           basename(path = features.loc))
    }
    if (!file.exists(matrix.loc)) {
      stop("Expression matrix file missing. Expecting ", 
           basename(path = matrix.loc))
    }
    data <- readMM(file = matrix.loc)
    cell.barcodes <- read.table(file = barcode.loc, header = FALSE, 
                                sep = "\t", row.names = NULL)
    if (ncol(x = cell.barcodes) > 1) {
      cell.names <- cell.barcodes[, cell.column]
    }
    else {
      cell.names <- readLines(con = barcode.loc)
    }
    if (all(grepl(pattern = "\\-1$", x = cell.names)) & 
        strip.suffix) {
      cell.names <- as.vector(x = as.character(x = sapply(X = cell.names, 
                                                          FUN = ExtractField, field = 1, delim = "-")))
    }
    if (is.null(x = names(x = data.dir))) {
      if (length(x = data.dir) < 2) {
        colnames(x = data) <- cell.names
      }
      else {
        colnames(x = data) <- paste0(i, "_", cell.names)
      }
    }
    else {
      colnames(x = data) <- paste0(names(x = data.dir)[i], 
                                   "_", cell.names)
    }
    feature.names <- read.delim(file = ifelse(test = pre_ver_3, 
                                              yes = gene.loc, no = features.loc), header = FALSE, 
                                stringsAsFactors = FALSE)
    if (any(is.na(x = feature.names[, gene.column]))) {
      warning("Some features names are NA. Replacing NA names with ID from the opposite column requested", 
              call. = FALSE, immediate. = TRUE)
      na.features <- which(x = is.na(x = feature.names[, 
                                                       gene.column]))
      replacement.column <- ifelse(test = gene.column == 
                                     2, yes = 1, no = 2)
      feature.names[na.features, gene.column] <- feature.names[na.features, 
                                                               replacement.column]
    }
    if (unique.features) {
      fcols = ncol(x = feature.names)
      if (fcols < gene.column) {
        stop(paste0("gene.column was set to ", gene.column, 
                    " but feature.tsv.gz (or genes.tsv) only has ", 
                    fcols, " columns.", " Try setting the gene.column argument to a value <= to ", 
                    fcols, "."))
      }
      rownames(x = data) <- make.unique(names = feature.names[, 
                                                              gene.column])
    }
    if (ncol(x = feature.names) > 2) {
      data_types <- factor(x = feature.names$V3)
      lvls <- levels(x = data_types)
      if (length(x = lvls) > 1 && length(x = full.data) == 
          0) {
        message("10X data contains more than one type and is being returned as a list containing matrices of each type.")
      }
      expr_name <- "Gene Expression"
      if (expr_name %in% lvls) {
        lvls <- c(expr_name, lvls[-which(x = lvls == 
                                           expr_name)])
      }
      data <- lapply(X = lvls, FUN = function(l) {
        return(data[data_types == l, , drop = FALSE])
      })
      names(x = data) <- lvls
    }
    else {
      data <- list(data)
    }
    full.data[[length(x = full.data) + 1]] <- data
  }
  list_of_data <- list()
  for (j in 1:length(x = full.data[[1]])) {
    list_of_data[[j]] <- do.call(cbind, lapply(X = full.data, 
                                               FUN = `[[`, j))
    list_of_data[[j]] <- as(object = list_of_data[[j]], 
                            Class = "dgCMatrix")
  }
  names(x = list_of_data) <- names(x = full.data[[1]])
  if (length(x = list_of_data) == 1) {
    return(list_of_data[[1]])
  }
  else {
    return(list_of_data)
  }
}

PlotORA<- function(DEG,nterm=20,cell="?",only.pos=F,weighted_log2FC=0.5){ #weighted_log2FC
  cat("Now ",cell)
  DEGup=DEG[DEG$weighted_log2FC>weighted_log2FC,]
  genelist= DEGup$gene
  EC_out.GO= enrichGO(genelist,OrgDb = 'org.Hs.eg.db',keyType = 'SYMBOL',ont = 'BP')
  EC_out.GO.sim= clusterProfiler::simplify(EC_out.GO,cutoff=0.7)  
  if(nrow(fortify(EC_out.GO.sim))){ ##May have result but no enrichedResult. fortify to convert the enrichedResult to a dataframe. 
    p1=clusterProfiler::dotplot(EC_out.GO.sim,showCategory=nterm)+
      scale_y_discrete(labels = function(x) str_wrap(x, width = 81)) #forcce break line if y label length is >81
    ggsave(p1,filename = paste0(output_dir,project,"_",cell,"_up_GOBP.pdf"),width = 10,height = 10)
  }
  tmp= genelist %>% bitr(fromType = "SYMBOL",toType = "ENTREZID",OrgDb = 'org.Hs.eg.db') 
  genelist2= bitr_kegg(tmp$ENTREZID,fromType = "ncbi-geneid",toType = "kegg",organism = 'hsa')
  EC_out.KEGG= enrichKEGG(genelist2$kegg,organism = "hsa",keyType ="kegg")
  if(nrow(fortify(EC_out.KEGG))){
    p1=clusterProfiler::dotplot(EC_out.KEGG,showCategory=nterm)
    ggsave(p1,filename = paste0(output_dir,project,"_",cell,"_up_KEGG.pdf"),width = 10,height = 10)
  }
  EC_out.hallmark= enricher(genelist, minGSSize = 10, TERM2GENE = geneset.use[,c(1,2)])
  if(nrow(EC_out.hallmark)){
    p1=clusterProfiler::dotplot(EC_out.hallmark,showCategory=30)
    ggsave(p1,filename = paste0(output_dir,project,"_",cell,"_up_Hallmarks.pdf"),width = 10,height = 10)
  }
  
  if(only.pos==FALSE){
    DEGdown=DEG[DEG$weighted_log2FC< -weighted_log2FC,]
    genelist= DEGdown$gene
    EC_out.GO= enrichGO(genelist,OrgDb = 'org.Hs.eg.db',keyType = 'SYMBOL',ont = 'BP')
    EC_out.GO.sim= clusterProfiler::simplify(EC_out.GO,cutoff=0.7)
    if(nrow(fortify(EC_out.GO.sim))){
      p1=clusterProfiler::dotplot(EC_out.GO.sim,showCategory=nterm)+
        scale_y_discrete(labels = function(x) str_wrap(x, width = 81)) #forcce break line if y label length is >81
      ggsave(p1,filename = paste0(output_dir,project,"_",cell,"_down_GOBP.pdf"),width = 10,height = 10)
    }
    tmp= genelist %>% bitr(fromType = "SYMBOL",toType = "ENTREZID",OrgDb = 'org.Hs.eg.db') 
    genelist2= bitr_kegg(tmp$ENTREZID,fromType = "ncbi-geneid",toType = "kegg",organism = 'hsa')
    EC_out.KEGG= enrichKEGG(genelist2$kegg,organism = "hsa",keyType ="kegg")
    if(nrow(fortify(EC_out.KEGG))){
      p1=clusterProfiler::dotplot(EC_out.KEGG,showCategory=nterm)
      ggsave(p1,filename = paste0(output_dir,project,"_",cell,"_down_KEGG.pdf"),width = 10,height = 10)
    }
    EC_out.hallmark= enricher(genelist, minGSSize = 10, TERM2GENE = geneset.use[,c(1,2)],pvalueCutoff = 0.05,qvalueCutoff = 0.2)
    if(nrow(fortify(EC_out.hallmark))){
      p1=clusterProfiler::dotplot(EC_out.hallmark,showCategory=30)
      ggsave(p1,filename = paste0(output_dir,project,"_",cell,"_down_Hallmarks.pdf"),width = 10,height = 10)
    }
  }
}
PlotORA_avg<- function(DEG,nterm=20,cell="?",only.pos=F,avg_log2FC=0.5,GO_only= FALSE,species='human',keyType='SYMBOL'){ #avg_log2FC
  library(GSEABase)
  library(msigdbr)
  library(clusterProfiler)
  cat("Now ",cell)
  DEGup=DEG[DEG$avg_log2FC>avg_log2FC,]
  genelist= DEGup$gene
  if(tolower(species)=='human'){
    OrgDb= 'org.Hs.eg.db'
    organism='hsa'
  }else if(tolower(species)=='mouse'){
    OrgDb= 'org.Mm.eg.db'
    organism='mmu'
  }

  EC_out.GO= enrichGO(genelist,OrgDb = OrgDb,keyType = keyType,ont = 'BP')
  EC_out.GO.sim= clusterProfiler::simplify(EC_out.GO,cutoff=0.7)  
  if(nrow(fortify(EC_out.GO.sim))){ ##May have result but no enrichedResult. fortify to convert the enrichedResult to a dataframe. 
    p1=clusterProfiler::dotplot(EC_out.GO.sim,showCategory=nterm)+
      scale_y_discrete(labels = function(x) str_wrap(x, width = 81)) #forcce break line if y label length is >81
    ggsave(p1,filename = paste0(output_dir,project,"_",cell,"_up_GOBP.pdf"),width = 10,height = 10)
  }
  
  if(!GO_only){
    tmp= genelist %>% bitr(fromType = keyType,toType = "ENTREZID",OrgDb = OrgDb)
    genelist2= bitr_kegg(tmp$ENTREZID,fromType = "ncbi-geneid",toType = "kegg",organism = organism)
    EC_out.KEGG= enrichKEGG(genelist2$kegg,organism = organism,keyType ="kegg")
    if(nrow(fortify(EC_out.KEGG))){
      p1=clusterProfiler::dotplot(EC_out.KEGG,showCategory=nterm)
      ggsave(p1,filename = paste0(output_dir,project,"_",cell,"_up_KEGG.pdf"),width = 10,height = 10)
    }
    EC_out.hallmark= enricher(genelist, minGSSize = 10, TERM2GENE = geneset.use[,c(1,2)])
    if(nrow(EC_out.hallmark)){
      p1=clusterProfiler::dotplot(EC_out.hallmark,showCategory=30)
      ggsave(p1,filename = paste0(output_dir,project,"_",cell,"_up_Hallmarks.pdf"),width = 10,height = 10)
    }
  }
  
  if(only.pos==FALSE){
    DEGdown=DEG[DEG$avg_log2FC< -avg_log2FC,]
    genelist= DEGdown$gene
    EC_out.GO= enrichGO(genelist,OrgDb = OrgDb,keyType = keyType,ont = 'BP')
    EC_out.GO.sim= clusterProfiler::simplify(EC_out.GO,cutoff=0.7)
    if(nrow(fortify(EC_out.GO.sim))){
      p1=clusterProfiler::dotplot(EC_out.GO.sim,showCategory=nterm)+
        scale_y_discrete(labels = function(x) str_wrap(x, width = 81)) #forcce break line if y label length is >81
      ggsave(p1,filename = paste0(output_dir,project,"_",cell,"_down_GOBP.pdf"),width = 10,height = 10)
    }
    if(!GO_only){
      tmp= genelist %>% bitr(fromType = keyType,toType = "ENTREZID",OrgDb = OrgDb) 
      genelist2= bitr_kegg(tmp$ENTREZID,fromType = "ncbi-geneid",toType = "kegg",organism = organism)
      EC_out.KEGG= enrichKEGG(genelist2$kegg,organism = organism,keyType ="kegg")
      if(nrow(fortify(EC_out.KEGG))){
        p1=clusterProfiler::dotplot(EC_out.KEGG,showCategory=nterm)
        ggsave(p1,filename = paste0(output_dir,project,"_",cell,"_down_KEGG.pdf"),width = 10,height = 10)
      }
      EC_out.hallmark= enricher(genelist, minGSSize = 10, TERM2GENE = geneset.use[,c(1,2)],pvalueCutoff = 0.05,qvalueCutoff = 0.2)
      if(nrow(fortify(EC_out.hallmark))){
        p1=clusterProfiler::dotplot(EC_out.hallmark,showCategory=30)
        ggsave(p1,filename = paste0(output_dir,project,"_",cell,"_down_Hallmarks.pdf"),width = 10,height = 10)
      }
    }
  }
}

############general#########
threetile=function(x){ifelse(x>quantile(x,0.6666),"high",ifelse(x>quantile(x,0.3333),"medium","low"))}

capitalize_1st <- function(charcter_string,original=c('upper','lower')){
  if(tolower(original)=='lower'){
    CapStr <- function(y) {
      c <- strsplit(y, " ")[[1]]
      paste(toupper(substring(c, 1,1)), substring(c, 2),
            sep="", collapse=" ")
    }
    sapply(charcter_string, CapStr)
  } else if(tolower(original)=='upper'){
    CapStr <- function(y) {
      c <- strsplit(y, " ")[[1]]
      paste(substring(c, 1,1), tolower(substring(c, 2)),
            sep="", collapse=" ")
    }
    sapply(charcter_string, CapStr)
  } else {
    stop('please specify whether your original string is upper or lower cases')
  }
}

capitalize_1st_in_mix <- function(charcter_string){
  charcter_string= tolower(charcter_string)
  CapStr <- function(y) {
    c <- strsplit(y, " ")[[1]]
    letter_location= gregexpr(pattern = '[a-zA-Z]', text = c,ignore.case = T)
    letter_location=letter_location[[1]][1]
    if(letter_location==1){
      paste(toupper(substring(c, 1,1)), substring(c, 2),
            sep="", collapse=" ")
    }else{
      paste(substring(c, 1,letter_location-1),toupper(substring(c, letter_location,letter_location)), substring(c, letter_location+1),
            sep="", collapse=" ")
    }}
  sapply(charcter_string, CapStr)
}
