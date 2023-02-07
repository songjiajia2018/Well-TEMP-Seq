library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(pheatmap)

cre.seu <- function(rds, project, time, type){
    seu <- readRDS(rds)
    seu <- as.matrix(seu)
    seu <- CreateSeuratObject(counts=seu, project=project, min.cells = 5)
    seu$time <- time
    seu$type <- type
    return(seu)
}
time0_total <- cre.seu('H0_total.rds', 'total_0', '0d', 'total')
time1_total <- cre.seu('H1_total.rds', 'total_1', '1d', 'total')
time2_total <- cre.seu('H2_total.rds', 'total_2', '2d', 'total')
time3_total <- cre.seu('H3_total.rds', 'total_3', '3d', 'total')


## create big seurat object for total RNA
group_id <- c('total_time0','total_time1','total_time2','total_time3')
HCL.total <- merge(x=time0_total, y=c(time1_total,time2_total,time3_total), 
                   add.cell.ids = group_id, project = "total")
HCL.total[["percent.mt"]] <- PercentageFeatureSet(HCL.total, pattern = "^MT-")
HCL.total[["percent.rbp"]] <- PercentageFeatureSet(HCL.total, pattern = "^RP[SL]")

filter_pbmc_filter <- subset(HCL.total, subset = nFeature_RNA > 200 & nFeature_RNA <5000 & percent.mt < 15) # 
filter_pbmc_norm <- NormalizeData(filter_pbmc_filter, normalization.method = "LogNormalize", scale.factor = 10000)
filter_pbmc_norm <- FindVariableFeatures(filter_pbmc_norm, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(filter_pbmc_norm)
filter_pbmc_norm <- ScaleData(filter_pbmc_norm, features = all.genes)
filter_pbmc_norm <- RunPCA(filter_pbmc_norm, features = VariableFeatures(object = filter_pbmc_norm))
filter_pbmc_Neighbors <- FindNeighbors(filter_pbmc_norm, dims = 1:10)
filter_pbmc_Nei_cluster <- FindClusters(filter_pbmc_Neighbors, resolution = 0.2)
filter_pbmc_Nei_cluster_umap <- RunUMAP(filter_pbmc_Nei_cluster, dims = 1:10)


split.list <- SplitObject(filter_pbmc_filter, split.by = "time")
#extract cell barcode for total,new,old
time0_cell_barcode_filter <- colnames(split.list[['0d']]) %>% gsub("^total_time0_","",.)
time1_cell_barcode_filter <- colnames(split.list[['1d']]) %>% gsub("^total_time1_","",.)
time2_cell_barcode_filter <- colnames(split.list[['2d']]) %>% gsub("^total_time2_","",.)
time3_cell_barcode_filter <- colnames(split.list[['3d']]) %>% gsub("^total_time3_","",.)
gene_filter <- rownames(split.list[['0d']])

##read new RNA data and filter by total RNA cell barcode
time0_new <- cre.seu('H0_pre_new.rds', 'new_0', '0d', 'new')
time1_new <- cre.seu('H1_pre_new.rds', 'new_1', '1d', 'new')
time2_new <- cre.seu('H2_pre_new.rds', 'new_2', '2d', 'new')
time3_new <- cre.seu('H3_pre_new.rds', 'new_3', '3d', 'new')

time0_new<-time0_new[,time0_cell_barcode_filter]
time1_new<-time1_new[,time1_cell_barcode_filter]
time2_new<-time2_new[,time2_cell_barcode_filter]
time3_new<-time3_new[,time3_cell_barcode_filter]

group_id <- c('new_time0','new_time1','new_time2','new_time3')
HCL.new <- merge(x=time0_new, y=c(time1_new,time2_new,time3_new), add.cell.ids = group_id, project = "new")
HCL.new[["percent.mt"]] <- PercentageFeatureSet(HCL.new, pattern = "^MT-")
HCL.new[["percent.rbp"]] <- PercentageFeatureSet(HCL.new, pattern = "^RP[SL]")

new_filter_pbmc_norm <- NormalizeData(HCL.new, normalization.method = "LogNormalize", scale.factor = 10000)
new_filter_pbmc_norm <- FindVariableFeatures(new_filter_pbmc_norm, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(new_filter_pbmc_norm)
new_filter_pbmc_norm <- ScaleData(new_filter_pbmc_norm, features = all.genes)
new_filter_pbmc_norm <- RunPCA(new_filter_pbmc_norm, features = VariableFeatures(object = new_filter_pbmc_norm))
new_filter_pbmc_Neighbors <- FindNeighbors(new_filter_pbmc_norm, dims = 1:10)
new_filter_pbmc_Nei_cluster <- FindClusters(new_filter_pbmc_Neighbors, resolution = 0.8)
new_filter_pbmc_Nei_cluster_umap <- RunUMAP(new_filter_pbmc_Nei_cluster, dims = 1:10)

### FindMarkers for total RNA
cluster1.markers <- FindMarkers(filter_pbmc_norm, ident.2 ='total_0', ident.1 ='total_1', min.pct = 0.25,logfc.threshold = 0.25,test.use="wilcox")
cluster2.markers <- FindMarkers(filter_pbmc_norm, ident.2 ='total_0', ident.1 ='total_2', min.pct = 0.25,logfc.threshold = 0.25,test.use="wilcox")
cluster3.markers <- FindMarkers(filter_pbmc_norm, ident.2 ='total_0', ident.1 ='total_3', min.pct = 0.25,logfc.threshold = 0.25,test.use="wilcox")

cutoff <- 0.58
cluster1.markers_filter_pos<-subset(cluster1.markers,p_val_adj < 0.05&avg_log2FC > cutoff)
cluster2.markers_filter_pos<-subset(cluster2.markers,p_val_adj < 0.05&avg_log2FC > cutoff)
cluster3.markers_filter_pos<-subset(cluster3.markers,p_val_adj < 0.05&avg_log2FC > cutoff)
cluster1.markers_filter_neg<-subset(cluster1.markers,p_val_adj < 0.05&avg_log2FC < -cutoff)
cluster2.markers_filter_neg<-subset(cluster2.markers,p_val_adj < 0.05&avg_log2FC < -cutoff)
cluster3.markers_filter_neg<-subset(cluster3.markers,p_val_adj < 0.05&avg_log2FC < -cutoff)


pos_DE_table<-union(x=rownames(cluster1.markers_filter_pos),y=rownames(cluster2.markers_filter_pos))
pos_DE_table<-union(x=pos_DE_table,y=rownames(cluster3.markers_filter_pos))
pos_DE_table<-as.data.frame(pos_DE_table)
pos_DE_table<-pos_DE_table[!pos_DE_table$pos_DE_table%>%grepl('MT-',perl=TRUE,.),]
pos_DE_table<-pos_DE_table[!grepl('^RP[SL]',perl=TRUE,pos_DE_table)]
pos_DE_table<-as.data.frame(pos_DE_table)

neg_DE_table<-union(x=rownames(cluster1.markers_filter_neg),y=rownames(cluster2.markers_filter_neg))
neg_DE_table<-union(x=neg_DE_table,y=rownames(cluster3.markers_filter_neg))
neg_DE_table<-as.data.frame(neg_DE_table)
neg_DE_table<-neg_DE_table[!neg_DE_table$neg_DE_table%>%grepl('MT-',perl=TRUE,.),]
neg_DE_table<-neg_DE_table[!grepl('^RP[SL]',perl=TRUE,neg_DE_table)]
neg_DE_table<-as.data.frame(neg_DE_table)

write.table(pos_DE_table, file ="total_posDE_0.58.txt", sep ="\t", quote =FALSE)
write.table(neg_DE_table, file ="total_negDE_0.58.txt", sep ="\t", quote =FALSE)



### FindMarkers for new RNA
cluster1.markers <- FindMarkers(new_filter_pbmc_norm, ident.2 ='new_0', ident.1 ='new_1', min.pct = 0.25,logfc.threshold = 0.25, test.use="t")
cluster2.markers <- FindMarkers(new_filter_pbmc_norm, ident.2 ='new_0', ident.1 ='new_2', min.pct = 0.25,logfc.threshold = 0.25, test.use="t")
cluster3.markers <- FindMarkers(new_filter_pbmc_norm, ident.2 ='new_0', ident.1 ='new_3', min.pct = 0.25,logfc.threshold = 0.25, test.use="t")

cutoff <- 0.58
cluster1.markers_filter_pos<-subset(cluster1.markers,p_val_adj < 0.05&avg_log2FC > cutoff)
cluster2.markers_filter_pos<-subset(cluster2.markers,p_val_adj < 0.05&avg_log2FC > cutoff)
cluster3.markers_filter_pos<-subset(cluster3.markers,p_val_adj < 0.05&avg_log2FC > cutoff)
cluster1.markers_filter_neg<-subset(cluster1.markers,p_val_adj < 0.05&avg_log2FC < -cutoff)
cluster2.markers_filter_neg<-subset(cluster2.markers,p_val_adj < 0.05&avg_log2FC < -cutoff)
cluster3.markers_filter_neg<-subset(cluster3.markers,p_val_adj < 0.05&avg_log2FC < -cutoff)

pos_DE_table<-union(x=rownames(cluster1.markers_filter_pos),y=rownames(cluster2.markers_filter_pos))
pos_DE_table<-union(x=pos_DE_table,y=rownames(cluster3.markers_filter_pos))
pos_DE_table<-as.data.frame(pos_DE_table)
pos_DE_table<-pos_DE_table[!pos_DE_table$pos_DE_table%>%grepl('MT-',perl=TRUE,.),]
pos_DE_table<-pos_DE_table[!grepl('^RP[SL]',perl=TRUE,pos_DE_table)]
pos_DE_table<-as.data.frame(pos_DE_table)

neg_DE_table<-union(x=rownames(cluster1.markers_filter_neg),y=rownames(cluster2.markers_filter_neg))
neg_DE_table<-union(x=neg_DE_table,y=rownames(cluster3.markers_filter_neg))
neg_DE_table<-as.data.frame(neg_DE_table)
neg_DE_table<-neg_DE_table[!neg_DE_table$neg_DE_table%>%grepl('MT-',perl=TRUE,.),]
neg_DE_table<-neg_DE_table[!grepl('^RP[SL]',perl=TRUE,neg_DE_table)]
neg_DE_table<-as.data.frame(neg_DE_table)


write.table(pos_DE_table, file ="new_posDE_0.58.txt", sep ="\t", quote =FALSE, row.names = FALSE, col.names = FALSE)
write.table(neg_DE_table, file ="new_negDE_0.58.txt", sep ="\t", quote =FALSE, row.names = FALSE, col.names = FALSE)