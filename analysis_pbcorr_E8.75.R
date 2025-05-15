library(Seurat)

# pseudo-bulk correlation analysis of mouse E8.75 anterior and posterior cell types

obj = readRDS(file='../objects/GSE123046_E8.75_AP.rds')

obj@meta.data$group <- paste0(obj@meta.data$Sample, "_", obj@meta.data$CellType)

obj.avg <- AverageExpression(obj, assays = "RNA", group.by = "group", return.seurat = TRUE)
obj.avg <- FindVariableFeatures(obj.avg)
obj.avg <- ScaleData(obj.avg)
obj.avg <- RunPCA(obj.avg, npcs = 10)

mat <- obj.avg@reductions$pca@cell.embeddings

cor_mat <- cor(t(mat))  
