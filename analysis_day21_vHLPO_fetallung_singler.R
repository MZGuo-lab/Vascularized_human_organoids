# ============================================
# Load libraries
# ============================================
library(Seurat)
library(SingleR)
library(SingleCellExperiment)
library(ggplot2)

# ============================================
# Prepare reference
# ============================================

ref = readRDS(file="../objects/fetallung_scRNAseq.rds")

DefaultAssay(ref) <- "RNA"
ref.sce <- as.SingleCellExperiment(ref)

# ============================================
# Prepare query
# ============================================

query = readRDS(file="../objects/Day21.obj.rds")

query = subset(query, DataID=="Day21_vHLPO" & Celltype != "Low-quality")
query@meta.data = droplevels(query@meta.data)

DefaultAssay(query) <- "RNA"
query_use <- NormalizeData(query)
query.sce <- as.SingleCellExperiment(query_use)

# ============================================
# Run SingleR predictions
# ============================================

pred_broad_type <- SingleR(
  test = query.sce,
  ref = ref.sce,
  labels = as.character(ref.sce@colData$broad_celltype),
  de.method = "wilcox"
)

# ============================================
# Assign predictions to query metadata
# ============================================
all(rownames(query_use@meta.data) == rownames(pred_broad_type))

query_use@meta.data$singleR_broadtype <- pred_broad_type[rownames(query_use@meta.data), "labels"]

query_use@meta.data$singleR_lineage <- "Mesen&PNS"
query_use@meta.data$singleR_lineage[query_use@meta.data$singleR_broadtype %in% c("Distal epithelial", "Proximal epithelial")] <- "Epithelial"
query_use@meta.data$singleR_lineage[query_use@meta.data$singleR_broadtype %in% c("Lymph endothelial", "Vas endothelial")] <- "Endothelial"
query_use@meta.data$singleR_lineage[query_use@meta.data$singleR_broadtype %in% c("T & ILC", "B", "NK", "Other myeloid", "Meg-ery")] <- "Immune"

query_use@meta.data$Celltype2 <- paste0(query_use@meta.data$Celltype, "-", query_use@meta.data$seurat_clusters)

# ============================================
# Confusion matrix: Day21 vHLPO clusters vs SingleR predicted lineages
# ============================================
viz <- table(query_use@meta.data$Celltype2, query_use@meta.data$singleR_lineage)

vizp = viz/rowSums(viz)

print(vizp)


###########################################################################################
## Mesen & PNS
###########################################################################################

ref_use <- subset(ref, broad_celltype %in% c("Fibroblast", "Mesothelial", "Myofibro & SMC", "Chondrocyte", "PNS"))
ref_use@meta.data <- droplevels(ref_use@meta.data)
DefaultAssay(ref_use) <- "RNA"
ref.sce <- as.SingleCellExperiment(ref_use)

query_use = subset(query, Celltype %in% c("FB", "Pericyte", "pFB"))
query_use@meta.data = droplevels(query_use@meta.data)
query_use@meta.data <- droplevels(query_use@meta.data)

DefaultAssay(query_use) <- "RNA"
query_use <- NormalizeData(query_use)
query.sce <- as.SingleCellExperiment(query_use)

pred_celltype <- SingleR(
  test = query.sce,
  ref = ref.sce,
  labels = as.character(ref.sce@colData$new_celltype),
  de.method = "wilcox"
)

query_use@meta.data$single_celltype <- pred_celltype[rownames(query_use@meta.data), "labels"]
query_use@meta.data$Celltype2 <- paste0(query_use@meta.data$Celltype, "-", query_use@meta.data$seurat_clusters)

viz <- table(query_use@meta.data$Celltype2, query_use@meta.data$single_celltype)

vizp = viz/rowSums(viz)

print(vizp)


###########################################################################################
## Epi
###########################################################################################

ref_use <- subset(ref, broad_celltype %in% c("Distal epithelial", "Proximal epithelial"))
ref_use@meta.data <- droplevels(ref_use@meta.data)
DefaultAssay(ref_use) <- "RNA"
ref.sce <- as.SingleCellExperiment(ref_use)

query_use = subset(query, Celltype %in% c("Progenitor1", "Progenitor2"))
query_use@meta.data = droplevels(query_use@meta.data)
query_use@meta.data <- droplevels(query_use@meta.data)

DefaultAssay(query_use) <- "RNA"
query_use <- NormalizeData(query_use)
query.sce <- as.SingleCellExperiment(query_use)

pred_celltype <- SingleR(
  test = query.sce,
  ref = ref.sce,
  labels = as.character(ref.sce@colData$new_celltype),
  de.method = "wilcox"
)

query_use@meta.data$single_celltype <- pred_celltype[rownames(query_use@meta.data), "labels"]
query_use@meta.data$Celltype2 <- paste0(query_use@meta.data$Celltype, "-", query_use@meta.data$seurat_clusters)

viz <- table(query_use@meta.data$Celltype2, query_use@meta.data$single_celltype)

vizp = viz/rowSums(viz)

print(vizp)


###########################################################################################
## EC
###########################################################################################

ref_use <- subset(ref, broad_celltype %in% c("Vas endothelial", "Lymph endothelial"))
ref_use@meta.data <- droplevels(ref_use@meta.data)
DefaultAssay(ref_use) <- "RNA"
ref.sce <- as.SingleCellExperiment(ref_use)

query_use = subset(query, Celltype %in% c("EC"))
query_use@meta.data = droplevels(query_use@meta.data)
query_use@meta.data <- droplevels(query_use@meta.data)

DefaultAssay(query_use) <- "RNA"
query_use <- NormalizeData(query_use)
query.sce <- as.SingleCellExperiment(query_use)

pred_celltype <- SingleR(
  test = query.sce,
  ref = ref.sce,
  labels = as.character(ref.sce@colData$new_celltype),
  de.method = "wilcox"
)

query_use@meta.data$single_celltype <- pred_celltype[rownames(query_use@meta.data), "labels"]
query_use@meta.data$Celltype2 <- paste0(query_use@meta.data$Celltype, "-", query_use@meta.data$seurat_clusters)

viz <- table(query_use@meta.data$Celltype2, query_use@meta.data$single_celltype)

vizp = viz/rowSums(viz)

print(vizp)

