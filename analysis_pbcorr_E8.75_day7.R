library(Seurat)
library(ggplot2)
library(reshape2)

# pseudo-bulk correlation analysis of mouse E8.75 anterior and posterior gut, EC, mes cells and day 7 Epi, Ec, and Mesen cell clusters

# Load data
ref <- readRDS("../objects/GSE123046_E8.75_AP.Human.rds")
ref <- subset(ref, CellType %in% c("Endothelial", "Gut tube", "Mes"))
ref@meta.data <- droplevels(ref@meta.data)

afg_b1 <- readRDS("../objects/Day7_vAFG_B1.rds")
mhg_b3 <- readRDS("../objects/Day7_vMHG_B3.rds")

# Subset for relevant clusters
afg_b1 <- subset(afg_b1, seurat_clusters %in% 1:7)
afg_b1@meta.data <- droplevels(afg_b1@meta.data)
mhg_b3@meta.data$type_sample <- paste0(mhg_b3@meta.data$DataID, "_C", mhg_b3@meta.data$seurat_clusters)
afg_b1@meta.data$type_sample <- paste0(afg_b1@meta.data$DataID, "_C", afg_b1@meta.data$seurat_clusters)
ap <- merge(afg_b1, mhg_b3)

# Pseudobulk
ref$group <- gsub(" ", "", paste0(ref$Sample, "_", ref$CellType))
ref$group <- gsub("Endothelial", "EC", ref$group)
ref$group <- gsub("Guttube", "Gut", ref$group)

ap$group <- ap$type_sample
DefaultAssay(ref) <- "RNA"
DefaultAssay(ap) <- "RNA"

ref <- SetIdent(ref, value = "group")
ref.avg <- AverageExpression(ref, assay = "RNA", return.seurat = TRUE)
ref.avg <- FindVariableFeatures(ref.avg, nfeatures = 2000)

ap <- SetIdent(ap, value = "group")
ap.avg <- AverageExpression(ap, assay = "RNA", return.seurat = TRUE)
ap.avg <- FindVariableFeatures(ap.avg, nfeatures = 2000)

# Shared HVGs
hvg.common <- intersect(ref.avg@assays$RNA@var.features, ap.avg@assays$RNA@var.features)
hvg.common <- intersect(hvg.common, rownames(ap.avg@assays$RNA@data))
hvg.common <- intersect(hvg.common, rownames(ref.avg@assays$RNA@data))

# Scale
ref.avg <- ScaleData(ref.avg, features = hvg.common)
ap.avg <- ScaleData(ap.avg, features = hvg.common)

# Combine and correlate
ref.data <- ref.avg@assays$RNA@scale.data
ap.data <- ap.avg@assays$RNA@scale.data
combined <- cbind(ref.data[hvg.common, ], ap.data[hvg.common, ])

cor_mat = cor(combined)



