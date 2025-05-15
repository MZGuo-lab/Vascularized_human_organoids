library(Seurat)
library(ggplot2)
library(reshape2)
library(harmony)
set.seed(42)

# day 3 samples
samples = c("Day3_B0","Day3_B1","Day3_B2","Day3_B3")

# Loading individual samples from h5 files and perform cell-prefiltering
counts_all = cells_all = NULL
for (i in 1:length(samples)) {

  counts = obj = meta = i.sample = NULL
  npcs.use = 50

  i.sample = samples[i]

  cat(i.sample, "\n")

  counts = Read10X_h5(filename = paste0(i.sample, "_filtered_feature_bc_matrix.h5"))

  print(dim(counts))

  meta = data.frame(DataID=i.sample, Barcode=colnames(counts))
  rownames(meta) = paste0(i.sample, "_", colnames(counts))
  colnames(counts) = paste0(i.sample, "_", colnames(counts))

  obj = CreateSeuratObject(counts=counts, project = i.sample, meta.data = meta)

  print(obj)

  mt.genes = grep("^MT-", rownames(counts), value=T)
  print(mt.genes)

  obj@meta.data$pMT = PercentageFeatureSet(obj, features=mt.genes)[, 1]

  obj@meta.data$QC = "Cell"
  obj@meta.data[rownames(subset(obj@meta.data, nFeature_RNA>=500 & pMT< 25 & nCount_RNA < 100000)), "QC"] = "Selected"

  cat("\nSelection\n")
  print(table(obj@meta.data$QC))

  obj = subset(obj, QC=="Selected")

  print(obj)

  obj = NormalizeData(obj)

  obj = CellCycleScoring(obj, s.features = Seurat::cc.genes.updated.2019$s.genes,
                         g2m.features = Seurat::cc.genes.updated.2019$g2m.genes)

  obj = SCTransform(obj, vars.to.regress=c("S.Score","G2M.Score","pMT"))

  obj = RunPCA(obj, npcs=npcs.use)
  obj = RunUMAP(obj, reduction = "pca", dims=1:npcs.use)

  DimPlot(obj, reduction = "umap", pt.size = 0.1)

  obj = FindNeighbors(obj, dims=1:npcs.use)

  graph.name = names(obj@graphs)[1]
  obj = FindClusters(obj, algorithm = 4, method="igraph", graph.name = graph.name)

  DimPlot(obj, reduction = "umap", pt.size = 0.001, group.by = "seurat_clusters", label=T, label.size = 2)

  saveRDS(obj, file=paste0(i.sample, ".rds"))

  if (is.null(counts_all)) {
    counts_all = obj@assays$RNA@counts
  } else {
    counts_all = cbind(counts_all, obj@assays$RNA@counts[rownames(counts_all), ])
  }

  if (is.null(cells_all)) {
    cells_all = obj@meta.data
  } else {
    cells_all = rbind(cells_all, obj@meta.data)
  }

  cat("\n")
}

saveRDS(counts_all, file="day3_counts.rds")
saveRDS(cells_all, file="day3_cells.rds")


# integration using harmony
npcs = 200
obj = CreateSeuratObject(counts=counts_all, project = "Day3", meta.data = cells_all)
obj = NormalizeData(obj)

obj = SCTransform(obj, vars.to.regress=c("S.Score","G2M.Score","pMT"))

obj = RunPCA(obj, npcs=npcs)

hmat = HarmonyMatrix(obj@reductions$pca@cell.embeddings[, 1:npcs],
                     meta_data = obj@meta.data, vars_use = "DataID",
                     do_pca = FALSE)
obj@reductions$harmony = CreateDimReducObject(embeddings = hmat, assay= DefaultAssay(obj), key="harmony_")

obj = RunUMAP(obj, dims=1:npcs, reduction = "harmony", return.model = TRUE)

DimPlot(obj, group.by="DataID")

obj = FindNeighbors(obj, reduction = "harmony", dims=1:npcs, nn.method = "annoy", annoy.metric = "cosine")

graph.name = names(obj@graphs)[1]
obj = FindClusters(obj, graph.name=graph.name, resolution = 0.8, algorithm=4, method="igraph")

DimPlot(obj, reduction = "umap", pt.size = 0.001,
        group.by = c("seurat_clusters","DataID"), label=T, label.size = 2) + NoLegend()

obj$Celltype = "Def. endoderm"
obj$Celltype[which(obj$seurat_clusters %in% c(2,3,7))] = "Meso-endoderm progenitor"
obj$Celltype[which(obj$seurat_clusters %in% c(4,13))] = "Lateral plate mesoderm"
obj$Celltype[which(obj$seurat_clusters %in% c(8))] = "Cardiac mesendoderm"
obj$Celltype[which(obj$seurat_clusters %in% c(11))] = "Low quality"

saveRDS(obj, file="Day3_all.rds")
