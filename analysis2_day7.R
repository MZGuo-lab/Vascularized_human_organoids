library(Seurat) ## version: SeuratObject_4.1.2, Seurat_4.2.0
library(ggplot2)
library(reshape2)
library(reticulate) ## for calling Leiden clustering algorithm in python
set.seed(42)

data_path="../h5/" # path to the folder with h5 files, which can be downloaded from GSE250399

# day 7 samples
samples = c("Day7_vAFG_B0","Day7_vAFG_B1",
            "Day7_vMHG_B1","Day7_vMHG_B2","Day7_vMHG_B3")

# Loading individual samples from h5 files and perform cell-prefiltering
for (i in 1:length(samples)) {

  counts = obj = meta = i.sample = NULL
  npcs.use = 50

  i.sample = samples[i]

  cat(i.sample, "\n")

  counts = Read10X_h5(filename = paste0(data_path, i.sample, "_filtered_feature_bc_matrix.h5"))

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

  obj = FindClusters(obj, algorithm = 4, method="igraph", graph.name =graph.name)

  g = DimPlot(obj, reduction = "umap", pt.size = 0.001, group.by = "seurat_clusters", label=T, label.size = 2) + NoLegend()
  ggsave(file=paste0(i.sample, ".cluster.tiff"), width=4.5, height=4, dpi=300, units="in", compression="lzw")

  saveRDS(obj, file=paste0(i.sample, ".rds"))

  cat("\n")
}
