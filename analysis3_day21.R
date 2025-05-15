library(Seurat) ## version: SeuratObject_4.1.2, Seurat_4.2.0
library(ggplot2)
library(reshape2)
library(harmony) ## version: 0.1.0
library(reticulate) ## for calling Leiden clustering algorithm in python
set.seed(42)

data_path="../h5/" # path to the folder with h5 files, which can be downloaded from GSE250399

# day 21 samples
samples = c("Day21_vHCO","Day21_vHIO","Day21_vHLPO")

# Loading individual samples from h5 files and perform cell-prefiltering
objlist = NULL
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

  obj = FindClusters(obj, algorithm = 4, method="igraph", graph.name =graph.name, resolution=1)

  g = DimPlot(obj, reduction = "umap", pt.size = 0.001, group.by = "seurat_clusters", label=T, label.size = 2) + NoLegend()
  ggsave(file=paste0(i.sample, ".cluster.tiff"), width=4.5, height=4, dpi=300, units="in", compression="lzw")

  saveRDS(obj, file=paste0(i.sample, ".rds"))

  objlist[[i.sample]] = obj

  cat("\n")
}


#################################################################################################
# Integration of Day 21 vHIO/vHCO
#################################################################################################

samples = c("Day21_vHCO","Day21_vHIO")

counts = NULL
cells = NULL

for (sample in samples) {
  i.obj = objlist[[sample]]
  if (is.null(counts)) {
    counts = i.obj@assays$RNA@counts
    cells = i.obj@meta.data
  } else {
    counts = cbind(counts, i.obj@assays$RNA@counts[rownames(counts), ])
    cells = rbind(cells, i.obj@meta.data)
  }
}

obj = CreateSeuratObject(counts=counts, meta.data=cells)
obj = NormalizeData(obj)
obj = SCTransform(obj, vars.to.regress=c("S.Score","G2M.Score", "pMT"))

npcs = 200

obj = RunPCA(obj, npcs=npcs)

hmat = HarmonyMatrix(obj@reductions$pca@cell.embeddings[, 1:npcs],
                     meta_data = obj@meta.data, vars_use = "DataID", do_pca = FALSE)
obj@reductions$harmony = CreateDimReducObject(embeddings = hmat,
                                              assay= DefaultAssay(obj), key="harmony_")

obj = RunUMAP(obj, dims=1:npcs, reduction = "harmony", return.model = TRUE)

obj = FindNeighbors(obj, reduction = "harmony", dims=1:npcs, nn.method = "annoy", annoy.metric = "cosine")
obj = FindClusters(obj, graph.name=graph.name, resolution = 1, algorithm=4, method="igraph")

saveRDS(obj, file="Day21_vHCOvHIO.rds")


#################################################################################################
# combine all day 21 objects
#################################################################################################

meta_use = c("orig.ident","nCount_RNA", "nFeature_RNA", "DataID", "pMT", "S.Score","G2M.Score" ,"Phase","seurat_clusters")

counts = NULL
cells = NULL

counts1 = objlist$Day21_vHLPO@assays$RNA@counts
cells1 = objlist$Day21_vHLPO@meta.data[, meta_use]

counts2 = obj@assays$RNA@counts
cells2 = obj@meta.data[, meta_use]

counts = cbind(counts2, counts1)
cells = rbind(cells2, cells1)

obj = CreateSeuratObject(counts=counts, meta.data=cells)
obj = NormalizeData(obj)

saveRDS(obj, file="Day21.obj.rds")
