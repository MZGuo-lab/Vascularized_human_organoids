library(Seurat)
library(harmony)
library(ggplot2)

obj = readRDS(file="../objects/Day21.obj.rds")

obj1 = subset(obj, DataID %in% c("Day21_vHCO", "Day21_vHIO") & seurat_clusters %in% c(9,13))
obj2 = subset(obj, DataID %in% c("Day21_vHLPO") & seurat_clusters %in% c(8))

genes = rownames(obj1@assays$RNA@counts)
counts = cbind(obj1@assays$RNA@counts, obj2@assays$RNA@counts[genes, ])

cells = rbind(obj1@meta.data, obj2@meta.data)

all(rownames(cells) == colnames(counts))

obj = CreateSeuratObject(counts, meta.data=cells)
obj = NormalizeData(obj)

obj = SCTransform(obj, vars.to.regress = c("S.Score","G2M.Score", "pMT"))
obj = RunPCA(obj, npcs=200)
ElbowPlot(obj, ndims=200)

npcs = 50
batch_var = "DataID"

hmat = HarmonyMatrix(obj@reductions$pca@cell.embeddings[, 1:npcs], meta_data = obj@meta.data, 
                     vars_use = batch_var, do_pca = FALSE)
obj@reductions$harmony = CreateDimReducObject(embeddings = hmat, 
                                              assay= DefaultAssay(obj), key="harmony_")

obj = RunUMAP(obj, dims=1:npcs, reduction = "harmony", min.dist=0.3, 
              seed.use = 42L, return.model = TRUE)

saveRDS(obj, file=paste0("Day21_endo.rds"))


