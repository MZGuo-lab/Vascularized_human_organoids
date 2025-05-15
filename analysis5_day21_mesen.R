library(Seurat)
library(harmony)
library(ggplot2)

obj = readRDS(file="../objects/Day21.obj.rds")

obj1 = subset(obj, DataID %in% c("Day21_vHCO", "Day21_vHIO")) 
obj1 = subset(obj1, seurat_clusters %in% c(6,9,13,14), invert=T)
obj2 = subset(obj, DataID %in% c("Day21_vHLPO") & seurat_clusters %in% c(1,2,4,5,6,9,11,18,14))

genes = rownames(obj1@assays$RNA@counts)
counts = cbind(obj1@assays$RNA@counts, 
               obj2@assays$RNA@counts[genes, ])

cells = rbind(obj1@meta.data, obj2@meta.data)

all(rownames(cells) == colnames(counts))

obj = CreateSeuratObject(counts, meta.data=cells)
obj = NormalizeData(obj)

npcs=200

obj = SCTransform(obj, vars.to.regress = c("S.Score","G2M.Score", "pMT"))
obj = RunPCA(obj, npcs=npcs)

hmat = HarmonyMatrix(obj@reductions$pca@cell.embeddings, meta_data = obj@meta.data, 
                     vars_use = "DataID", do_pca = FALSE)
obj@reductions$harmony = CreateDimReducObject(embeddings = hmat, 
                                              assay= DefaultAssay(obj), key="harmony_")

obj = RunUMAP(obj, dims=1:npcs, reduction = "harmony", return.model = TRUE)

saveRDS(obj, file=paste0("Day21_mesen.rds"))
