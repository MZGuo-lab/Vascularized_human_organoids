library(Seurat)
library(SingleR)
library(SingleCellExperiment)
library(ggplot2)

ref = readRDS(file="../objects/fetaltissues_mesenchyme_atlas.rds")
ref = subset(ref, Age == "d80" & Tissue %in% c("Ileum","Jejunum","Lung-distal"))
# exclude cells with annotations of proliferation
ref = subset(ref, Cell_type %in% c("Proliferative","Lung proliferative","SFRP1+ proli."), invert=T)
ref@meta.data = droplevels(ref@meta.data)

DefaultAssay(ref) = "RNA"
ref = NormalizeData(ref)
ref.sce = as.SingleCellExperiment(ref)

query = readRDS(file = "../objects/Day21.obj.rds")
query = subset(query, Celltype %in% c("FB","pFB","Pericyte","SMC"))
query@meta.data = droplevels(query@meta.data)

DefaultAssay(query) = "RNA"
query = NormalizeData(query)
query.sce = as.SingleCellExperiment(query)

pred <- SingleR(test=query.sce, ref=ref.sce, 
                 labels=as.character(ref.sce@colData$Organ), 
                 de.method="wilcox")

all(rownames(query@meta.data) == rownames(pred))
query@meta.data$pruned.labels = pred[rownames(query@meta.data), "pruned.labels"]
table(query@meta.data$pruned.labels, query@meta.data$DataID)

saveRDS(query, file="Day21_mesen.singleR_pred.rds")

