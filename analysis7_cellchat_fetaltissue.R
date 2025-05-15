library(Seurat)
library(CellChat)
library(ggplot2)                  
library(patchwork)
library(ggalluvial)
library(igraph)
library(dplyr)

options(stringsAsFactors = FALSE)
set.seed(42)

obj = readRDS(file='../objects/fetaltissues.rds')

print(with(obj@meta.data, table(Tissue, Age)))

obj = subset(obj, Age=="d80" & Tissue %in% c("Colon","Lung-distal","Lung-airway","Jejunum","Ileum","Duodenum"))
obj@meta.data = droplevels(obj@meta.data)

print(with(obj@meta.data, table(Tissue, Age)))

print(with(obj@meta.data, table(celltype, Tissue)))

cat("\nMerge subtype\n")
obj@meta.data$celltype = gsub("Mesenchyme subtype[[:print:]]+", "Mesenchyme", obj@meta.data$celltype)
obj@meta.data$celltype = gsub("Macrophage/monocyte[[:print:]]+", "Macrophage/monocyte", obj@meta.data$celltype)
obj@meta.data$celltype = gsub("T cell/NK cell[[:print:]]+", "T cell/NK cell", obj@meta.data$celltype)

print(with(obj@meta.data, table(celltype, Tissue)))

cat("\nRemove undefined or contaminated cells\n")
obj = subset(obj, celltype != "Contaminant" &  celltype != "Undefined")

print(with(obj@meta.data, table(celltype, Tissue)))

## start cellchat analysis
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)

dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use <- CellChatDB 

samples = names(table(obj@meta.data$Tissue))

cellchat.list = list()

for (i in 1:length(samples)) {
  
  i.sample = i.obj = i.input = i.meta = cellchat = NULL
  
  i.sample = samples[i]
  
  cat("\nTissue:", i.sample,"\n")
  
  i.obj = subset(obj, Tissue == i.sample)
  
  print(i.obj)
  
  i.input = i.obj@assays$RNA@data
  i.meta = droplevels(i.obj@meta.data)
  
  print(sort(table(i.meta$celltype)))
  
  cellchat <- createCellChat(object = i.input, meta = i.meta, group.by = "celltype")
  
  levels(cellchat@idents) 
  groupSize <- as.numeric(table(cellchat@idents)) 
  
  cellchat@DB <- CellChatDB.use
  
  cellchat <- CellChat::subsetData(cellchat)  
  
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
  
  cellchat <- computeCommunProbPathway(cellchat)

  cellchat <- filterCommunication(cellchat, min.cells = 5)
  
  cellchat <- aggregateNet(cellchat)
  
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
  
  saveRDS(cellchat, file=paste0(i.sample, ".cellchat.rds"))
  
  cellchat.list[[i.sample]] = cellchat
  
}

saveRDS(cellchat.list, file="cellchat.list.rds")


lrs = NULL

for (i in 1:length(cellchat.list)) {
  
  i.cellchat = i.df = NULL
  
  i.name = names(cellchat.list)[i]
  i.cellchat = cellchat.list[[i]]
  
  i.df = subsetCommunication(i.cellchat)
  i.df$Dataset = i.name
  
  if (is.null(lrs)) {
    lrs = i.df
  } else {
    lrs = rbind(lrs, i.df)
  }
  
}

write.table(lrs, file="lrs.txt", sep="\t", col.names = T, row.names = F, quote=F)

