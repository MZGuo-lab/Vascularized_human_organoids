library(Seurat)
library(CellChat)
library(ggplot2)
library(patchwork)
library(igraph)
library(dplyr)

options(stringsAsFactors = FALSE)
set.seed(42)

### load the day21 data and exclude low quality cells
obj = readRDS(file="../objects/Day21.obj.rds")

table(obj@meta.data$Celltype, obj@meta.data$DataID)

obj = subset(obj, Celltype != "Low-quality")
obj@meta.data = droplevels(obj@meta.data)

obj@meta.data$celltype = as.character(obj@meta.data$Celltype)

cell_dist = table(obj@meta.data$celltype, obj@meta.data$DataID)

print(cell_dist)

## start cellchat analysis
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)

dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use <- CellChatDB

samples = names(table(obj@meta.data$DataID))

cellchat.list = list()

for (i in 1:length(samples)) {

  i.sample = i.obj = i.input = i.meta = cellchat = NULL

  i.sample = samples[i]

  i.obj = subset(obj, DataID == i.sample)

  i.input = i.obj@assays$RNA@data
  i.meta = droplevels(i.obj@meta.data)

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


# combine the interactions
lrs = NULL

for (i in 1:length(cellchat.list)) {

  i.cellchat = i.df = NULL
  i.name = names(cellchat.list)[i]
  i.cellchat = cellchat.list[[i]]
  i.idents = table(i.cellchat@idents)
  i.df = subsetCommunication(i.cellchat)
  i.df$Dataset = i.name

  if (is.null(lrs)) {
    lrs = i.df
  } else {
    lrs = rbind(lrs, i.df)
  }
}

write.table(lrs, file="lrs.txt", sep="\t", col.names = T, row.names = F, quote=F)
