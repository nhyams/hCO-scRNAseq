library(CellChat)
library(patchwork)
library(Seurat)
library(future)
options(future.globals.maxSize = 10e9)

threeandfour <- readRDS("C:/Users/nateh/OneDrive/Documents/Mei Lab/RNAseq/scRNAseq/Integrated Datasets/3 + 4 cell type organoids/Isogenic_Organoids_simple.rds")

DimPlot(threeandfour, reduction = 'umap')

CellChatinput <- threeandfour[["RNA"]]$data
labels <- Idents(threeandfour)
meta <- data.frame(labels = labels, row.names = names(labels))

CellChatobject <- createCellChat(object = threeandfour, group.by = "ident", assay = "RNA")

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use <- subsetDB(CellChatDB, search = 'Cell-Cell Contact', key = 'annotation')

CellChatobject@DB <- CellChatDB.use

#Preprocessing
CellChatobject <- subsetData(CellChatobject)
future::plan('multisession', workers = 4)
CellChatobject <- identifyOverExpressedGenes(CellChatobject)
CellChatobject <- identifyOverExpressedInteractions(CellChatobject)

ptm = Sys.time()
CellChatobject <- computeCommunProb(CellChatobject, type = 'triMean')

CellChatobject <- filterCommunication(CellChatobject, min.cells = 10)

df.net <- subsetCommunication(CellChatobject, slot.name = 'netP')

CellChatobject <- computeCommunProbPathway(CellChatobject)
CellChatobject <- aggregateNet(CellChatobject)

groupSize <- as.numeric(table(CellChatobject@idents))
par(mfrow = c(1,2), xpd = T)
pathways.show <- c("NOTCH")
par(mfrow=c(1,1))
plot <- netVisual_aggregate(CellChatobject, 
                    signaling = pathways.show)

setwd("C:/Users/nateh/OneDrive/Documents/Mei Lab/RNAseq/scRNAseq/Integrated Datasets/3 + 4 cell type organoids/CellChat/CellChat Simplified Annotations/")

png("notch.png", width = 6, height = 8, units = "in", res = 600)
netVisual_aggregate(CellChatobject, signaling = pathways.show)
dev.off()


#Visualize
circle_plot <- netVisual_circle(CellChatobject@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Cardiac Organoid Interactions")
circle_plot <- grid.echo()
grid.draw(circle_plot)
netVisual_bubble(CellChatobject, signaling = pathways.show)
netVisual_embedding(CellChatobject, type = 'functional', layout = 'fr', signaling = pathways.show)
pathways.show <- c("AKT")
par(mfrow=c(1,1))
netVisual_aggregate(CellChatobject, signaling = pathways.show, target.cells = c('Pericytes'))

netVisual_diffInteraction(CellChatobject,
                          comparison = c(1,2),
                          measure = 'weight',
                          weight.scale = F)


