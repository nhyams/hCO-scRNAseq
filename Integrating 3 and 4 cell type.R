library(Seurat)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(gridExtra)
library(harmony)
library(sctransform)
library(RPresto)
library(presto)
library(Nebulosa)
library(SeuratDisk)
library(DoubletFinder)


L1.data <- Read10X(data.dir = "C:/Users/nateh/OneDrive/Documents/Mei Lab/RNAseq/scRNAseq/L1/L1 Data")
L1 <- CreateSeuratObject(counts = L1.data, 
                         project = "L1", 
                         min.cells = 3, 
                         min.features = 200)


L2.data <- Read10X(data.dir = "C:/Users/nateh/OneDrive/Documents/Mei Lab/RNAseq/scRNAseq/L2/L2 Data")
L2 <- CreateSeuratObject(counts = L2.data, 
                         project = "L2", 
                         min.cells = 3, 
                         min.features = 200)

L3.data <- Read10X(data.dir = "C:/Users/nateh/OneDrive/Documents/Mei Lab/RNAseq/scRNAseq/L3/L3 Data")
L3 <- CreateSeuratObject(counts = L3.data, 
                         project = "L3", 
                         min.cells = 3, 
                         min.features = 200)

L4.data <- Read10X(data.dir = "C:/Users/nateh/OneDrive/Documents/Mei Lab/RNAseq/scRNAseq/L4/Reanalyzed Data/L4 Data")
L4 <- CreateSeuratObject(counts = L4.data, 
                         project = "L4", 
                         min.cells = 3, 
                         min.features = 200)

L1[["percent.mt"]] <- PercentageFeatureSet(L1, pattern = "^MT-")
VlnPlot(L1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
L2[["percent.mt"]] <- PercentageFeatureSet(L2, pattern = "^MT-")
VlnPlot(L2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
L3[["percent.mt"]] <- PercentageFeatureSet(L3, pattern = "^MT-")
VlnPlot(L3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
L4[["percent.mt"]] <- PercentageFeatureSet(L4, pattern = "^MT-")
VlnPlot(L4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#### Merge datasets together
merged_datasets <- merge(L1, y = list(L2, L3, L4),
                         add.cell.ids = c("3 cell type", "3 cell type", "4 cell type", "4 cell type"))
merged_datasets$sample <- rownames(merged_datasets@meta.data)
#Split sample column
merged_datasets@meta.data <- separate(merged_datasets@meta.data, col = "sample", 
                                      into = c("Organoid Type", "Barcode", "batch"), sep = "_")
#Run Seurat Pipeline
merged_datasets <- NormalizeData(merged_datasets)
merged_datasets <- FindVariableFeatures(merged_datasets)
merged_datasets <- ScaleData(merged_datasets)
merged_datasets <- RunPCA(merged_datasets)
merged_datasets <- RunUMAP(merged_datasets, dims = 1:20)
merged_datasets <- RunTSNE(merged_datasets)
DimPlot(merged_datasets, reduction = 'umap')
DimPlot(merged_datasets, reduction = 'tsne')
DimPlot(merged_datasets, reduction = 'pca')

## Subset out low quality cells
merged_datasets <- subset(merged_datasets,
                          nFeature_RNA > 200 &
                            nFeature_RNA < 8000 &
                            nCount_RNA > 800 &
                            percent.mt < 25)
DimPlot(merged_datasets, reduction = 'umap')


### Harmony integrates our libraries together
threeandfour <- RunHarmony(merged_datasets, group.by.vars = 'batch', 
                           reduction.use = 'pca',
                           reduction.save = 'harmony')
threeandfour <- RunUMAP(threeandfour, reduction = 'harmony', dims = 1:30)
threeandfour <- FindNeighbors(threeandfour, reduction = "harmony")
threeandfour <- FindClusters(threeandfour, resolution = 0.3)
DimPlot(threeandfour, reduction = "umap", label = T)
DimPlot(threeandfour, reduction = "umap", group.by = "batch")
DimPlot(threeandfour, reduction = "umap", group.by = "Organoid Type")

p <- DimPlot(threeandfour, reduction = "umap", split.by = "Organoid Type", ) + labs(
  title = '3 vs 4 Cell Type Isogenic Cardiac Organoids',
  x = 'UMAP 1',
  y = 'UMAP 2'
) + 
  theme(plot.title = element_text(hjust = 0.5, size = 16))

DimPlot(threeandfour, features("TNNT2"))
FeaturePlot(threeandfour, features = c("PDGFRB", "RGS5"), blend = T)
FeaturePlot(threeandfour, features = c("TNNT2", "NPPA"), blend = T, min.cutoff = 'q10')

#FindMarkers
threeandfour <- JoinLayers(threeandfour)
allmarkers <- FindAllMarkers(threeandfour,
               logfc.threshold = 0.25,
               min.pct = 0.1,
               only.pos = T)

### Using Nebulosa to map out specific features of interest
plot_density(threeandfour, features = c("HEY2"), reduction = "umap")
#Atrial Markers
plot_density(threeandfour, features = c("NPPA"), reduction = "umap")
#Fibroblasts
plot_density(threeandfour, features = c("PDGFRA"), reduction = "umap")
#Endothelial Cells
plot_density(threeandfour, features = c("PECAM1"), reduction = "umap")

plot_density(threeandfour, features = c("HIF1A", "HCN4", "KLF4"), reduction = "umap")
#Proliferating CMs
plot_density(threeandfour, features = c("CDK1"), reduction = "umap")
#Cardiac Conduction System
plot_density(threeandfour, features = c("SHOX2", "KCNJ2", "TBX5"), reduction = "umap")
#Metabolically Active CMs
plot_density(threeandfour, features = c("PHGDH", "CHAC1", "PSAT1"), reduction = "umap")
#Cluster 4
plot_density(threeandfour, features = c("RYR2", "CACNA1C", "CTNNA3", "DLG2"), reduction = "umap")
#SIRPA characterization
plot_density(integrated_data, features = c("SIRPA", "CD200", "COL1A1", "TNNI1"), reduction = "umap")
# Pericytes
plot_density(threeandfour, features = c("PDGFRB", "RGS5", "HIGD1B"), reduction = "umap")

FeaturePlot(threeandfour, features = c("HEY2", "MYL2", "IRX3", "GJA5"))
# Atrial Markers 2.0
FeaturePlot(threeandfour, features = c("NPPA"))
# 1st heart field vs 2nd heart field
FeaturePlot(threeandfour, features = c("NKX2-5", "TBX5", "HCN4", "ISL1"), reduction = "umap")

plot_density(threeandfour, features = c("ACTA2"), reduction = "umap")


# Violin Plots
VlnPlot(threeandfour, c(features = "MYH7", "MYL2", "TNNT2", "HEY2", "IRX4"))
VlnPlot(threeandfour, c(features = "MYL4", "NPPA"))
VlnPlot(integrated_data, c(features = "COL1A1", "COL3A1", "FAP", "VIM", "PDGFRA"))
VlnPlot(integrated_data, c(features = "PECAM1", "FLT1", "CDH5"))
VlnPlot(threeandfour, c(Features = "CDK1"))
VlnPlot(threeandfour, c(Features = "SHOX2", "KCNJ3", "TBX5"))
VlnPlot(integrated_data, c(Features = "PHGDH", "CHAC1", "PSAT1"))
VlnPlot(integrated_data, c(Features = "RYR2", "CACNA1C", "CTNNA3"))
VlnPlot(integrated_data, c(Features = "PHGDH", "CACNA1C", "CTNNA3"))
VlnPlot(threeandfour, c(Features = "PDGFRB", "RGS5"))
VlnPlot(integrated_data, c(Features = "SIRPA", "CD200", "COL1A1", "TNNI1"))
VlnPlot(threeandfour, c(Features = "ESM1", "CD93"))

threeandfour <- JoinLayers(threeandfour)
## based on iterative findings of our samples
new_cluster_ids <- c(
  "0" = "Ventricular-like CMs",
  "1" = "Ventricular-like CMs",
  "2" = "Atrial-like CMs",
  "3" = "Ventricular CMs",
  "4" = "Cardiac Pacemaker Cells",
  "5" = "Fibroblasts 1",
  "6" = "Ventricular CMs",
  "7" = "Pericytes",
  "8" = "Endothelial Cells",
  "9" = "Cardiac Fibroblasts 3",
  "10" = "Proliferative CMs",
  "11" = "Cardiac Pacemaker Cells",
  '12' = 'Cardiac Fibroblasts 2',
  '13' = 'Cardiac Fibroblasts 1',
  '14' = 'Fibroblasts 2')
### Annotations for cell types
threeandfour <- RenameIdents(threeandfour, new_cluster_ids)
plot <- DimPlot(threeandfour, reduction = "umap") + ggtitle("Isogenic Human Cardiac Organoids") + 
  theme(element_text(family = "Arial", size = 12),
  plot.title = element_text(size = 15, face = "bold"),
  axis.title = element_text(size = 12, face = "bold"),  
  axis.text = element_text(size = 12),                       
  legend.title = element_text(size = 12, face = "bold"),    
  legend.text = element_text(size = 12))
print(plot)
Idents(threeandfour)
LabelClusters(plot, id = "ident", fontface = "bold")

DimPlot(threeandfour, reduction = 'umap', split.by = 'Organoid Type')
DimPlot(threeandfour, reduction = 'umap', group.by = 'Organoid Type') +  theme(element_text(family = "Arial", size = 12),
                                                                               plot.title = element_text(size = 16, face = "bold"),
                                                                               axis.title = element_text(size = 12, face = "bold"),  
                                                                               axis.text = element_text(size = 12),                       
                                                                               legend.title = element_text(size = 12, face = "bold"),    
                                                                               legend.text = element_text(size = 12))
threeandfour <- JoinLayers(threeandfour)


saveRDS(threeandfour, file = 'C:/Users/nateh/OneDrive/Documents/Mei Lab/RNAseq/scRNAseq/Integrated Datasets/3 + 4 cell type organoids/With L4 Included/Isogenic_Organoids.rds')

simple_cluster_ids <- c(
  "0" = "Ventricular-like CMs",
  "1" = "Ventricular-like CMs",
  "2" = "Ventricular-like CMs",
  "3" = "Ventricular-like CMs",
  "4" = "Atrial-like CMs",
  "5" = "Ventricular-like CMs",
  "6" = "Ventricular-like CMs",
  "7" = "Cardiac Pacemaker Cells",
  "8" = "Cardiac Pacemaker Cells",
  "9" = "Pericytes",
  "10" = "Endothelial Cells",
  "11" = "Cardiac Fibroblasts",
  '12' = 'Cardiac Fibroblasts',
  '13' = 'Cardiac Fibroblasts')

organoids <- readRDS("C:/Users/nateh/OneDrive/Documents/Mei Lab/RNAseq/scRNAseq/Integrated Datasets/3 + 4 cell type organoids/integrated1.rds")
Idents(organoids) <- organoids$`seurat_clusters`
organoids <- RenameIdents(organoids, simple_cluster_ids)
DimPlot(organoids)

## How 3 cell vs 4 cell impacts each cell cluster
#### Adding in a new column of data to compare condition impacts by cell cluster
threeandfour <- AddMetaData(object = threeandfour, metadata = Idents(threeandfour), col.name = "Seurat Annotations")
threeandfour$celltype.cnd <- paste0(threeandfour@meta.data$`Seurat Annotations`,'_',threeandfour@meta.data$`Organoid Type`)
Idents(threeandfour) <- threeandfour$`celltype.cnd`
Idents(threeandfour) <- threeandfour$`Seurat Annotations`
DimPlot(threeandfour, reduction = 'umap', label = T)
Vent1markers <- FindConservedMarkers(threeandfour, ident.1 = "Ventricular CMs", grouping.var = 'Organoid Type')
Vent2markers <- FindConservedMarkers(threeandfour, ident.1 = "Ventricular CMs 2", grouping.var = 'Organoid Type')
Vent3markers <- FindConservedMarkers(threeandfour, ident.1 = "Ventricular CMs 3", grouping.var = 'Organoid Type')
PM1markers <- FindConservedMarkers(threeandfour, ident.1 = "Cardiac Pacemaker Cells 1", grouping.var = 'Organoid Type')
PM2markers <- FindConservedMarkers(threeandfour, ident.1 = "Cardiac Pacemaker Cells 2", grouping.var = 'Organoid Type')
ICM1markers <- FindConservedMarkers(threeandfour, ident.1 = "Immature CMs 1", grouping.var = 'Organoid Type')
ICM2markers <- FindConservedMarkers(threeandfour, ident.1 = "Immature CMs 2", grouping.var = 'Organoid Type')

ECmarkers <- FindConservedMarkers(threeandfour, ident.1 = "Endothelial Cells", grouping.var = 'Organoid Type')
cFBmarkers <- FindConservedMarkers(threeandfour, ident.1 = "Cardiac Fibroblasts", grouping.var = 'Organoid Type')
Atrialmarkers <- FindConservedMarkers(threeandfour, ident.1 = "Atrial CMs", grouping.var = 'Organoid Type')

Idents(threeandfour) <- threeandfour$celltype.cnd
Vent1cnd <- FindMarkers(threeandfour, ident.1 = "Ventricular CMs_3 cell type", ident.2 = "Ventricular CMs_4 cell type")
cfbcnd <- FindMarkers(threeandfour, ident.1 = "Cardiac Fibroblasts_3 cell type", ident.2 = "Cardiac Fibroblasts_4 cell type")
### Cluster isolation of significantly expressed genes- crucial to our understanding- Ventricular 1
VentDEGs <- FindMarkers(threeandfour, ident.1 = "Ventricular CMs_3 cell type", ident.2 = "Ventricular CMs_4 cell type", logfc.threshold = 0.25, min.pct = 0.1, test.use = 'wilcox')
SigVent1DEGs <- VentDEGs[VentDEGs$p_val_adj < 0.01 & abs(VentDEGs$avg_log2FC) > 0.5,]
SigVent1DEGs$significance <- "Not Significant"
SigVent1DEGs$significance[SigVent1DEGs$p_val_adj < 0.001 & SigVent1DEGs$avg_log2FC < -0.5] <- "Downregulated"
SigVent1DEGs$significance[SigVent1DEGs$p_val_adj < 0.001 & SigVent1DEGs$avg_log2FC > 0.5] <- "Upregulated"
SigVent1DEGs <- rownames_to_column(SigVent1DEGs, var = 'gene')
write.csv(SigVent1DEGs, file = "C:/Users/nateh/OneDrive/Documents/Mei Lab/RNAseq/scRNAseq/Integrated Datasets/3 + 4 cell type organoids/Vent1clustersigDEGs2.csv")

### Cluster isolation of significantly expressed genes- crucial to our understanding- Ventricular 2
Vent2DEGs <- FindMarkers(threeandfour, ident.1 = "Ventricular CMs 2_3 cell type", ident.2 = "Ventricular CMs 2_4 cell type", logfc.threshold = 0.25, min.pct = 0.1, test.use = 'wilcox')
SigVent2DEGs <- Vent2DEGs[Vent2DEGs$p_val_adj < 0.01 & abs(Vent2DEGs$avg_log2FC) > 0.5,]
SigVent2DEGs$significance <- "Not Significant"
SigVent2DEGs$significance[SigVent2DEGs$p_val_adj < 0.001 & SigVent2DEGs$avg_log2FC < -0.5] <- "Downregulated"
SigVent2DEGs$significance[SigVent2DEGs$p_val_adj < 0.001 & SigVent2DEGs$avg_log2FC > 0.5] <- "Upregulated"
SigVent2DEGs <- rownames_to_column(SigVent2DEGs, var = 'gene')
write.csv(SigVent2DEGs, file = "C:/Users/nateh/OneDrive/Documents/Mei Lab/RNAseq/scRNAseq/Integrated Datasets/3 + 4 cell type organoids/Vent2clustersigDEGs.csv")

### Cluster isolation of significantly expressed genes- crucial to our understanding- Atrial
AtrialDEGs <- FindMarkers(threeandfour, ident.1 = "Atrial CMs_3 cell type", ident.2 = "Atrial CMs_4 cell type", logfc.threshold = 0.25, min.pct = 0.1, test.use = 'wilcox')
AtrialDEGs <- AtrialDEGs[AtrialDEGs$p_val_adj < 0.01 & abs(AtrialDEGs$avg_log2FC) > 0.5,]
AtrialDEGs$significance <- "Not Significant"
AtrialDEGs$significance[AtrialDEGs$p_val_adj < 0.001 & AtrialDEGs$avg_log2FC < -0.5] <- "Downregulated"
AtrialDEGs$significance[AtrialDEGs$p_val_adj < 0.001 & AtrialDEGs$avg_log2FC > 0.5] <- "Upregulated"
AtrialDEGs <- rownames_to_column(AtrialDEGs, var = 'gene')
write.csv(AtrialDEGs, file = "C:/Users/nateh/OneDrive/Documents/Mei Lab/RNAseq/scRNAseq/Integrated Datasets/3 + 4 cell type organoids/AtrialclustersigDEGs.csv")

### Cluster isolation of significantly expressed genes- crucial to our understanding- cFBs
cFBDEGs <- FindMarkers(threeandfour, ident.1 = "Cardiac Fibroblasts_3 cell type", ident.2 = "Cardiac Fibroblasts_4 cell type", logfc.threshold = 0.25, min.pct = 0.1, test.use = 'wilcox')
cFBDEGs <- cFBDEGs[cFBDEGs$p_val_adj < 0.01 & abs(cFBDEGs$avg_log2FC) > 0.5,]
cFBDEGs$significance <- "Not Significant"
cFBDEGs$significance[cFBDEGs$p_val_adj < 0.001 & cFBDEGs$avg_log2FC < -0.5] <- "Downregulated"
cFBDEGs$significance[cFBDEGs$p_val_adj < 0.001 & cFBDEGs$avg_log2FC > 0.5] <- "Upregulated"
cFBDEGs <- rownames_to_column(cFBDEGs, var = 'gene')
write.csv(cFBDEGs, file = "C:/Users/nateh/OneDrive/Documents/Mei Lab/RNAseq/scRNAseq/Integrated Datasets/3 + 4 cell type organoids/cFBclustersigDEGs.csv")

### Cluster isolation of significantly expressed genes- crucial to our understanding- ECs
ECDEGs <- FindMarkers(threeandfour, ident.1 = "Endothelial Cells_3 cell type", ident.2 = "Endothelial Cells_4 cell type", logfc.threshold = 0.1, min.pct = 0.1, test.use = 'wilcox')
ECDEGs <- ECDEGs[ECDEGs$p_val_adj < 0.01 & abs(ECDEGs$avg_log2FC) > 0.5,]
ECDEGs$significance <- "Not Significant"
ECDEGs$significance[ECDEGs$p_val_adj < 0.01 & ECDEGs$avg_log2FC < -0.5] <- "Downregulated"
ECDEGs$significance[ECDEGs$p_val_adj < 0.01 & ECDEGs$avg_log2FC > 0.5] <- "Upregulated"
ECDEGs <- rownames_to_column(ECDEGs, var = 'gene')
write.csv(ECDEGs, file = "C:/Users/nateh/OneDrive/Documents/Mei Lab/RNAseq/scRNAseq/Integrated Datasets/3 + 4 cell type organoids/ECclustersigDEGs2.csv")

### Cluster isolation of significantly expressed genes- crucial to our understanding- EC 2
EC2DEGs <- FindMarkers(threeandfour, ident.1 = "Endothelial Cells 2_3 cell type", ident.2 = "Endothelial Cells 2_4 cell type", logfc.threshold = 0.25, min.pct = 0.1, test.use = 'wilcox')
EC2DEGs <- EC2DEGs[ECDEGs$p_val_adj < 0.01 & abs(EC2DEGs$avg_log2FC) > 0.5,]
EC2DEGs$significance <- "Not Significant"
EC2DEGs$significance[EC2DEGs$p_val_adj < 0.001 & EC2DEGs$avg_log2FC < -0.5] <- "Downregulated"
EC2DEGs$significance[EC2DEGs$p_val_adj < 0.001 & EC2DEGs$avg_log2FC > 0.5] <- "Upregulated"
ECDEGs <- rownames_to_column(ECDEGs, var = 'gene')
write.csv(ECDEGs, file = "C:/Users/nateh/OneDrive/Documents/Mei Lab/RNAseq/scRNAseq/Integrated Datasets/3 + 4 cell type organoids/ECclustersigDEGs.csv")

### Cluster isolation of significantly expressed genes- crucial to our understanding- Ventricular CMs 3
VCM3DEGs <- FindMarkers(threeandfour, ident.1 = "Ventricular CMs 3_3 cell type", ident.2 = "Ventricular CMs 3_4 cell type", logfc.threshold = 0.25, min.pct = 0.1, test.use = 'wilcox')
VCM3DEGs <- VCM3DEGs[VCM3DEGs$p_val_adj < 0.01 & abs(VCM3DEGs$avg_log2FC) > 0.5,]
VCM3DEGs$significance <- "Not Significant"
VCM3DEGs$significance[VCM3DEGs$p_val_adj < 0.001 & VCM3DEGs$avg_log2FC < -0.5] <- "Downregulated"
VCM3DEGs$significance[VCM3DEGs$p_val_adj < 0.001 & VCM3DEGs$avg_log2FC > 0.5] <- "Upregulated"
VCM3DEGs <- rownames_to_column(VCM3DEGs, var = 'gene')
write.csv(VCM3DEGs, file = "C:/Users/nateh/OneDrive/Documents/Mei Lab/RNAseq/scRNAseq/Integrated Datasets/3 + 4 cell type organoids/VCM3clustersigDEGs.csv")

### Cluster isolation of significantly expressed genes- crucial to our understanding- Pacemaker 1
PM1DEGs <- FindMarkers(threeandfour, ident.1 = "Cardiac Pacemaker Cells 1_3 cell type", ident.2 = "Cardiac Pacemaker Cells 1_4 cell type", logfc.threshold = 0.25, min.pct = 0.1, test.use = 'wilcox')
PM1DEGs <- PM1DEGs[PM1DEGs$p_val_adj < 0.01 & abs(PM1DEGs$avg_log2FC) > 0.5,]
PM1DEGs$significance <- "Not Significant"
PM1DEGs$significance[PM1DEGs$p_val_adj < 0.001 & PM1DEGs$avg_log2FC < -0.5] <- "Downregulated"
PM1DEGs$significance[PM1DEGs$p_val_adj < 0.001 & PM1DEGs$avg_log2FC > 0.5] <- "Upregulated"
PM1DEGs <- rownames_to_column(PM1DEGs, var = 'gene')
write.csv(PM1DEGs, file = "C:/Users/nateh/OneDrive/Documents/Mei Lab/RNAseq/scRNAseq/Integrated Datasets/3 + 4 cell type organoids/PM1clustersigDEGs.csv")

FeaturePlot(threeandfour, features = c("TNNI1", "FABP3", "RPL10"), min.cutoff = 'q10', split.by = 'Organoid Type')
FeaturePlot(threeandfour, features = c("ITGA11", "PDGFRB", "FAU"), min.cutoff = 'q10', split.by = 'Organoid Type')


DEgenes$significance <- "Not Significant"
DEgenes$significance[DEgenes$p_val_adj < 0.05 & abs(DEgenes$avg_log2FC) > 1] <- "Significant"


filtered_DEgenes <- subset(DEgenes, avg_log2FC > -2.5 & p_val_adj < 200)
write.csv(filtered_DEgenes, file = "C:/Users/nateh/OneDrive/Documents/Mei Lab/RNAseq/scRNAseq/Integrated Datasets/3 + 4 cell type organoids/top_2000_DEgenes.csv")

# Plot the filtered data
ggplot(filtered_DEgenes, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c("grey", "blue")) +
  theme_minimal() +
  xlab("Log2 Fold Change") +
  ylab("-Log10 Adjusted P-value") +
  ggtitle("Differentially Expressed Genes: 3 vs 4 cell type") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text_repel(data = subset(filtered_DEgenes, p_val_adj < 50 & abs(avg_log2FC) > 1),
                  aes(label = rownames(subset(filtered_DEgenes, p_val_adj < 50 & abs(avg_log2FC) > 1))), 
                  size = 3)

ventricular_marker_response <- FindMarkers(threeandfour, ident.1 = 'Ventricular CMs 1_3 cell type', ident.2 = 'Ventricular CMs 1_4 cell type')

# Create a new column for significance
top_500_DEgenes$significance <- "Not Significant"
top_500_DEgenes$significance[top_500_DEgenes$p_val_adj < 0.05 & abs(top_500_DEgenes$avg_log2FC) > 1.5] <- "Significant"

# Create the volcano plot
ggplot(top_500_DEgenes, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  xlab("Log2 Fold Change") +
  ylab("-Log10 Adjusted P-value") +
  ggtitle("Differentially Expressed Genes: 3 vs 4 cell type") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text_repel(data = subset(top_500_DEgenes, p_val_adj < 0.01 & abs(avg_log2FC) > 1),
                  aes(label = rownames(subset(top_500_DEgenes, p_val_adj < 0.01 & abs(avg_log2FC) > 1))), 
                  size = 3)

write.csv(DEgenes, file = "top_2000_DEgenes.csv")

top_markers <- rownames(head(DEgenes, n = 10))  # Get the top 10 differentially expressed genes
