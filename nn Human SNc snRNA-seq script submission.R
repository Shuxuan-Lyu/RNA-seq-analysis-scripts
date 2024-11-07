library(Seurat)
library(harmony)
library(cowplot)
library(COSG)
library(dplyr)
library(ggplot2)
library(CellChat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(GOSemSim)
library(enrichplot)
library(WGCNA)
library(hdWGCNA)
library(patchwork)


rm(list = ls())
options(stringsAsFactors = F)


# loading nn SNc dataset
nn_scRNA_seq.data <- Read10X(data.dir = "matrix/")
nn_all <- CreateSeuratObject(counts = nn_scRNA_seq.data, project = "nn_scRNA_seq", min.cells = 3)

length(levels(nn_all@meta.data$orig.ident))
write.table(levels(nn_all@meta.data$orig.ident), file = "nn_all_orig.ident.txt", col.names = TRUE, sep = "\t", quote = FALSE)
nn_SNc_meta.data <- read.csv(file = "nn SNc meta.data.csv", header = TRUE)

nn_SNc_batch_ID <- nn_SNc_meta.data$batch
nn_SNc <- subset(nn_all, idents = nn_SNc_batch_ID)

nn_SNc@meta.data$batch_ID <- nn_SNc@active.ident
nn_SNc@meta.data$Sample_ID <- nn_SNc_meta.data[match(nn_SNc@meta.data$batch_ID, nn_SNc_meta.data$batch), "Sample_ID"]
length (unique(nn_SNc$Sample_ID))

nn_SNc@meta.data$Sample_ID <- as.factor(nn_SNc@meta.data$Sample_ID)
levels(nn_SNc@meta.data$Sample_ID)

Idents(nn_SNc) <- "Sample_ID"
nn_SNc <- RenameIdents(nn_SNc, `3298` = "HC", `3322` = "HC", `3345` = "HC", `3346` = "HC", `3482` = "HC", `4956` = "HC", 
                           `5610` = "HC", `6173` = "HC")

nn_SNc <- RenameIdents(nn_SNc, `2544` = "LBD", `2561` = "LBD", `2569` = "LBD", `1963` = "PD", `2142` = "PD", `3873` = "PD", 
                           `3887` = "PD", `4560` = "PD", `4568` = "PD", `4775` = "PD")

levels(nn_SNc@active.ident)
nn_SNc$status <- nn_SNc@active.ident
nn_SNc
## An object of class Seurat 
## 38594 features across 387558 samples within 1 assay 
## Active assay: RNA (38594 features, 0 variable features)

nn_SNc <- PercentageFeatureSet(nn_SNc, pattern = "^MT-", col.name = "percent.mt")
VlnPlot(nn_SNc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(nn_SNc, features = "percent.mt", ncol = 1, group.by = "status", y.max = 15, pt.size = 0)

Idents(nn_SNc) <- "Sample_ID"


# extracting hDANs and ODCs from different SNc samples
# healthy control
## 3322
SNc_3322 <- subset(nn_SNc, idents = "3322")
VlnPlot(SNc_3322, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
SNc_3322 <- NormalizeData(SNc_3322, normalization.method = "LogNormalize", scale.factor = 10000)
SNc_3322 <- FindVariableFeatures(SNc_3322, selection.method = "vst", nfeatures = 3000)
SNc_3322
## An object of class Seurat 
## 38594 features across 49759 samples within 1 assay 
## Active assay: RNA (38594 features, 3000 variable features)

all.genes <- rownames(SNc_3322)
SNc_3322 <- ScaleData(SNc_3322, features = all.genes, verbose = FALSE)
SNc_3322 <- RunPCA(SNc_3322, npcs = 30, verbose = FALSE)

table(SNc_3322$batch_ID)
VizDimLoadings(SNc_3322, dims = 1:2, reduction = "pca")
DimPlot(object = SNc_3322, reduction = "pca", pt.size = .1, group.by = "batch_ID")
VlnPlot(object = SNc_3322, features = "PC_1", group.by = "batch_ID", pt.size = .1)

plot1 <- FeatureScatter(SNc_3322, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(SNc_3322, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

## Run Harmony
Idents(SNc_3322) <- "batch_ID"
SNc_3322$batch_ID <- SNc_3322@active.ident
table(SNc_3322$batch_ID)

SNc_3322 <- RunHarmony(SNc_3322, group.by.vars = "batch_ID", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(SNc_3322, 'harmony')
harmony_embeddings[1:5, 1:5]

DimPlot(object = SNc_3322, reduction = "harmony", pt.size = .1, group.by = "batch_ID")
VlnPlot(object = SNc_3322, features = "harmony_1", group.by = "batch_ID", pt.size = .1)

SNc_3322
## An object of class Seurat 
## 38594 features across 49759 samples within 1 assay 
## Active assay: RNA (38594 features, 3000 variable features)
## 2 dimensional reductions calculated: pca, harmony

## U-MAP and Clustering with harmony
SNc_3322 <- FindNeighbors(SNc_3322, reduction = "harmony", dims = 1:30)
SNc_3322 <- FindClusters(SNc_3322, resolution = 0.1)
SNc_3322 <- RunUMAP(SNc_3322, reduction = "harmony", dims = 1:30)

plot1 <- DimPlot(SNc_3322, reduction = "umap", label = TRUE)
plot2 <- DimPlot(SNc_3322, reduction = "umap", group.by = "batch_ID")
plot_grid(plot1, plot2, ncol = 1)
DimPlot(SNc_3322, reduction = "umap", split.by = "batch_ID", label = TRUE, ncol = 4)

## Assigning cell type identity to clusters
FeaturePlot(SNc_3322, features = c("NR4A2", "TH", "SLC6A3", "SLC18A2"), label = TRUE)
VlnPlot(SNc_3322, features = c("NR4A2", "TH", "SLC6A3", "SLC18A2"), ncol = 1)
FeaturePlot(SNc_3322, features = c("OPALIN", "MOG", "PLP1", "MBP", "CLDN11"), label = TRUE)
FeaturePlot(SNc_3322, features = c("VCAN", "BCAN", "PDGFRA", "OLIG2"), label = TRUE)
FeaturePlot(SNc_3322, features = c("ENO2", "RBFOX3", "SNAP25", "SYT1"), label = TRUE)
VlnPlot(SNc_3322, features = c("ENO2", "RBFOX3", "SNAP25", "SYT1"), ncol = 1)
FeaturePlot(SNc_3322, features = c("AQP4", "GFAP", "SLC1A3", "SLC1A2", "SLC4A4", "NTSR2"), label = TRUE)
FeaturePlot(SNc_3322, features = c("CX3CR1", "P2RY12", "CSF1R"), label = TRUE)
FeaturePlot(SNc_3322, features = c("FLT1", "CLDN5", "PTPRB"), label = TRUE)
FeaturePlot(SNc_3322, features = c("DCN", "COL1A1", "KCNJ8", "ABCC9"), label = TRUE)

SNc_3322 <- RenameIdents(SNc_3322, `0` = "Oligodendrocytes", `1` = "Non-DA neurons", `2` = "Non-DA neurons", `3` = "Microglia", `4` = "Non-DA neurons",
                         `5` = "Endothelial cells", `6` = "Astrocytes", `7` = "Non-DA neurons", `8` = "Non-DA neurons", `9` = "Non-DA neurons",
                         `10` = "Fibroblast", `11` = "Non-DA neurons", `12` = "Non-DA neurons", `13` = "DA neurons", `14` = "Non-DA neurons", `15` = "OPCs")

SNc_3322$Major_celltype <- SNc_3322@active.ident
table(SNc_3322$Major_celltype)

hDANs_3322 <- subset(SNc_3322, idents = "DA neurons")
hDANs_3322.data <- as.data.frame(hDANs_3322@assays$RNA@counts)

hODCs_3322 <- subset(SNc_3322, idents = "Oligodendrocytes")
hODCs_3322.data <- as.data.frame(hODCs_3322@assays$RNA@counts)

## 3298
SNc_3298 <- subset(nn_SNc, idents = "3298")
SNc_3298 <- NormalizeData(SNc_3298, normalization.method = "LogNormalize", scale.factor = 10000)
SNc_3298 <- FindVariableFeatures(SNc_3298, selection.method = "vst", nfeatures = 3000)
SNc_3298
## An object of class Seurat 
## 38594 features across 11579 samples within 1 assay 
## Active assay: RNA (38594 features, 3000 variable features)

all.genes <- rownames(SNc_3298)
SNc_3298 <- ScaleData(SNc_3298, features = all.genes, verbose = FALSE)
SNc_3298 <- RunPCA(SNc_3298, npcs = 30, verbose = FALSE)

Idents(SNc_3298) <- "batch_ID"
SNc_3298$batch_ID <- SNc_3298@active.ident
table(SNc_3298$batch_ID)

SNc_3298 <- RunHarmony(SNc_3298, group.by.vars = "batch_ID", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(SNc_3298, 'harmony')
harmony_embeddings[1:5, 1:5]

SNc_3298
## An object of class Seurat 
## 38594 features across 11579 samples within 1 assay 
## Active assay: RNA (38594 features, 3000 variable features)
## 2 dimensional reductions calculated: pca, harmony

SNc_3298 <- FindNeighbors(SNc_3298, reduction = "harmony", dims = 1:30)
SNc_3298 <- FindClusters(SNc_3298, resolution = 0.1)
SNc_3298 <- RunUMAP(SNc_3298, reduction = "harmony", dims = 1:30)

DimPlot(SNc_3298, reduction = "umap", label = TRUE)
DimPlot(SNc_3298, reduction = "umap", split.by = "batch_ID", label = TRUE, ncol = 4)

FeaturePlot(SNc_3298, features = c("NR4A2", "TH", "SLC6A3", "SLC18A2"), label = TRUE)
VlnPlot(SNc_3298, features = c("NR4A2", "TH", "SLC6A3", "SLC18A2"), ncol = 1)
FeaturePlot(SNc_3298, features = c("OPALIN", "MOG", "PLP1", "MBP", "CLDN11"), label = TRUE)
FeaturePlot(SNc_3298, features = c("VCAN", "BCAN", "PDGFRA", "OLIG2"), label = TRUE)
FeaturePlot(SNc_3298, features = c("ENO2", "RBFOX3", "SNAP25", "SYT1"), label = TRUE)
VlnPlot(SNc_3298, features = c("ENO2", "RBFOX3", "SNAP25", "SYT1"), ncol = 1)
FeaturePlot(SNc_3298, features = c("AQP4", "GFAP", "SLC1A3", "SLC1A2", "SLC4A4", "NTSR2"), label = TRUE)
FeaturePlot(SNc_3298, features = c("CX3CR1", "P2RY12", "CSF1R"), label = TRUE)
FeaturePlot(SNc_3298, features = c("FLT1", "CLDN5", "PTPRB"), label = TRUE)
FeaturePlot(SNc_3298, features = c("DCN", "COL1A1", "KCNJ8", "ABCC9"), label = TRUE)

SNc_3298 <- RenameIdents(SNc_3298, `0` = "Oligodendrocytes", `1` = "Astrocytes", `2` = "OPCs", `3` = "Oligodendrocytes", `4` = "Microglia",
                         `5` = "Endothelial cells", `6` = "Non-DA neurons", `7` = "Fibroblast", `8` = "DA neurons", `9` = "Non-DA neurons")

SNc_3298$Major_celltype <- SNc_3298@active.ident
table(SNc_3298$Major_celltype)

hDANs_3298 <- subset(SNc_3298, idents = "DA neurons")
hDANs_3298.data <- as.data.frame(hDANs_3298@assays$RNA@counts)

hODCs_3298 <- subset(SNc_3298, idents = "Oligodendrocytes")
hODCs_3298.data <- as.data.frame(hODCs_3298@assays$RNA@counts)

## 4956
SNc_4956 <- subset(nn_SNc, idents = "4956")
SNc_4956 <- NormalizeData(SNc_4956, normalization.method = "LogNormalize", scale.factor = 10000)
SNc_4956 <- FindVariableFeatures(SNc_4956, selection.method = "vst", nfeatures = 3000)
SNc_4956
## An object of class Seurat 
## 38594 features across 25277 samples within 1 assay 
## Active assay: RNA (38594 features, 3000 variable features)

all.genes <- rownames(SNc_4956)
SNc_4956 <- ScaleData(SNc_4956, features = all.genes, verbose = FALSE)
SNc_4956 <- RunPCA(SNc_4956, npcs = 30, verbose = FALSE)

Idents(SNc_4956) <- "batch_ID"
SNc_4956$batch_ID <- SNc_4956@active.ident
table(SNc_4956$batch_ID)

SNc_4956 <- RunHarmony(SNc_4956, group.by.vars = "batch_ID", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(SNc_4956, 'harmony')
harmony_embeddings[1:5, 1:5]

SNc_4956
## An object of class Seurat 
## 38594 features across 25277 samples within 1 assay 
## Active assay: RNA (38594 features, 3000 variable features)
## 2 dimensional reductions calculated: pca, harmony

SNc_4956 <- FindNeighbors(SNc_4956, reduction = "harmony", dims = 1:30)
SNc_4956 <- FindClusters(SNc_4956, resolution = 0.2)
SNc_4956 <- RunUMAP(SNc_4956, reduction = "harmony", dims = 1:30)

DimPlot(SNc_4956, reduction = "umap", label = TRUE)
DimPlot(SNc_4956, reduction = "umap", split.by = "batch_ID", label = TRUE, ncol = 4)

FeaturePlot(SNc_4956, features = c("NR4A2", "TH", "SLC6A3", "SLC18A2"), label = TRUE)
VlnPlot(SNc_4956, features = c("NR4A2", "TH", "SLC6A3", "SLC18A2"), ncol = 1)
FeaturePlot(SNc_4956, features = c("OPALIN", "MOG", "PLP1", "MBP", "CLDN11"), label = TRUE)
FeaturePlot(SNc_4956, features = c("VCAN", "BCAN", "PDGFRA", "OLIG2"), label = TRUE)
FeaturePlot(SNc_4956, features = c("ENO2", "RBFOX3", "SNAP25", "SYT1"), label = TRUE)
VlnPlot(SNc_4956, features = c("ENO2", "RBFOX3", "SNAP25", "SYT1"), ncol = 1)
FeaturePlot(SNc_4956, features = c("AQP4", "GFAP", "SLC1A3", "SLC1A2", "SLC4A4", "NTSR2"), label = TRUE)
FeaturePlot(SNc_4956, features = c("CX3CR1", "P2RY12", "CSF1R"), label = TRUE)
FeaturePlot(SNc_4956, features = c("FLT1", "CLDN5", "PTPRB"), label = TRUE)
FeaturePlot(SNc_4956, features = c("DCN", "COL1A1", "KCNJ8", "ABCC9"), label = TRUE)

SNc_4956 <- RenameIdents(SNc_4956, `0` = "DA neurons", `1` = "Oligodendrocytes", `2` = "Astrocytes", `3` = "Microglia", `4` = "Oligodendrocytes",
                         `5` = "DA neurons", `6` = "Oligodendrocytes", `7` = "Non-DA neurons", `8` = "Non-DA neurons", `9` = "Endothelial cells",
                         `10` = "Microglia", `11` = "Non-DA neurons", `12` = "Non-DA neurons", `13` = "OPCs")

SNc_4956$Major_celltype <- SNc_4956@active.ident
table(SNc_4956$Major_celltype)

hDANs_4956 <- subset(SNc_4956, idents = "DA neurons")
hDANs_4956.data <- as.data.frame(hDANs_4956@assays$RNA@counts)

hODCs_4956 <- subset(SNc_4956, idents = "Oligodendrocytes")
hODCs_4956.data <- as.data.frame(hODCs_4956@assays$RNA@counts)

## 3345
SNc_3345 <- subset(nn_SNc, idents = "3345")
SNc_3345 <- NormalizeData(SNc_3345, normalization.method = "LogNormalize", scale.factor = 10000)
SNc_3345 <- FindVariableFeatures(SNc_3345, selection.method = "vst", nfeatures = 3000)
SNc_3345
## An object of class Seurat 
## 38594 features across 29230 samples within 1 assay 
## Active assay: RNA (38594 features, 3000 variable features)

all.genes <- rownames(SNc_3345)
SNc_3345 <- ScaleData(SNc_3345, features = all.genes, verbose = FALSE)
SNc_3345 <- RunPCA(SNc_3345, npcs = 30, verbose = FALSE)

Idents(SNc_3345) <- "batch_ID"
SNc_3345$batch_ID <- SNc_3345@active.ident
table(SNc_3345$batch_ID)

SNc_3345 <- RunHarmony(SNc_3345, group.by.vars = "batch_ID", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(SNc_3345, 'harmony')
harmony_embeddings[1:5, 1:5]

SNc_3345
## An object of class Seurat 
## 38594 features across 29230 samples within 1 assay 
## Active assay: RNA (38594 features, 3000 variable features)
## 2 dimensional reductions calculated: pca, harmony

SNc_3345 <- FindNeighbors(SNc_3345, reduction = "harmony", dims = 1:30)
SNc_3345 <- FindClusters(SNc_3345, resolution = 0.2)
SNc_3345 <- RunUMAP(SNc_3345, reduction = "harmony", dims = 1:30)

DimPlot(SNc_3345, reduction = "umap", label = TRUE)
DimPlot(SNc_3345, reduction = "umap", split.by = "batch_ID", label = TRUE, ncol = 4)

FeaturePlot(SNc_3345, features = c("NR4A2", "TH", "SLC6A3", "SLC18A2"), label = TRUE)
VlnPlot(SNc_3345, features = c("NR4A2", "TH", "SLC6A3", "SLC18A2"), ncol = 1)
FeaturePlot(SNc_3345, features = c("OPALIN", "MOG", "PLP1", "MBP", "CLDN11"), label = TRUE)
FeaturePlot(SNc_3345, features = c("VCAN", "BCAN", "PDGFRA", "OLIG2"), label = TRUE)
FeaturePlot(SNc_3345, features = c("ENO2", "RBFOX3", "SNAP25", "SYT1"), label = TRUE)
VlnPlot(SNc_3345, features = c("ENO2", "RBFOX3", "SNAP25", "SYT1"), ncol = 1)
FeaturePlot(SNc_3345, features = c("AQP4", "GFAP", "SLC1A3", "SLC1A2", "SLC4A4", "NTSR2"), label = TRUE)
FeaturePlot(SNc_3345, features = c("CX3CR1", "P2RY12", "CSF1R"), label = TRUE)
FeaturePlot(SNc_3345, features = c("FLT1", "CLDN5", "PTPRB"), label = TRUE)
FeaturePlot(SNc_3345, features = c("DCN", "COL1A1", "KCNJ8", "ABCC9"), label = TRUE)

SNc_3345 <- RenameIdents(SNc_3345, `0` = "Oligodendrocytes", `1` = "Oligodendrocytes", `2` = "Microglia", `3` = "Astrocytes", `4` = "DA neurons",
                         `5` = "Non-DA neurons", `6` = "OPCs", `7` = "Endothelial cells", `8` = "Non-DA neurons", `9` = "Non-DA neurons",
                         `10` = "Microglia", `11` = "Fibroblast", `12` = "Microglia")

SNc_3345$Major_celltype <- SNc_3345@active.ident
table(SNc_3345$Major_celltype)

hDANs_3345 <- subset(SNc_3345, idents = "DA neurons")
hDANs_3345.data <- as.data.frame(hDANs_3345@assays$RNA@counts)

hODCs_3345 <- subset(SNc_3345, idents = "Oligodendrocytes")
hODCs_3345.data <- as.data.frame(hODCs_3345@assays$RNA@counts)                      

## 3346
SNc_3346 <- subset(nn_SNc, idents = "3346")
SNc_3346 <- NormalizeData(SNc_3346, normalization.method = "LogNormalize", scale.factor = 10000)
SNc_3346 <- FindVariableFeatures(SNc_3346, selection.method = "vst", nfeatures = 3000)
SNc_3346
## An object of class Seurat 
## 38594 features across 23747 samples within 1 assay 
## Active assay: RNA (38594 features, 3000 variable features)

all.genes <- rownames(SNc_3346)
SNc_3346 <- ScaleData(SNc_3346, features = all.genes, verbose = FALSE)
SNc_3346 <- RunPCA(SNc_3346, npcs = 30, verbose = FALSE)

Idents(SNc_3346) <- "batch_ID"
SNc_3346$batch_ID <- SNc_3346@active.ident
table(SNc_3346$batch_ID)

SNc_3346 <- RunHarmony(SNc_3346, group.by.vars = "batch_ID", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(SNc_3346, 'harmony')
harmony_embeddings[1:5, 1:5]

SNc_3346
## An object of class Seurat 
## 38594 features across 23747 samples within 1 assay 
## Active assay: RNA (38594 features, 3000 variable features)
## 2 dimensional reductions calculated: pca, harmony

SNc_3346 <- FindNeighbors(SNc_3346, reduction = "harmony", dims = 1:30)
SNc_3346 <- FindClusters(SNc_3346, resolution = 0.1)
SNc_3346 <- RunUMAP(SNc_3346, reduction = "harmony", dims = 1:30)

DimPlot(SNc_3346, reduction = "umap", label = TRUE)
DimPlot(SNc_3346, reduction = "umap", split.by = "batch_ID", label = TRUE, ncol = 4)

FeaturePlot(SNc_3346, features = c("NR4A2", "TH", "SLC6A3", "SLC18A2"), label = TRUE)
VlnPlot(SNc_3346, features = c("NR4A2", "TH", "SLC6A3", "SLC18A2"), ncol = 1)
FeaturePlot(SNc_3346, features = c("OPALIN", "MOG", "PLP1", "MBP", "CLDN11"), label = TRUE)
FeaturePlot(SNc_3346, features = c("VCAN", "BCAN", "PDGFRA", "OLIG2"), label = TRUE)
FeaturePlot(SNc_3346, features = c("ENO2", "RBFOX3", "SNAP25", "SYT1"), label = TRUE)
VlnPlot(SNc_3346, features = c("ENO2", "RBFOX3", "SNAP25", "SYT1"), ncol = 1)
FeaturePlot(SNc_3346, features = c("AQP4", "GFAP", "SLC1A3", "SLC1A2", "SLC4A4", "NTSR2"), label = TRUE)
FeaturePlot(SNc_3346, features = c("CX3CR1", "P2RY12", "CSF1R"), label = TRUE)
FeaturePlot(SNc_3346, features = c("FLT1", "CLDN5", "PTPRB"), label = TRUE)
FeaturePlot(SNc_3346, features = c("DCN", "COL1A1", "KCNJ8", "ABCC9"), label = TRUE)

SNc_3346 <- RenameIdents(SNc_3346, `0` = "Oligodendrocytes", `1` = "Non-DA neurons", `2` = "Non-DA neurons", `3` = "Microglia", `4` = "DA neurons",
                         `5` = "Endothelial cells", `6` = "Astrocytes", `7` = "Non-DA neurons", `8` = "OPCs", `9` = "Non-DA neurons",
                         `10` = "Non-DA neurons", `11` = "Non-DA neurons", `12` = "Non-DA neurons", `13` = "Non-DA neurons")

SNc_3346$Major_celltype <- SNc_3346@active.ident
table(SNc_3346$Major_celltype)

hDANs_3346 <- subset(SNc_3346, idents = "DA neurons")
hDANs_3346.data <- as.data.frame(hDANs_3346@assays$RNA@counts)

hODCs_3346 <- subset(SNc_3346, idents = "Oligodendrocytes")
hODCs_3346.data <- as.data.frame(hODCs_3346@assays$RNA@counts) 

## 3482
SNc_3482 <- subset(nn_SNc, idents = "3482")
SNc_3482 <- NormalizeData(SNc_3482, normalization.method = "LogNormalize", scale.factor = 10000)
SNc_3482 <- FindVariableFeatures(SNc_3482, selection.method = "vst", nfeatures = 3000)
SNc_3482
## An object of class Seurat 
## 38594 features across 24856 samples within 1 assay 
## Active assay: RNA (38594 features, 3000 variable features)

all.genes <- rownames(SNc_3482)
SNc_3482 <- ScaleData(SNc_3482, features = all.genes, verbose = FALSE)
SNc_3482 <- RunPCA(SNc_3482, npcs = 30, verbose = FALSE)

Idents(SNc_3482) <- "batch_ID"
SNc_3482$batch_ID <- SNc_3482@active.ident
table(SNc_3482$batch_ID)

SNc_3482 <- RunHarmony(SNc_3482, group.by.vars = "batch_ID", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(SNc_3482, 'harmony')
harmony_embeddings[1:5, 1:5]

SNc_3482
## An object of class Seurat 
## 38594 features across 24856 samples within 1 assay 
## Active assay: RNA (38594 features, 3000 variable features)
## 2 dimensional reductions calculated: pca, harmony

SNc_3482 <- FindNeighbors(SNc_3482, reduction = "harmony", dims = 1:30)
SNc_3482 <- FindClusters(SNc_3482, resolution = 0.1)
SNc_3482 <- RunUMAP(SNc_3482, reduction = "harmony", dims = 1:30)

DimPlot(SNc_3482, reduction = "umap", label = TRUE)
DimPlot(SNc_3482, reduction = "umap", split.by = "batch_ID", label = TRUE, ncol = 4)

FeaturePlot(SNc_3482, features = c("NR4A2", "TH", "SLC6A3", "SLC18A2"), label = TRUE)
VlnPlot(SNc_3482, features = c("NR4A2", "TH", "SLC6A3", "SLC18A2"), ncol = 1)
FeaturePlot(SNc_3482, features = c("OPALIN", "MOG", "PLP1", "MBP", "CLDN11"), label = TRUE)
FeaturePlot(SNc_3482, features = c("VCAN", "BCAN", "PDGFRA", "OLIG2"), label = TRUE)
FeaturePlot(SNc_3482, features = c("ENO2", "RBFOX3", "SNAP25", "SYT1"), label = TRUE)
VlnPlot(SNc_3482, features = c("ENO2", "RBFOX3", "SNAP25", "SYT1"), ncol = 1)
FeaturePlot(SNc_3482, features = c("AQP4", "GFAP", "SLC1A3", "SLC1A2", "SLC4A4", "NTSR2"), label = TRUE)
FeaturePlot(SNc_3482, features = c("CX3CR1", "P2RY12", "CSF1R"), label = TRUE)
FeaturePlot(SNc_3482, features = c("FLT1", "CLDN5", "PTPRB"), label = TRUE)
FeaturePlot(SNc_3482, features = c("DCN", "COL1A1", "KCNJ8", "ABCC9"), label = TRUE)

SNc_3482 <- RenameIdents(SNc_3482, `0` = "Oligodendrocytes", `1` = "Astrocytes", `2` = "DA neurons", `3` = "Non-DA neurons", `4` = "Microglia",
                         `5` = "Endothelial cells", `6` = "OPCs", `7` = "Non-DA neurons", `8` = "Non-DA neurons", `9` = "Non-DA neurons",
                         `10` = "Non-DA neurons", `11` = "Fibroblast", `12` = "Non-DA neurons", `13` = "Non-DA neurons", `14` = "Non-DA neurons",
                         `15` = "Non-DA neurons")

SNc_3482$Major_celltype <- SNc_3482@active.ident
table(SNc_3482$Major_celltype)

hDANs_3482 <- subset(SNc_3482, idents = "DA neurons")
hDANs_3482.data <- as.data.frame(hDANs_3482@assays$RNA@counts)

hODCs_3482 <- subset(SNc_3482, idents = "Oligodendrocytes")
hODCs_3482.data <- as.data.frame(hODCs_3482@assays$RNA@counts) 

## 5610
SNc_5610 <- subset(nn_SNc, idents = "5610")
SNc_5610 <- NormalizeData(SNc_5610, normalization.method = "LogNormalize", scale.factor = 10000)
SNc_5610 <- FindVariableFeatures(SNc_5610, selection.method = "vst", nfeatures = 3000)
SNc_5610
## An object of class Seurat 
## 38594 features across 6970 samples within 1 assay 
## Active assay: RNA (38594 features, 3000 variable features)

all.genes <- rownames(SNc_5610)
SNc_5610 <- ScaleData(SNc_5610, features = all.genes, verbose = FALSE)
SNc_5610 <- RunPCA(SNc_5610, npcs = 30, verbose = FALSE)

Idents(SNc_5610) <- "batch_ID"
SNc_5610$batch_ID <- SNc_5610@active.ident
table(SNc_5610$batch_ID)

SNc_5610 <- RunHarmony(SNc_5610, group.by.vars = "batch_ID", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(SNc_5610, 'harmony')
harmony_embeddings[1:5, 1:5]

SNc_5610
## An object of class Seurat 
## 38594 features across 6970 samples within 1 assay 
## Active assay: RNA (38594 features, 3000 variable features)
## 2 dimensional reductions calculated: pca, harmony

SNc_5610 <- FindNeighbors(SNc_5610, reduction = "harmony", dims = 1:30)
SNc_5610 <- FindClusters(SNc_5610, resolution = 0.4)
SNc_5610 <- RunUMAP(SNc_5610, reduction = "harmony", dims = 1:30)

DimPlot(SNc_5610, reduction = "umap", label = TRUE)
DimPlot(SNc_5610, reduction = "umap", split.by = "batch_ID", label = TRUE, ncol = 4)

FeaturePlot(SNc_5610, features = c("NR4A2", "TH", "SLC6A3", "SLC18A2"), label = TRUE)
VlnPlot(SNc_5610, features = c("NR4A2", "TH", "SLC6A3", "SLC18A2"), ncol = 1)
FeaturePlot(SNc_5610, features = c("OPALIN", "MOG", "PLP1", "MBP", "CLDN11"), label = TRUE)
FeaturePlot(SNc_5610, features = c("VCAN", "BCAN", "PDGFRA", "OLIG2"), label = TRUE)
FeaturePlot(SNc_5610, features = c("ENO2", "RBFOX3", "SNAP25", "SYT1"), label = TRUE)
VlnPlot(SNc_5610, features = c("ENO2", "RBFOX3", "SNAP25", "SYT1"), ncol = 1)
FeaturePlot(SNc_5610, features = c("AQP4", "GFAP", "SLC1A3", "SLC1A2", "SLC4A4", "NTSR2"), label = TRUE)
FeaturePlot(SNc_5610, features = c("CX3CR1", "P2RY12", "CSF1R"), label = TRUE)
FeaturePlot(SNc_5610, features = c("FLT1", "CLDN5", "PTPRB"), label = TRUE)
FeaturePlot(SNc_5610, features = c("DCN", "COL1A1", "KCNJ8", "ABCC9"), label = TRUE)

SNc_5610 <- RenameIdents(SNc_5610, `0` = "Oligodendrocytes", `1` = "Non-DA neurons", `2` = "Astrocytes", `3` = "Oligodendrocytes", `4` = "Oligodendrocytes",
                         `5` = "DA neurons", `6` = "Microglia", `7` = "Non-DA neurons", `8` = "Non-DA neurons", `9` = "OPCs", `10` = "DA neurons", 
                         `11` = "Astrocytes", `12` = "Fibroblast", `13` = "DA neurons", `14` = "Endothelial cells", `15` = "Non-DA neurons",
                         `16` = "Non-DA neurons")

SNc_5610$Major_celltype <- SNc_5610@active.ident
table(SNc_5610$Major_celltype)

hDANs_5610 <- subset(SNc_5610, idents = "DA neurons")
hDANs_5610.data <- as.data.frame(hDANs_5610@assays$RNA@counts)

hODCs_5610 <- subset(SNc_5610, idents = "Oligodendrocytes")
hODCs_5610.data <- as.data.frame(hODCs_5610@assays$RNA@counts) 

## 6173
SNc_6173 <- subset(nn_SNc, idents = "6173")
SNc_6173 <- NormalizeData(SNc_6173, normalization.method = "LogNormalize", scale.factor = 10000)
SNc_6173 <- FindVariableFeatures(SNc_6173, selection.method = "vst", nfeatures = 3000)
SNc_6173
## An object of class Seurat 
## 38594 features across 13330 samples within 1 assay 
## Active assay: RNA (38594 features, 3000 variable features)

all.genes <- rownames(SNc_6173)
SNc_6173 <- ScaleData(SNc_6173, features = all.genes, verbose = FALSE)
SNc_6173 <- RunPCA(SNc_6173, npcs = 30, verbose = FALSE)

Idents(SNc_6173) <- "batch_ID"
SNc_6173$batch_ID <- SNc_6173@active.ident
table(SNc_6173$batch_ID)

SNc_6173 <- RunHarmony(SNc_6173, group.by.vars = "batch_ID", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(SNc_6173, 'harmony')
harmony_embeddings[1:5, 1:5]

SNc_6173
## An object of class Seurat 
## 38594 features across 13330 samples within 1 assay 
## Active assay: RNA (38594 features, 3000 variable features)
## 2 dimensional reductions calculated: pca, harmony

SNc_6173 <- FindNeighbors(SNc_6173, reduction = "harmony", dims = 1:30)
SNc_6173 <- FindClusters(SNc_6173, resolution = 0.1)
SNc_6173 <- RunUMAP(SNc_6173, reduction = "harmony", dims = 1:30)

DimPlot(SNc_6173, reduction = "umap", label = TRUE)
DimPlot(SNc_6173, reduction = "umap", split.by = "batch_ID", label = TRUE, ncol = 4)

FeaturePlot(SNc_6173, features = c("NR4A2", "TH", "SLC6A3", "SLC18A2"), label = TRUE)
VlnPlot(SNc_6173, features = c("NR4A2", "TH", "SLC6A3", "SLC18A2"), ncol = 1)
FeaturePlot(SNc_6173, features = c("OPALIN", "MOG", "PLP1", "MBP", "CLDN11"), label = TRUE)
FeaturePlot(SNc_6173, features = c("VCAN", "BCAN", "PDGFRA", "OLIG2"), label = TRUE)
FeaturePlot(SNc_6173, features = c("ENO2", "RBFOX3", "SNAP25", "SYT1"), label = TRUE)
VlnPlot(SNc_6173, features = c("ENO2", "RBFOX3", "SNAP25", "SYT1"), ncol = 1)
FeaturePlot(SNc_6173, features = c("AQP4", "GFAP", "SLC1A3", "SLC1A2", "SLC4A4", "NTSR2"), label = TRUE)
FeaturePlot(SNc_6173, features = c("CX3CR1", "P2RY12", "CSF1R"), label = TRUE)
FeaturePlot(SNc_6173, features = c("FLT1", "CLDN5", "PTPRB"), label = TRUE)
FeaturePlot(SNc_6173, features = c("DCN", "COL1A1", "KCNJ8", "ABCC9"), label = TRUE)

SNc_6173 <- RenameIdents(SNc_6173, `0` = "Oligodendrocytes", `1` = "DA neurons", `2` = "Non-DA neurons", `3` = "Microglia", `4` = "Non-DA neurons",
                         `5` = "Astrocytes", `6` = "OPCs", `7` = "Non-DA neurons", `8` = "Endothelial cells", `9` = "Fibroblast",
                         `10` = "Non-DA neurons", `11` = "Non-DA neurons")

SNc_6173$Major_celltype <- SNc_6173@active.ident
table(SNc_6173$Major_celltype)

hDANs_6173 <- subset(SNc_6173, idents = "DA neurons")
hDANs_6173.data <- as.data.frame(hDANs_6173@assays$RNA@counts)

hODCs_6173 <- subset(SNc_6173, idents = "Oligodendrocytes")
hODCs_6173.data <- as.data.frame(hODCs_6173@assays$RNA@counts) 




# PD
## 4560
SNc_4560 <- subset(nn_SNc, idents = "4560")
SNc_4560 <- NormalizeData(SNc_4560, normalization.method = "LogNormalize", scale.factor = 10000)
SNc_4560 <- FindVariableFeatures(SNc_4560, selection.method = "vst", nfeatures = 3000)
SNc_4560
## An object of class Seurat 
## 38594 features across 38335 samples within 1 assay 
## Active assay: RNA (38594 features, 3000 variable features)

all.genes <- rownames(SNc_4560)
SNc_4560 <- ScaleData(SNc_4560, features = all.genes, verbose = FALSE)
SNc_4560 <- RunPCA(SNc_4560, npcs = 30, verbose = FALSE)

Idents(SNc_4560) <- "batch_ID"
SNc_4560$batch_ID <- SNc_4560@active.ident
table(SNc_4560$batch_ID)

SNc_4560 <- RunHarmony(SNc_4560, group.by.vars = "batch_ID", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(SNc_4560, 'harmony')
harmony_embeddings[1:5, 1:5]

SNc_4560
## An object of class Seurat 
## 38594 features across 38335 samples within 1 assay 
## Active assay: RNA (38594 features, 3000 variable features)
## 2 dimensional reductions calculated: pca, harmony

SNc_4560 <- FindNeighbors(SNc_4560, reduction = "harmony", dims = 1:30)
SNc_4560 <- FindClusters(SNc_4560, resolution = 0.1)
SNc_4560 <- RunUMAP(SNc_4560, reduction = "harmony", dims = 1:30)

DimPlot(SNc_4560, reduction = "umap", label = TRUE)
DimPlot(SNc_4560, reduction = "umap", split.by = "batch_ID", label = TRUE, ncol = 4)

FeaturePlot(SNc_4560, features = c("NR4A2", "TH", "SLC6A3", "SLC18A2"), label = TRUE)
VlnPlot(SNc_4560, features = c("NR4A2", "TH", "SLC6A3", "SLC18A2"), ncol = 1)
FeaturePlot(SNc_4560, features = c("OPALIN", "MOG", "PLP1", "MBP", "CLDN11"), label = TRUE)
FeaturePlot(SNc_4560, features = c("VCAN", "BCAN", "PDGFRA", "OLIG2"), label = TRUE)
FeaturePlot(SNc_4560, features = c("ENO2", "RBFOX3", "SNAP25", "SYT1"), label = TRUE)
VlnPlot(SNc_4560, features = c("ENO2", "RBFOX3", "SNAP25", "SYT1"), ncol = 1)
FeaturePlot(SNc_4560, features = c("AQP4", "GFAP", "SLC1A3", "SLC1A2", "SLC4A4", "NTSR2"), label = TRUE)
FeaturePlot(SNc_4560, features = c("CX3CR1", "P2RY12", "CSF1R"), label = TRUE)
FeaturePlot(SNc_4560, features = c("FLT1", "CLDN5", "PTPRB"), label = TRUE)
FeaturePlot(SNc_4560, features = c("DCN", "COL1A1", "KCNJ8", "ABCC9"), label = TRUE)

SNc_4560 <- RenameIdents(SNc_4560, `0` = "Oligodendrocytes", `1` = "Non-DA neurons", `2` = "Oligodendrocytes", `3` = "Astrocytes", `4` = "Microglia",
                         `5` = "DA neurons", `6` = "Non-DA neurons", `7` = "Endothelial cells", `8` = "Non-DA neurons", `9` = "Non-DA neurons",
                         `10` = "OPCs", `11` = "Non-DA neurons", `12` = "Non-DA neurons", `13` = "Non-DA neurons", `14` = "Non-DA neurons", `15` = "Fibroblast")

SNc_4560$Major_celltype <- SNc_4560@active.ident
table(SNc_4560$Major_celltype)

hDANs_4560 <- subset(SNc_4560, idents = "DA neurons")
hDANs_4560.data <- as.data.frame(hDANs_4560@assays$RNA@counts)

hODCs_4560 <- subset(SNc_4560, idents = "Oligodendrocytes")
hODCs_4560.data <- as.data.frame(hODCs_4560@assays$RNA@counts) 

## 4568
SNc_4568 <- subset(nn_SNc, idents = "4568")
SNc_4568 <- NormalizeData(SNc_4568, normalization.method = "LogNormalize", scale.factor = 10000)
SNc_4568 <- FindVariableFeatures(SNc_4568, selection.method = "vst", nfeatures = 3000)
SNc_4568
## An object of class Seurat 
## 38594 features across 21625 samples within 1 assay 
## Active assay: RNA (38594 features, 3000 variable features)

all.genes <- rownames(SNc_4568)
SNc_4568 <- ScaleData(SNc_4568, features = all.genes, verbose = FALSE)
SNc_4568 <- RunPCA(SNc_4568, npcs = 30, verbose = FALSE)

Idents(SNc_4568) <- "batch_ID"
SNc_4568$batch_ID <- SNc_4568@active.ident
table(SNc_4568$batch_ID)

SNc_4568 <- RunHarmony(SNc_4568, group.by.vars = "batch_ID", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(SNc_4568, 'harmony')
harmony_embeddings[1:5, 1:5]

SNc_4568
## An object of class Seurat 
## 38594 features across 21625 samples within 1 assay 
## Active assay: RNA (38594 features, 3000 variable features)
## 2 dimensional reductions calculated: pca, harmony

SNc_4568 <- FindNeighbors(SNc_4568, reduction = "harmony", dims = 1:30)
SNc_4568 <- FindClusters(SNc_4568, resolution = 0.2)
SNc_4568 <- RunUMAP(SNc_4568, reduction = "harmony", dims = 1:30)

DimPlot(SNc_4568, reduction = "umap", label = TRUE)
DimPlot(SNc_4568, reduction = "umap", split.by = "batch_ID", label = TRUE, ncol = 4)

FeaturePlot(SNc_4568, features = c("NR4A2", "TH", "SLC6A3", "SLC18A2"), label = TRUE)
VlnPlot(SNc_4568, features = c("NR4A2", "TH", "SLC6A3", "SLC18A2"), ncol = 1)
FeaturePlot(SNc_4568, features = c("OPALIN", "MOG", "PLP1", "MBP", "CLDN11"), label = TRUE)
FeaturePlot(SNc_4568, features = c("VCAN", "BCAN", "PDGFRA", "OLIG2"), label = TRUE)
FeaturePlot(SNc_4568, features = c("ENO2", "RBFOX3", "SNAP25", "SYT1"), label = TRUE)
VlnPlot(SNc_4568, features = c("ENO2", "RBFOX3", "SNAP25", "SYT1"), ncol = 1)
FeaturePlot(SNc_4568, features = c("AQP4", "GFAP", "SLC1A3", "SLC1A2", "SLC4A4", "NTSR2"), label = TRUE)
FeaturePlot(SNc_4568, features = c("CX3CR1", "P2RY12", "CSF1R"), label = TRUE)
FeaturePlot(SNc_4568, features = c("FLT1", "CLDN5", "PTPRB"), label = TRUE)
FeaturePlot(SNc_4568, features = c("DCN", "COL1A1", "KCNJ8", "ABCC9"), label = TRUE)

SNc_4568 <- RenameIdents(SNc_4568, `0` = "Oligodendrocytes", `1` = "Oligodendrocytes", `2` = "Non-DA neurons", `3` = "Astrocytes", `4` = "Non-DA neurons",
                         `5` = "Endothelial cells", `6` = "Oligodendrocytes", `7` = "Microglia", `8` = "OPCs", `9` = "DA neurons", `10` = "Non-DA neurons", 
                         `11` = "Non-DA neurons", `12` = "Fibroblast", `13` = "Non-DA neurons", `14` = "Non-DA neurons", `15` = "Non-DA neurons")

SNc_4568$Major_celltype <- SNc_4568@active.ident
table(SNc_4568$Major_celltype)

hDANs_4568 <- subset(SNc_4568, idents = "DA neurons")
hDANs_4568.data <- as.data.frame(hDANs_4568@assays$RNA@counts)

hODCs_4568 <- subset(SNc_4568, idents = "Oligodendrocytes")
hODCs_4568.data <- as.data.frame(hODCs_4568@assays$RNA@counts) 

## 3873
SNc_3873 <- subset(nn_SNc, idents = "3873")
SNc_3873 <- NormalizeData(SNc_3873, normalization.method = "LogNormalize", scale.factor = 10000)
SNc_3873 <- FindVariableFeatures(SNc_3873, selection.method = "vst", nfeatures = 3000)
SNc_3873
## An object of class Seurat 
## 38594 features across 28521 samples within 1 assay 
## Active assay: RNA (38594 features, 3000 variable features)

all.genes <- rownames(SNc_3873)
SNc_3873 <- ScaleData(SNc_3873, features = all.genes, verbose = FALSE)
SNc_3873 <- RunPCA(SNc_3873, npcs = 30, verbose = FALSE)

Idents(SNc_3873) <- "batch_ID"
SNc_3873$batch_ID <- SNc_3873@active.ident
table(SNc_3873$batch_ID)

SNc_3873 <- RunHarmony(SNc_3873, group.by.vars = "batch_ID", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(SNc_3873, 'harmony')
harmony_embeddings[1:5, 1:5]

SNc_3873
## An object of class Seurat 
## 38594 features across 28521 samples within 1 assay 
## Active assay: RNA (38594 features, 3000 variable features)
## 2 dimensional reductions calculated: pca, harmony

SNc_3873 <- FindNeighbors(SNc_3873, reduction = "harmony", dims = 1:30)
SNc_3873 <- FindClusters(SNc_3873, resolution = 0.1)
SNc_3873 <- RunUMAP(SNc_3873, reduction = "harmony", dims = 1:30)

DimPlot(SNc_3873, reduction = "umap", label = TRUE)
DimPlot(SNc_3873, reduction = "umap", split.by = "batch_ID", label = TRUE, ncol = 4)

FeaturePlot(SNc_3873, features = c("NR4A2", "TH", "SLC6A3", "SLC18A2"), label = TRUE)
VlnPlot(SNc_3873, features = c("NR4A2", "TH", "SLC6A3", "SLC18A2"), ncol = 1)
FeaturePlot(SNc_3873, features = c("OPALIN", "MOG", "PLP1", "MBP", "CLDN11"), label = TRUE)
FeaturePlot(SNc_3873, features = c("VCAN", "BCAN", "PDGFRA", "OLIG2"), label = TRUE)
FeaturePlot(SNc_3873, features = c("ENO2", "RBFOX3", "SNAP25", "SYT1"), label = TRUE)
VlnPlot(SNc_3873, features = c("ENO2", "RBFOX3", "SNAP25", "SYT1"), ncol = 1)
FeaturePlot(SNc_3873, features = c("AQP4", "GFAP", "SLC1A3", "SLC1A2", "SLC4A4", "NTSR2"), label = TRUE)
FeaturePlot(SNc_3873, features = c("CX3CR1", "P2RY12", "CSF1R"), label = TRUE)
FeaturePlot(SNc_3873, features = c("FLT1", "CLDN5", "PTPRB"), label = TRUE)
FeaturePlot(SNc_3873, features = c("DCN", "COL1A1", "KCNJ8", "ABCC9"), label = TRUE)

SNc_3873 <- RenameIdents(SNc_3873, `0` = "Oligodendrocytes", `1` = "Non-DA neurons", `2` = "Astrocytes", `3` = "Non-DA neurons", `4` = "Non-DA neurons",
                         `5` = "Microglia", `6` = "OPCs", `7` = "Non-DA neurons", `8` = "Endothelial cells", `9` = "Non-DA neurons", `10` = "DA neurons", 
                         `11` = "Fibroblast")

SNc_3873$Major_celltype <- SNc_3873@active.ident
table(SNc_3873$Major_celltype)

hDANs_3873 <- subset(SNc_3873, idents = "DA neurons")
hDANs_3873.data <- as.data.frame(hDANs_3873@assays$RNA@counts)

hODCs_3873 <- subset(SNc_3873, idents = "Oligodendrocytes")
hODCs_3873.data <- as.data.frame(hODCs_3873@assays$RNA@counts) 

## 4775
SNc_4775 <- subset(nn_SNc, idents = "4775")
SNc_4775 <- NormalizeData(SNc_4775, normalization.method = "LogNormalize", scale.factor = 10000)
SNc_4775 <- FindVariableFeatures(SNc_4775, selection.method = "vst", nfeatures = 3000)
SNc_4775
## An object of class Seurat 
## 38594 features across 25539 samples within 1 assay 
## Active assay: RNA (38594 features, 3000 variable features)

all.genes <- rownames(SNc_4775)
SNc_4775 <- ScaleData(SNc_4775, features = all.genes, verbose = FALSE)
SNc_4775 <- RunPCA(SNc_4775, npcs = 30, verbose = FALSE)

Idents(SNc_4775) <- "batch_ID"
SNc_4775$batch_ID <- SNc_4775@active.ident
table(SNc_4775$batch_ID)

SNc_4775 <- RunHarmony(SNc_4775, group.by.vars = "batch_ID", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(SNc_4775, 'harmony')
harmony_embeddings[1:5, 1:5]

SNc_4775
## An object of class Seurat 
## 38594 features across 25539 samples within 1 assay 
## Active assay: RNA (38594 features, 3000 variable features)
## 2 dimensional reductions calculated: pca, harmony

SNc_4775 <- FindNeighbors(SNc_4775, reduction = "harmony", dims = 1:30)
SNc_4775 <- FindClusters(SNc_4775, resolution = 0.2)
SNc_4775 <- RunUMAP(SNc_4775, reduction = "harmony", dims = 1:30)

DimPlot(SNc_4775, reduction = "umap", label = TRUE)
DimPlot(SNc_4775, reduction = "umap", split.by = "batch_ID", label = TRUE, ncol = 4)

FeaturePlot(SNc_4775, features = c("NR4A2", "TH", "SLC6A3", "SLC18A2"), label = TRUE)
VlnPlot(SNc_4775, features = c("NR4A2", "TH", "SLC6A3", "SLC18A2"), ncol = 1)
FeaturePlot(SNc_4775, features = c("OPALIN", "MOG", "PLP1", "MBP", "CLDN11"), label = TRUE)
FeaturePlot(SNc_4775, features = c("VCAN", "BCAN", "PDGFRA", "OLIG2"), label = TRUE)
FeaturePlot(SNc_4775, features = c("ENO2", "RBFOX3", "SNAP25", "SYT1"), label = TRUE)
VlnPlot(SNc_4775, features = c("ENO2", "RBFOX3", "SNAP25", "SYT1"), ncol = 1)
FeaturePlot(SNc_4775, features = c("AQP4", "GFAP", "SLC1A3", "SLC1A2", "SLC4A4", "NTSR2"), label = TRUE)
FeaturePlot(SNc_4775, features = c("CX3CR1", "P2RY12", "CSF1R"), label = TRUE)
FeaturePlot(SNc_4775, features = c("FLT1", "CLDN5", "PTPRB"), label = TRUE)
FeaturePlot(SNc_4775, features = c("DCN", "COL1A1", "KCNJ8", "ABCC9"), label = TRUE)

SNc_4775 <- RenameIdents(SNc_4775, `0` = "Oligodendrocytes", `1` = "Oligodendrocytes", `2` = "Oligodendrocytes", `3` = "Microglia", `4` = "Astrocytes",
                         `5` = "Non-DA neurons", `6` = "OPCs", `7` = "Non-DA neurons", `8` = "Non-DA neurons", `9` = "Endothelial cells", 
                         `10` = "Non-DA neurons", `11` = "DA neurons", `12` = "Fibroblast", `13` = "Non-DA neurons")

SNc_4775$Major_celltype <- SNc_4775@active.ident
table(SNc_4775$Major_celltype)

hDANs_4775 <- subset(SNc_4775, idents = "DA neurons")
hDANs_4775.data <- as.data.frame(hDANs_4775@assays$RNA@counts)

hODCs_4775 <- subset(SNc_4775, idents = "Oligodendrocytes")
hODCs_4775.data <- as.data.frame(hODCs_4775@assays$RNA@counts) 

## 3887
SNc_3887 <- subset(nn_SNc, idents = "3887")
SNc_3887 <- NormalizeData(SNc_3887, normalization.method = "LogNormalize", scale.factor = 10000)
SNc_3887 <- FindVariableFeatures(SNc_3887, selection.method = "vst", nfeatures = 3000)
SNc_3887
## An object of class Seurat 
## 38594 features across 20374 samples within 1 assay 
## Active assay: RNA (38594 features, 3000 variable features)

all.genes <- rownames(SNc_3887)
SNc_3887 <- ScaleData(SNc_3887, features = all.genes, verbose = FALSE)
SNc_3887 <- RunPCA(SNc_3887, npcs = 30, verbose = FALSE)

Idents(SNc_3887) <- "batch_ID"
SNc_3887$batch_ID <- SNc_3887@active.ident
table(SNc_3887$batch_ID)

SNc_3887 <- RunHarmony(SNc_3887, group.by.vars = "batch_ID", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(SNc_3887, 'harmony')
harmony_embeddings[1:5, 1:5]

SNc_3887
## An object of class Seurat 
## 38594 features across 20374 samples within 1 assay 
## Active assay: RNA (38594 features, 3000 variable features)
## 2 dimensional reductions calculated: pca, harmony

SNc_3887 <- FindNeighbors(SNc_3887, reduction = "harmony", dims = 1:30)
SNc_3887 <- FindClusters(SNc_3887, resolution = 0.2)
SNc_3887 <- RunUMAP(SNc_3887, reduction = "harmony", dims = 1:30)

DimPlot(SNc_3887, reduction = "umap", label = TRUE)
DimPlot(SNc_3887, reduction = "umap", split.by = "batch_ID", label = TRUE, ncol = 4)

FeaturePlot(SNc_3887, features = c("NR4A2", "TH", "SLC6A3", "SLC18A2"), label = TRUE)
VlnPlot(SNc_3887, features = c("NR4A2", "TH", "SLC6A3", "SLC18A2"), ncol = 1)
FeaturePlot(SNc_3887, features = c("OPALIN", "MOG", "PLP1", "MBP", "CLDN11"), label = TRUE)
FeaturePlot(SNc_3887, features = c("VCAN", "BCAN", "PDGFRA", "OLIG2"), label = TRUE)
FeaturePlot(SNc_3887, features = c("ENO2", "RBFOX3", "SNAP25", "SYT1"), label = TRUE)
VlnPlot(SNc_3887, features = c("ENO2", "RBFOX3", "SNAP25", "SYT1"), ncol = 1)
FeaturePlot(SNc_3887, features = c("AQP4", "GFAP", "SLC1A3", "SLC1A2", "SLC4A4", "NTSR2"), label = TRUE)
FeaturePlot(SNc_3887, features = c("CX3CR1", "P2RY12", "CSF1R"), label = TRUE)
FeaturePlot(SNc_3887, features = c("FLT1", "CLDN5", "PTPRB"), label = TRUE)
FeaturePlot(SNc_3887, features = c("DCN", "COL1A1", "KCNJ8", "ABCC9"), label = TRUE)

SNc_3887 <- RenameIdents(SNc_3887, `0` = "Oligodendrocytes", `1` = "Non-DA neurons", `2` = "Non-DA neurons", `3` = "Microglia", `4` = "Oligodendrocytes",
                         `5` = "Non-DA neurons", `6` = "Astrocytes", `7` = "OPCs", `8` = "Non-DA neurons", `9` = "Endothelial cells", 
                         `10` = "Non-DA neurons", `11` = "Fibroblast", `12` = "DA neurons", `13` = "Non-DA neurons")

SNc_3887$Major_celltype <- SNc_3887@active.ident
table(SNc_3887$Major_celltype)

hDANs_3887 <- subset(SNc_3887, idents = "DA neurons")
hDANs_3887.data <- as.data.frame(hDANs_3887@assays$RNA@counts)

hODCs_3887 <- subset(SNc_3887, idents = "Oligodendrocytes")
hODCs_3887.data <- as.data.frame(hODCs_3887@assays$RNA@counts) 

## 2142
SNc_2142 <- subset(nn_SNc, idents = "2142")
SNc_2142 <- NormalizeData(SNc_2142, normalization.method = "LogNormalize", scale.factor = 10000)
SNc_2142 <- FindVariableFeatures(SNc_2142, selection.method = "vst", nfeatures = 3000)
SNc_2142
## An object of class Seurat 
## 38594 features across 11459 samples within 1 assay 
## Active assay: RNA (38594 features, 3000 variable features)

all.genes <- rownames(SNc_2142)
SNc_2142 <- ScaleData(SNc_2142, features = all.genes, verbose = FALSE)
SNc_2142 <- RunPCA(SNc_2142, npcs = 30, verbose = FALSE)

Idents(SNc_2142) <- "batch_ID"
SNc_2142$batch_ID <- SNc_2142@active.ident
table(SNc_2142$batch_ID)

SNc_2142 <- RunHarmony(SNc_2142, group.by.vars = "batch_ID", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(SNc_2142, 'harmony')
harmony_embeddings[1:5, 1:5]

SNc_2142
## An object of class Seurat 
## 38594 features across 11459 samples within 1 assay 
## Active assay: RNA (38594 features, 3000 variable features)
## 2 dimensional reductions calculated: pca, harmony

SNc_2142 <- FindNeighbors(SNc_2142, reduction = "harmony", dims = 1:30)
SNc_2142 <- FindClusters(SNc_2142, resolution = 0.5)
SNc_2142 <- RunUMAP(SNc_2142, reduction = "harmony", dims = 1:30)

DimPlot(SNc_2142, reduction = "umap", label = TRUE)
DimPlot(SNc_2142, reduction = "umap", split.by = "batch_ID", label = TRUE, ncol = 4)

FeaturePlot(SNc_2142, features = c("NR4A2", "TH", "SLC6A3", "SLC18A2"), label = TRUE)
VlnPlot(SNc_2142, features = c("NR4A2", "TH", "SLC6A3", "SLC18A2"), ncol = 1)
FeaturePlot(SNc_2142, features = c("OPALIN", "MOG", "PLP1", "MBP", "CLDN11"), label = TRUE)
FeaturePlot(SNc_2142, features = c("VCAN", "BCAN", "PDGFRA", "OLIG2"), label = TRUE)
FeaturePlot(SNc_2142, features = c("ENO2", "RBFOX3", "SNAP25", "SYT1"), label = TRUE)
VlnPlot(SNc_2142, features = c("ENO2", "RBFOX3", "SNAP25", "SYT1"), ncol = 1)
FeaturePlot(SNc_2142, features = c("AQP4", "GFAP", "SLC1A3", "SLC1A2", "SLC4A4", "NTSR2"), label = TRUE)
FeaturePlot(SNc_2142, features = c("CX3CR1", "P2RY12", "CSF1R"), label = TRUE)
FeaturePlot(SNc_2142, features = c("FLT1", "CLDN5", "PTPRB"), label = TRUE)
FeaturePlot(SNc_2142, features = c("DCN", "COL1A1", "KCNJ8", "ABCC9"), label = TRUE)

SNc_2142 <- RenameIdents(SNc_2142, `0` = "Oligodendrocytes", `1` = "Oligodendrocytes", `2` = "Oligodendrocytes", `3` = "Oligodendrocytes", `4` = "Microglia",
                         `5` = "Astrocytes", `6` = "OPCs", `7` = "Endothelial cells", `8` = "Non-DA neurons", `9` = "Microglia", `10` = "DA neurons")

SNc_2142$Major_celltype <- SNc_2142@active.ident
table(SNc_2142$Major_celltype)

hDANs_2142 <- subset(SNc_2142, idents = "DA neurons")
hDANs_2142.data <- as.data.frame(hDANs_2142@assays$RNA@counts)

hODCs_2142 <- subset(SNc_2142, idents = "Oligodendrocytes")
hODCs_2142.data <- as.data.frame(hODCs_2142@assays$RNA@counts) 

## 1963
SNc_1963 <- subset(nn_SNc, idents = "1963")
SNc_1963 <- NormalizeData(SNc_1963, normalization.method = "LogNormalize", scale.factor = 10000)
SNc_1963 <- FindVariableFeatures(SNc_1963, selection.method = "vst", nfeatures = 3000)
SNc_1963
## An object of class Seurat 
## 38594 features across 15030 samples within 1 assay 
## Active assay: RNA (38594 features, 3000 variable features)

all.genes <- rownames(SNc_1963)
SNc_1963 <- ScaleData(SNc_1963, features = all.genes, verbose = FALSE)
SNc_1963 <- RunPCA(SNc_1963, npcs = 30, verbose = FALSE)

Idents(SNc_1963) <- "batch_ID"
SNc_1963$batch_ID <- SNc_1963@active.ident
table(SNc_1963$batch_ID)

SNc_1963 <- RunHarmony(SNc_1963, group.by.vars = "batch_ID", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(SNc_1963, 'harmony')
harmony_embeddings[1:5, 1:5]

SNc_1963
## An object of class Seurat 
## 38594 features across 15030 samples within 1 assay 
## Active assay: RNA (38594 features, 3000 variable features)
## 2 dimensional reductions calculated: pca, harmony

SNc_1963 <- FindNeighbors(SNc_1963, reduction = "harmony", dims = 1:30)
SNc_1963 <- FindClusters(SNc_1963, resolution = 0.5)
SNc_1963 <- RunUMAP(SNc_1963, reduction = "harmony", dims = 1:30)

DimPlot(SNc_1963, reduction = "umap", label = TRUE)
DimPlot(SNc_1963, reduction = "umap", split.by = "batch_ID", label = TRUE, ncol = 4)

FeaturePlot(SNc_1963, features = c("NR4A2", "TH", "SLC6A3", "SLC18A2"), label = TRUE)
VlnPlot(SNc_1963, features = c("NR4A2", "TH", "SLC6A3", "SLC18A2"), ncol = 1)
FeaturePlot(SNc_1963, features = c("OPALIN", "MOG", "PLP1", "MBP", "CLDN11"), label = TRUE)
FeaturePlot(SNc_1963, features = c("VCAN", "BCAN", "PDGFRA", "OLIG2"), label = TRUE)
FeaturePlot(SNc_1963, features = c("ENO2", "RBFOX3", "SNAP25", "SYT1"), label = TRUE)
VlnPlot(SNc_1963, features = c("ENO2", "RBFOX3", "SNAP25", "SYT1"), ncol = 1)
FeaturePlot(SNc_1963, features = c("AQP4", "GFAP", "SLC1A3", "SLC1A2", "SLC4A4", "NTSR2"), label = TRUE)
FeaturePlot(SNc_1963, features = c("CX3CR1", "P2RY12", "CSF1R"), label = TRUE)
FeaturePlot(SNc_1963, features = c("FLT1", "CLDN5", "PTPRB"), label = TRUE)
FeaturePlot(SNc_1963, features = c("DCN", "COL1A1", "KCNJ8", "ABCC9"), label = TRUE)

SNc_1963 <- RenameIdents(SNc_1963, `0` = "Oligodendrocytes", `1` = "Oligodendrocytes", `2` = "Oligodendrocytes", `3` = "Oligodendrocytes", `4` = "Microglia",
                         `5` = "Astrocytes", `6` = "Astrocytes", `7` = "OPCs", `8` = "Oligodendrocytes", `9` = "Fibroblast", `10` = "Endothelial cells", 
                         `11` = "Neurons", `12` = "Microglia", `13` = "Pericytes", `14` = "Neurons")

SNc_1963$Major_celltype <- SNc_1963@active.ident
table(SNc_1963$Major_celltype)

hDANs_1963 <- subset(SNc_1963, idents = "Neurons")
hDANs_1963 <- subset(hDANs_1963, subset = TH > 0, slot = 'counts')
hDANs_1963 <- subset(hDANs_1963, subset = SLC6A3 > 0, slot = 'counts')
hDANs_1963 <- subset(hDANs_1963, subset = SLC18A2 > 0, slot = 'counts')
hDANs_1963.data <- as.data.frame(hDANs_1963@assays$RNA@counts)

hODCs_1963 <- subset(SNc_1963, idents = "Oligodendrocytes")
hODCs_1963.data <- as.data.frame(hODCs_1963@assays$RNA@counts) 


# hDANs data integration by status
## HC
SNc_hDANs_HC.data <- cbind(hDANs_3298.data, hDANs_3322.data, hDANs_3345.data, hDANs_3346.data, hDANs_3482.data, hDANs_4956.data,
                        hDANs_5610.data, hDANs_6173.data)

SNc_hDANs_HC <- CreateSeuratObject(counts = SNc_hDANs_HC.data, project = "SNc_hDANs_HC", min.cells = 3)

SNc_hDANs_HC@meta.data$batch_ID <- c(hDANs_3298@meta.data$batch_ID, hDANs_3322@meta.data$batch_ID, hDANs_3345@meta.data$batch_ID, hDANs_3346@meta.data$batch_ID, 
                            hDANs_3482@meta.data$batch_ID, hDANs_4956@meta.data$batch_ID, hDANs_5610@meta.data$batch_ID, hDANs_6173@meta.data$batch_ID)

SNc_hDANs_HC@meta.data$Sample_ID <- c(hDANs_3298@meta.data$Sample_ID, hDANs_3322@meta.data$Sample_ID, hDANs_3345@meta.data$Sample_ID, hDANs_3346@meta.data$Sample_ID, 
                                      hDANs_3482@meta.data$Sample_ID, hDANs_4956@meta.data$Sample_ID, hDANs_5610@meta.data$Sample_ID, hDANs_6173@meta.data$Sample_ID)

SNc_hDANs_HC@meta.data$status <- c(hDANs_3298@meta.data$status, hDANs_3322@meta.data$status, hDANs_3345@meta.data$status, hDANs_3346@meta.data$status, 
                                   hDANs_3482@meta.data$status, hDANs_4956@meta.data$status, hDANs_5610@meta.data$status, hDANs_6173@meta.data$status)
SNc_hDANs_HC
## An object of class Seurat 
## 31637 features across 15759 samples within 1 assay 
## Active assay: RNA (31637 features, 0 variable features)

SNc_hDANs_HC[["percent.mt"]] <- PercentageFeatureSet(SNc_hDANs_HC, pattern = "^MT-")
Idents(SNc_hDANs_HC) <- "Sample_ID"
SNc_hDANs_HC$Sample_ID <- SNc_hDANs_HC@active.ident

Idents(SNc_hDANs_HC) <- "status"
SNc_hDANs_HC$status <- SNc_hDANs_HC@active.ident

table(SNc_hDANs_HC$Sample_ID)
table(SNc_hDANs_HC$status)

VlnPlot(SNc_hDANs_HC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

SNc_hDANs_HC <- NormalizeData(SNc_hDANs_HC, normalization.method = "LogNormalize", scale.factor = 10000)
SNc_hDANs_HC <- FindVariableFeatures(SNc_hDANs_HC, selection.method = "vst", nfeatures = 3000)

## PD 
SNc_hDANs_PD.data <- cbind(hDANs_4775.data, hDANs_1963.data, hDANs_2142.data, hDANs_3873.data, hDANs_3887.data, hDANs_4560.data,
                           hDANs_4568.data)

SNc_hDANs_PD <- CreateSeuratObject(counts = SNc_hDANs_PD.data, project = "SNc_hDANs_PD", min.cells = 3)

SNc_hDANs_PD@meta.data$batch_ID <- c(hDANs_4775@meta.data$batch_ID, hDANs_1963@meta.data$batch_ID, hDANs_2142@meta.data$batch_ID, hDANs_3873@meta.data$batch_ID, 
                                     hDANs_3887@meta.data$batch_ID, hDANs_4560@meta.data$batch_ID, hDANs_4568@meta.data$batch_ID)

SNc_hDANs_PD@meta.data$Sample_ID <- c(hDANs_4775@meta.data$Sample_ID, hDANs_1963@meta.data$Sample_ID, hDANs_2142@meta.data$Sample_ID, hDANs_3873@meta.data$Sample_ID, 
                                      hDANs_3887@meta.data$Sample_ID, hDANs_4560@meta.data$Sample_ID, hDANs_4568@meta.data$Sample_ID)

SNc_hDANs_PD@meta.data$status <- c(hDANs_4775@meta.data$status, hDANs_1963@meta.data$status, hDANs_2142@meta.data$status, hDANs_3873@meta.data$status, 
                                   hDANs_3887@meta.data$status, hDANs_4560@meta.data$status, hDANs_4568@meta.data$status)
SNc_hDANs_PD
## An object of class Seurat 
## 28398 features across 2798 samples within 1 assay 
## Active assay: RNA (28398 features, 0 variable features)

SNc_hDANs_PD[["percent.mt"]] <- PercentageFeatureSet(SNc_hDANs_PD, pattern = "^MT-")
Idents(SNc_hDANs_PD) <- "Sample_ID"
SNc_hDANs_PD$Sample_ID <- SNc_hDANs_PD@active.ident

Idents(SNc_hDANs_PD) <- "status"
SNc_hDANs_PD$status <- SNc_hDANs_PD@active.ident

table(SNc_hDANs_PD$Sample_ID)
table(SNc_hDANs_PD$status)

VlnPlot(SNc_hDANs_PD, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

SNc_hDANs_PD <- NormalizeData(SNc_hDANs_PD, normalization.method = "LogNormalize", scale.factor = 10000)
SNc_hDANs_PD <- FindVariableFeatures(SNc_hDANs_PD, selection.method = "vst", nfeatures = 3000)

## Performing integration
SNc_hDANs.list = c(SNc_hDANs_HC, SNc_hDANs_PD)
features <- SelectIntegrationFeatures(object.list = SNc_hDANs.list)

hDANs.anchors <- FindIntegrationAnchors(object.list = SNc_hDANs.list, anchor.features = features)
SNc_hDANs.combined <- IntegrateData(anchorset = hDANs.anchors)

DefaultAssay(SNc_hDANs.combined) <- "integrated"

all.genes <- rownames(SNc_hDANs.combined)
SNc_hDANs.combined <- ScaleData(SNc_hDANs.combined, features = all.genes,verbose = FALSE)
SNc_hDANs.combined <- RunPCA(SNc_hDANs.combined, npcs = 30, verbose = FALSE)

VizDimLoadings(SNc_hDANs.combined, dims = 1:9, reduction = "pca")

SNc_hDANs.combined <- RunHarmony(SNc_hDANs.combined, group.by.vars = "Sample_ID", plot_convergence = TRUE, max.iter.harmony = 100)
harmony_embeddings <- Embeddings(SNc_hDANs.combined, 'harmony')
harmony_embeddings[1:5, 1:5]

SNc_hDANs.combined
## An object of class Seurat 
## 33754 features across 18557 samples within 2 assays 
## Active assay: integrated (2000 features, 2000 variable features)
## 1 other assay present: RNA
## 2 dimensional reductions calculated: pca, harmony

SNc_hDANs.combined <- FindNeighbors(SNc_hDANs.combined, reduction = "harmony", dims = 1:20)
SNc_hDANs.combined <- FindClusters(SNc_hDANs.combined, resolution = 0.2)
SNc_hDANs.combined <- RunUMAP(SNc_hDANs.combined, reduction = "harmony", dims = 1:20)

DimPlot(SNc_hDANs.combined, reduction = "umap", label = TRUE, repel = TRUE)

SNc_hDANs.combined <- RenameIdents(SNc_hDANs.combined, `0` = "hDAN_1", `1` = "hDAN_2", `2` = "hDAN_3", `3` = "hDAN_4", `4` = "hDAN_5", `5` = "hDAN_6",
                                `6` = "hDAN_7", `7` = "hDAN_8")

SNc_hDANs.combined$DAN_subtype <- SNc_hDANs.combined@active.ident
Idents(SNc_hDANs.combined) <- "DAN_subtype"

DimPlot(SNc_hDANs.combined, reduction = "umap", split.by = "status", label = TRUE)

## Cell number changes of different hDAN subtypes under PD conditions   
table(SNc_hDANs.combined$status)
##    HC     PD 
## 15759   2798 

table(Idents(SNc_hDANs.combined), SNc_hDANs.combined$status)
##          HC   PD
## hDAN_1 4700  495
## hDAN_2 4019  348
## hDAN_3 2016  715
## hDAN_4 1682  547
## hDAN_5 1846  149
## hDAN_6  796  326
## hDAN_7  310  187
## hDAN_8  390   31  

prop.table(table(Idents(SNc_hDANs.combined), SNc_hDANs.combined$status), margin = 2)
##                HC         PD
## hDAN_1 0.29824227 0.17691208
## hDAN_2 0.25502887 0.12437455
## hDAN_3 0.12792690 0.25553967
## hDAN_4 0.10673266 0.19549678
## hDAN_5 0.11713941 0.05325232
## hDAN_6 0.05051082 0.11651179
## hDAN_7 0.01967130 0.06683345
## hDAN_8 0.02474776 0.01107934

SNc_hDANs.combined$subtype.status <- paste(Idents(SNc_hDANs.combined), SNc_hDANs.combined$status, sep = "_")

## subset HC hDANs
Idents(SNc_hDANs.combined) <- "status"
SNc_hDANs.combined_HC <- subset(SNc_hDANs.combined, idents = "HC")
SNc_hDANs.combined_PD <- subset(SNc_hDANs.combined, idents = "PD")

Idents(SNc_hDANs.combined_HC) <- "DAN_subtype"
Idents(SNc_hDANs.combined_PD) <- "DAN_subtype"
Idents(SNc_hDANs.combined) <- "DAN_subtype"

DefaultAssay(SNc_hDANs.combined_HC) <- "RNA"
DefaultAssay(SNc_hDANs.combined_PD) <- "RNA"

all.genes <- rownames(SNc_hDANs.combined_HC)
SNc_hDANs.combined_HC <- ScaleData(SNc_hDANs.combined_HC, features = all.genes,verbose = FALSE)

hDANs_subtype.markers <- FindAllMarkers(SNc_hDANs.combined_HC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(hDANs_subtype.markers, file = "hDANs_subtype.markers.txt", col.names = TRUE, sep = "\t", quote = FALSE)

hDANs_subtype.markers_cosg <- cosg(SNc_hDANs.combined_HC, groups = "all", assay = "RNA", slot = "data", mu=1, remove_lowly_expressed = T, expressed_pct = 0.5,
                            n_genes_user = 50)
write.table(hDANs_subtype.markers_cosg$names, file = "hDANs_subtype.markers_cosg.txt", col.names = TRUE, sep = "\t", quote = FALSE)

FeaturePlot(SNc_hDANs.combined_HC, features = c("NR4A2", "TH", "SLC6A3", "SLC18A2", "SOX6", "CALB1", "AGTR1"), label = TRUE)

VlnPlot(SNc_hDANs.combined_HC, features = c("NR4A2", "TH", "SLC6A3", "SLC18A2", "SOX6", "CALB1", "AGTR1", "EYA4", "TUBB2A", "NTNG1", "PARD3B", 
                                           "CDH13", "RELN", "PTPRK", "NEAT1"), stack = TRUE, same.y.lims = TRUE, flip = TRUE, pt.size = 0)


# SNc_hDANs_HC subtype and hODCs CellChat
hODCs_HC.data <- cbind(hODCs_3298.data, hODCs_3322.data, hODCs_3345.data, hODCs_3346.data, hODCs_3482.data, hODCs_4956.data, hODCs_5610.data, hODCs_6173.data)

hODCs_HC <- CreateSeuratObject(counts = hODCs_HC.data, project = "hODCs_HC", min.cells = 3)

hODCs_HC@meta.data$batch_ID <- c(hODCs_3298@meta.data$batch_ID, hODCs_3322@meta.data$batch_ID, hODCs_3345@meta.data$batch_ID,
                                hODCs_3346@meta.data$batch_ID, hODCs_3482@meta.data$batch_ID, hODCs_4956@meta.data$batch_ID,
                                hODCs_5610@meta.data$batch_ID, hODCs_6173@meta.data$batch_ID)

hODCs_HC@meta.data$Sample_ID <- c(hODCs_3298@meta.data$Sample_ID, hODCs_3322@meta.data$Sample_ID, hODCs_3345@meta.data$Sample_ID,
                                      hODCs_3346@meta.data$Sample_ID, hODCs_3482@meta.data$Sample_ID, hODCs_4956@meta.data$Sample_ID,
                                      hODCs_5610@meta.data$Sample_ID, hODCs_6173@meta.data$Sample_ID)

hODCs_HC@meta.data$status <- c(hODCs_3298@meta.data$status, hODCs_3322@meta.data$status, hODCs_3345@meta.data$status,
                                   hODCs_3346@meta.data$status, hODCs_3482@meta.data$status, hODCs_4956@meta.data$status,
                                   hODCs_5610@meta.data$status, hODCs_6173@meta.data$status)
hODCs_HC
## An object of class Seurat 
## 31817 features across 76832 samples within 1 assay 
## Active assay: RNA (31817 features, 0 variable features)

Idents(hODCs_HC) <- "Sample_ID"
hODCs_HC$Sample_ID <- hODCs_HC@active.ident

Idents(hODCs_HC) <- "status"
hODCs_HC$status <- hODCs_HC@active.ident

table(hODCs_HC$Sample_ID)
table(hODCs_HC$status)

hODCs_HC[["percent.mt"]] <- PercentageFeatureSet(hODCs_HC, pattern = "^MT-")
VlnPlot(hODCs_HC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

hODCs_HC <- NormalizeData(hODCs_HC, normalization.method = "LogNormalize", scale.factor = 10000)
hODCs_HC <- FindVariableFeatures(hODCs_HC, selection.method = "vst", nfeatures = 2000)

HC_hDANs_ODCs_CC <- merge(hODCs_HC, y = SNc_hDANs.combined_HC, add.cell.ids = c("ODCs", "hDANs"), project = "HC_hDANs_ODCs_CC",
                    merge.data = TRUE)

HC_hDANs_ODCs_CC <- RenameIdents(HC_hDANs_ODCs_CC, `HC` = "ODCs")

HC_hDANs_ODCs_CC$cell_type <- HC_hDANs_ODCs_CC@active.ident
data.input <- GetAssayData(HC_hDANs_ODCs_CC, assay = "RNA", slot = "data")
identity <- subset(HC_hDANs_ODCs_CC@meta.data, select = "cell_type")

HC_hDANs_ODCs_CC <- createCellChat(object = data.input, meta = identity, group.by = "cell_type")
levels(HC_hDANs_ODCs_CC@idents)
groupSize <- as.numeric(table(HC_hDANs_ODCs_CC@idents))

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB
HC_hDANs_ODCs_CC@DB <- CellChatDB.use

HC_hDANs_ODCs_CC <- subsetData(HC_hDANs_ODCs_CC)
HC_hDANs_ODCs_CC <- identifyOverExpressedGenes(HC_hDANs_ODCs_CC)
HC_hDANs_ODCs_CC <- identifyOverExpressedInteractions(HC_hDANs_ODCs_CC)
HC_hDANs_ODCs_CC <- projectData(HC_hDANs_ODCs_CC, PPI.human)

HC_hDANs_ODCs_CC <- computeCommunProb(HC_hDANs_ODCs_CC, raw.use = TRUE)
HC_hDANs_ODCs_CC <- filterCommunication(HC_hDANs_ODCs_CC, min.cells = 10)

HC_hDANs_ODCs_CC <- computeCommunProbPathway(HC_hDANs_ODCs_CC)
HC_hDANs_ODCs_CC <- aggregateNet(HC_hDANs_ODCs_CC)

par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(HC_hDANs_ODCs_CC@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge = F, title.name = "Number of interactions")

netVisual_circle(HC_hDANs_ODCs_CC@net$weight, vertex.weight = groupSize, weight.scale = T,
                 label.edge = F, title.name = "Interaction strength")

mat_c <- HC_hDANs_ODCs_CC@net$count
mat_w <- HC_hDANs_ODCs_CC@net$weight

par(mfrow = c(3,3), xpd = TRUE)

for (i in 1:nrow(mat_c)) {
  mat_c_1 <- matrix(0, nrow = nrow(mat_c), ncol = ncol(mat_c), dimnames = dimnames(mat_c))
  mat_c_1[i, ] <- mat_c[i, ]
  netVisual_circle(mat_c_1, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat_c), title.name = rownames(mat_c)[i])
}

for (i in 1:nrow(mat_w)) {
  mat_w_1 <- matrix(0, nrow = nrow(mat_w), ncol = ncol(mat_w), dimnames = dimnames(mat_w))
  mat_w_1[i, ] <- mat_w[i, ]
  netVisual_circle(mat_w_1, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat_w), title.name = rownames(mat_w)[i])
}

write.table(mat_c, file = "HC_hDANs_ODCs_CC_count.txt", col.names = TRUE, sep = "\t", quote = FALSE)
write.table(mat_w, file = "PD_hDANs_ODCs_CC_weight.txt", col.names = TRUE, sep = "\t", quote = FALSE)

hDANs_subtype_CC <- read.csv(file = "SP table 3 hDANs subtype and and ODCs CellChat.csv", header = TRUE)

plot(hDANs_subtype_CC$ODCs.outgoing.count, hDANs_subtype_CC$Proportion.change)
cor(hDANs_subtype_CC$ODCs.outgoing.count, hDANs_subtype_CC$Proportion.change)
cor.test(hDANs_subtype_CC$ODCs.outgoing.count, hDANs_subtype_CC$Proportion.change)
## Pearson's product-moment correlation

## data:  hDANs_subtype_CC$ODCs.outgoing.count and hDANs_subtype_CC$Proportion.change
## t = 1.3318, df = 6, p-value = 0.2313
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##   -0.3421846  0.8845885
## sample estimates:
##   cor 
## 0.4776674 

result <- lm(Proportion.change~ODCs.outgoing.count, data = hDANs_subtype_CC)
result
summary(result)
plot(hDANs_subtype_CC$ODCs.outgoing.count, hDANs_subtype_CC$Proportion.change, xlim = c(0, 80), ylim = c(-1, 3),
     xaxs = "i", yaxs = "i", xaxp = c(0, 80, 4), yaxp = c(-1, 3, 4))
abline(result, col = "red", lwd = 2)

plot(hDANs_subtype_CC$ODCs.incoming.count, hDANs_subtype_CC$Proportion.change)
cor(hDANs_subtype_CC$ODCs.incoming.count, hDANs_subtype_CC$Proportion.change)
cor.test(hDANs_subtype_CC$ODCs.incoming.count, hDANs_subtype_CC$Proportion.change)
## Pearson's product-moment correlation

## data:  hDANs_subtype_CC$ODCs.incoming.count and hDANs_subtype_CC$Proportion.change
## t = 2.2256, df = 6, p-value = 0.06768
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##   -0.06121503  0.93437018
## sample estimates:
##   cor 
## 0.6724658 

result <- lm(Proportion.change~ODCs.incoming.count, data = hDANs_subtype_CC)
result
summary(result)
plot(hDANs_subtype_CC$ODCs.incoming.count, hDANs_subtype_CC$Proportion.change, xlim = c(0, 80), ylim = c(-1, 3),
     xaxs = "i", yaxs = "i", xaxp = c(0, 80, 4), yaxp = c(-1, 3, 4))
abline(result, col = "red", lwd = 2)

plot(hDANs_subtype_CC$ODCs.outgoing.weight, hDANs_subtype_CC$Proportion.change)
cor(hDANs_subtype_CC$ODCs.outgoing.weight, hDANs_subtype_CC$Proportion.change)
cor.test(hDANs_subtype_CC$ODCs.outgoing.weight, hDANs_subtype_CC$Proportion.change)
## Pearson's product-moment correlation

## data:  hDANs_subtype_CC$ODCs.outgoing.weight and hDANs_subtype_CC$Proportion.change
## t = 1.499, df = 6, p-value = 0.1845
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##   -0.2889767  0.8967938
## sample estimates:
##   cor 
## 0.5219915  

result <- lm(Proportion.change~ODCs.outgoing.weight, data = hDANs_subtype_CC)
result
summary(result)
plot(hDANs_subtype_CC$ODCs.outgoing.weight, hDANs_subtype_CC$Proportion.change, xlim = c(0, 2.8), ylim = c(-1, 3),
     xaxs = "i", yaxs = "i", xaxp = c(0, 2.8, 4), yaxp = c(-1, 3, 4))
abline(result, col = "red", lwd = 2)

plot(hDANs_subtype_CC$ODCs.incoming.weight, hDANs_subtype_CC$Proportion.change)
cor(hDANs_subtype_CC$ODCs.incoming.weight, hDANs_subtype_CC$Proportion.change)
cor.test(hDANs_subtype_CC$ODCs.incoming.weight, hDANs_subtype_CC$Proportion.change)
## Pearson's product-moment correlation

## data:  hDANs_subtype_CC$ODCs.incoming.weight and hDANs_subtype_CC$Proportion.change
## t = 1.9106, df = 6, p-value = 0.1046
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##   -0.1582027  0.9206846
## sample estimates:
##   cor 
## 0.615035  

result <- lm(Proportion.change~ODCs.incoming.weight, data = hDANs_subtype_CC)
result
summary(result)
plot(hDANs_subtype_CC$ODCs.incoming.weight, hDANs_subtype_CC$Proportion.change, xlim = c(0, 2.8), ylim = c(-1, 3),
     xaxs = "i", yaxs = "i", xaxp = c(0, 2.8, 4), yaxp = c(-1, 3, 4))
abline(result, col = "red", lwd = 2)

HC_hDANs_ODCs_CC <- netAnalysis_computeCentrality(HC_hDANs_ODCs_CC, slot.name = "netP")

write.table(HC_hDANs_ODCs_CC@netP$pathways, file = "pathways.txt", col.names = TRUE, sep = "\t", quote = FALSE)

pathways.show <- "JAM"

netVisual_aggregate(HC_hDANs_ODCs_CC, signaling = pathways.show, layout = "circle")
netVisual_heatmap(HC_hDANs_ODCs_CC, signaling = pathways.show, color.heatmap = "Reds")
netAnalysis_contribution(HC_hDANs_ODCs_CC, signaling = pathways.show)

netVisual_bubble(HC_hDANs_ODCs_CC, sources.use = 1, targets.use = c(2,3,4,5,6,7,8,9), remove.isolate = FALSE)
netVisual_bubble(HC_hDANs_ODCs_CC, sources.use = c(2,3,4,5,6,7,8,9), targets.use = 1, remove.isolate = FALSE)

netVisual_chord_gene(HC_hDANs_ODCs_CC, sources.use = 1, targets.use = c(2,3,4,5,6,7,8,9), signaling = pathways.show, legend.pos.x = 8)
netVisual_chord_gene(HC_hDANs_ODCs_CC, sources.use = c(2,3,4,5,6,7,8,9), targets.use = 1, signaling = pathways.show, legend.pos.x = 8)

netVisual_chord_gene(HC_hDANs_ODCs_CC, sources.use = 1, targets.use = c(2,3,4,5,6,7,8,9), slot.name = "netP", legend.pos.x = 10)
netVisual_chord_gene(HC_hDANs_ODCs_CC, sources.use = c(2,3,4,5,6,7,8,9), targets.use = 1, slot.name = "netP", legend.pos.x = 10)

FeaturePlot(SNc_hDANs.combined_HC, features = c("JAM2", "JAM3"), label = TRUE)
VlnPlot(SNc_hDANs.combined_HC, features = c("JAM2", "JAM3"), assay = "RNA", pt.size = 1, ncol = 1)


# SNc_hDANs_PD subtype and hODCs CellChat
hODCs_PD.data <- cbind(hODCs_4775.data, hODCs_1963.data, hODCs_2142.data, hODCs_3873.data, hODCs_3887.data, hODCs_4560.data, hODCs_4568.data)

hODCs_PD <- CreateSeuratObject(counts = hODCs_PD.data, project = "hODCs_PD", min.cells = 3)

hODCs_PD@meta.data$batch_ID <- c(hODCs_4775@meta.data$batch_ID, hODCs_1963@meta.data$batch_ID, hODCs_2142@meta.data$batch_ID,
                                   hODCs_3873@meta.data$batch_ID, hODCs_3887@meta.data$batch_ID, hODCs_4560@meta.data$batch_ID,
                                   hODCs_4568@meta.data$batch_ID)

hODCs_PD@meta.data$Sample_ID <- c(hODCs_4775@meta.data$Sample_ID, hODCs_1963@meta.data$Sample_ID, hODCs_2142@meta.data$Sample_ID,
                                    hODCs_3873@meta.data$Sample_ID, hODCs_3887@meta.data$Sample_ID, hODCs_4560@meta.data$Sample_ID,
                                    hODCs_4568@meta.data$Sample_ID)

hODCs_PD@meta.data$status <- c(hODCs_4775@meta.data$status, hODCs_1963@meta.data$status, hODCs_2142@meta.data$status,
                                 hODCs_3873@meta.data$status, hODCs_3887@meta.data$status, hODCs_4560@meta.data$status,
                                 hODCs_4568@meta.data$status)
hODCs_PD
## An object of class Seurat 
## 31514 features across 80686 samples within 1 assay 
## Active assay: RNA (31514 features, 0 variable features)

Idents(hODCs_PD) <- "Sample_ID"
hODCs_PD$Sample_ID <- hODCs_PD@active.ident

Idents(hODCs_PD) <- "status"
hODCs_PD$status <- hODCs_PD@active.ident

table(hODCs_PD$Sample_ID)
table(hODCs_PD$status)

hODCs_PD[["percent.mt"]] <- PercentageFeatureSet(hODCs_PD, pattern = "^MT-")
VlnPlot(hODCs_PD, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

hODCs_PD <- NormalizeData(hODCs_PD, normalization.method = "LogNormalize", scale.factor = 10000)
hODCs_PD <- FindVariableFeatures(hODCs_PD, selection.method = "vst", nfeatures = 2000)

PD_hDANs_ODCs_CC <- merge(hODCs_PD, y = SNc_hDANs.combined_PD, add.cell.ids = c("ODCs", "hDANs"), project = "PD_hDANs_ODCs_CC",
                         merge.data = TRUE)

PD_hDANs_ODCs_CC <- RenameIdents(PD_hDANs_ODCs_CC, `PD` = "ODCs")

PD_hDANs_ODCs_CC$cell_type <- PD_hDANs_ODCs_CC@active.ident
data.input <- GetAssayData(PD_hDANs_ODCs_CC, assay = "RNA", slot = "data")
identity <- subset(PD_hDANs_ODCs_CC@meta.data, select = "cell_type")

PD_hDANs_ODCs_CC <- createCellChat(object = data.input, meta = identity, group.by = "cell_type")
levels(PD_hDANs_ODCs_CC@idents)
groupSize <- as.numeric(table(PD_hDANs_ODCs_CC@idents))

PD_hDANs_ODCs_CC@DB <- CellChatDB.use

PD_hDANs_ODCs_CC <- subsetData(PD_hDANs_ODCs_CC)
PD_hDANs_ODCs_CC <- identifyOverExpressedGenes(PD_hDANs_ODCs_CC)
PD_hDANs_ODCs_CC <- identifyOverExpressedInteractions(PD_hDANs_ODCs_CC)
PD_hDANs_ODCs_CC <- projectData(PD_hDANs_ODCs_CC, PPI.human)

PD_hDANs_ODCs_CC <- computeCommunProb(PD_hDANs_ODCs_CC, raw.use = TRUE)
PD_hDANs_ODCs_CC <- filterCommunication(PD_hDANs_ODCs_CC, min.cells = 10)

PD_hDANs_ODCs_CC <- computeCommunProbPathway(PD_hDANs_ODCs_CC)
PD_hDANs_ODCs_CC <- aggregateNet(PD_hDANs_ODCs_CC)

PD_hDANs_ODCs_CC <- netAnalysis_computeCentrality(PD_hDANs_ODCs_CC, slot.name = "netP")


# CellChat integration
CellChat.list <- list(HC = HC_hDANs_ODCs_CC, PD = PD_hDANs_ODCs_CC)
I_CellChat <- mergeCellChat(CellChat.list, add.names = names(CellChat.list), cell.prefix = TRUE)

I_CellChat
## An object of class CellChat created from a merged object with multiple datasets 
## 928 signaling genes.
## 176075 cells. 
## CellChat analysis of single cell RNA-seq data! 

gg1 <- compareInteractions(I_CellChat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(I_CellChat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

par(mfrow = c(1,2), xpd = TRUE)
netVisual_diffInteraction(I_CellChat, weight.scale = T)
netVisual_diffInteraction(I_CellChat, weight.scale = T, measure = "weight")


# JAM3 expression 
SNc_hDANs.combined@meta.data$JAM3 <- SNc_hDANs.combined@assays$RNA@counts["JAM3",]

DefaultAssay(SNc_hDANs.combined) <- "RNA"
SNc_hDANs.combined <- ScaleData(SNc_hDANs.combined, features = rownames(SNc_hDANs.combined), verbose = FALSE)

VlnPlot(SNc_hDANs.combined, features = "rna_JAM3", split.by  = "status", split.plot = TRUE, assay = "RNA")

JAM3_pos_ratio <- as.data.frame(prop.table(table(SNc_hDANs.combined$JAM3 > 0 , SNc_hDANs.combined$subtype.status), margin = 2))
JAM3_pos_ratio <- subset(JAM3_pos_ratio, JAM3_pos_ratio$Var1 == "TRUE")
JAM3_pos_ratio <- JAM3_pos_ratio[,-1]
JAM3_pos_ratio$status <- rep(c("HC", "PD"), 8)
JAM3_pos_ratio$subtype <- c(rep("hDAN_1", 2), rep("hDAN_2", 2), rep("hDAN_3", 2), rep("hDAN_4", 2), rep("hDAN_5", 2), rep("hDAN_6", 2),
                            rep("hDAN_7", 2), rep("hDAN_8", 2))

JAM3_pos_ratio$status <- factor(x = JAM3_pos_ratio$status, levels = c("HC", "PD"))
JAM3_pos_ratio$subtype <- factor(x = JAM3_pos_ratio$subtype, levels = c("hDAN_1", "hDAN_2", "hDAN_3", "hDAN_4", "hDAN_5", "hDAN_6",
                                                                        "hDAN_7", "hDAN_8"))

ggplot() + geom_bar(data = JAM3_pos_ratio, aes(x = subtype, y = Freq, fill = status), position = position_dodge2(padding = 0.3), stat = "identity") + 
  theme(axis.title = element_blank(), axis.text.x = element_text(size = 20), axis.ticks.x = element_blank(),
        axis.ticks.length.y = unit(.25, "cm"), axis.line = element_line(color = "black", size = 0.4), 
        axis.text.y = element_text(size = 14), panel.background = element_rect(fill = "white"), legend.position = "top", 
        legend.title = element_blank(), legend.key.height = unit(.3, "cm"), legend.key.width = unit(.8, "cm"), 
        legend.direction = "vertical", legend.spacing.x = unit(.4, "cm"), legend.text = element_text(size = 12)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,1, by = 0.2), limits = c(0,1))

hDANs_subtype_JAM3_averExp <- AverageExpression(SNc_hDANs.combined_HC, slot = "counts", features = "JAM3")$RNA
hDANs_subtype_CC$JAM3_averExp <- hDANs_subtype_JAM3_averExp[1, match(hDANs_subtype_CC$X, colnames(hDANs_subtype_JAM3_averExp))]

plot(hDANs_subtype_CC$JAM3_averExp, hDANs_subtype_CC$Proportion.change)
cor(hDANs_subtype_CC$JAM3_averExp, hDANs_subtype_CC$Proportion.change)
cor.test(hDANs_subtype_CC$JAM3_averExp, hDANs_subtype_CC$Proportion.change)
## Pearson's product-moment correlation

## data:  hDANs_subtype_CC$JAM3_averExp and hDANs_subtype_CC$Proportion.change
## t = 1.2584, df = 6, p-value = 0.255
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##   -0.3653603  0.8786892
## sample estimates:
##   cor 
## 0.4569617 

result <- lm(Proportion.change~JAM3_averExp, data = hDANs_subtype_CC)
result
summary(result)

plot(hDANs_subtype_CC$JAM3_averExp, hDANs_subtype_CC$Proportion.change, xlim = c(0, 6), ylim = c(-1, 3),
     xaxs = "i", yaxs = "i", xaxp = c(0, 6, 4), yaxp = c(-1, 3, 4))
abline(result, col = "red", lwd = 2)

JAM3_pos_ratio <- subset(JAM3_pos_ratio, JAM3_pos_ratio$status == "HC")
hDANs_subtype_CC$JAM3_ExpRatio <- JAM3_pos_ratio[match(hDANs_subtype_CC$X, JAM3_pos_ratio$subtype), 2]

plot(hDANs_subtype_CC$JAM3_ExpRatio, hDANs_subtype_CC$Proportion.change)
cor(hDANs_subtype_CC$JAM3_ExpRatio, hDANs_subtype_CC$Proportion.change)
cor.test(hDANs_subtype_CC$JAM3_ExpRatio, hDANs_subtype_CC$Proportion.change)
## Pearson's product-moment correlation

## data:  hDANs_subtype_CC$JAM3_ExpRatio and hDANs_subtype_CC$Proportion.change
## t = 0.61921, df = 6, p-value = 0.5585
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##   -0.5555335  0.8098853
## sample estimates:
##   cor 
## 0.2450812 

result <- lm(Proportion.change~JAM3_ExpRatio, data = hDANs_subtype_CC)
result
summary(result)

plot(hDANs_subtype_CC$JAM3_ExpRatio, hDANs_subtype_CC$Proportion.change, xlim = c(0, 1), ylim = c(-1, 3),
     xaxs = "i", yaxs = "i", xaxp = c(0, 1, 4), yaxp = c(-1, 3, 4))
abline(result, col = "red", lwd = 2)


# JAM2 expression 
SNc_hDANs.combined@meta.data$JAM2 <- SNc_hDANs.combined@assays$RNA@counts["JAM2",]

VlnPlot(SNc_hDANs.combined, features = "rna_JAM2", split.by  = "status", split.plot = TRUE, assay = "RNA")

JAM2_pos_ratio <- as.data.frame(prop.table(table(SNc_hDANs.combined$JAM2 > 0 , SNc_hDANs.combined$subtype.status), margin = 2))
JAM2_pos_ratio <- subset(JAM2_pos_ratio, JAM2_pos_ratio$Var1 == "TRUE")
JAM2_pos_ratio <- JAM2_pos_ratio[,-1]
JAM2_pos_ratio$status <- rep(c("HC", "PD"), 8)
JAM2_pos_ratio$subtype <- c(rep("hDAN_1", 2), rep("hDAN_2", 2), rep("hDAN_3", 2), rep("hDAN_4", 2), rep("hDAN_5", 2), rep("hDAN_6", 2),
                            rep("hDAN_7", 2), rep("hDAN_8", 2))

JAM2_pos_ratio$status <- factor(x = JAM2_pos_ratio$status, levels = c("HC", "PD"))
JAM2_pos_ratio$subtype <- factor(x = JAM2_pos_ratio$subtype, levels = c("hDAN_1", "hDAN_2", "hDAN_3", "hDAN_4", "hDAN_5", "hDAN_6",
                                                                        "hDAN_7", "hDAN_8"))

ggplot() + geom_bar(data = JAM2_pos_ratio, aes(x = subtype, y = Freq, fill = status), position = position_dodge2(padding = 0.3), stat = "identity") + 
  theme(axis.title = element_blank(), axis.text.x = element_text(size = 20), axis.ticks.x = element_blank(),
        axis.ticks.length.y = unit(.25, "cm"), axis.line = element_line(color = "black", size = 0.4), 
        axis.text.y = element_text(size = 14), panel.background = element_rect(fill = "white"), legend.position = "top", 
        legend.title = element_blank(), legend.key.height = unit(.3, "cm"), legend.key.width = unit(.8, "cm"), 
        legend.direction = "vertical", legend.spacing.x = unit(.4, "cm"), legend.text = element_text(size = 12)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,0.4, by = 0.08), limits = c(0,0.4))

hDANs_subtype_JAM2_averExp <- AverageExpression(SNc_hDANs.combined_HC, slot = "counts", features = "JAM2")$RNA
hDANs_subtype_CC$JAM2_averExp <- hDANs_subtype_JAM2_averExp[1, match(hDANs_subtype_CC$X, colnames(hDANs_subtype_JAM2_averExp))]

plot(hDANs_subtype_CC$JAM2_averExp, hDANs_subtype_CC$Proportion.change)
cor(hDANs_subtype_CC$JAM2_averExp, hDANs_subtype_CC$Proportion.change)
cor.test(hDANs_subtype_CC$JAM2_averExp, hDANs_subtype_CC$Proportion.change)
## Pearson's product-moment correlation

## data:  hDANs_subtype_CC$JAM2_averExp and hDANs_subtype_CC$Proportion.change
## t = 3.4697, df = 6, p-value = 0.01331
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##   0.2645658 0.9656884
## sample estimates:
##   cor 
## 0.8169355  

result <- lm(Proportion.change~JAM2_averExp, data = hDANs_subtype_CC)
result
summary(result)

plot(hDANs_subtype_CC$JAM2_averExp, hDANs_subtype_CC$Proportion.change, xlim = c(0, 0.48), ylim = c(-1, 3),
     xaxs = "i", yaxs = "i", xaxp = c(0, 0.48, 4), yaxp = c(-1, 3, 4))
abline(result, col = "red", lwd = 2)

JAM2_pos_ratio <- subset(JAM2_pos_ratio, JAM3_pos_ratio$status == "HC")
hDANs_subtype_CC$JAM2_ExpRatio <- JAM2_pos_ratio[match(hDANs_subtype_CC$X, JAM2_pos_ratio$subtype), 2]

plot(hDANs_subtype_CC$JAM2_ExpRatio, hDANs_subtype_CC$Proportion.change)
cor(hDANs_subtype_CC$JAM2_ExpRatio, hDANs_subtype_CC$Proportion.change)
cor.test(hDANs_subtype_CC$JAM2_ExpRatio, hDANs_subtype_CC$Proportion.change)
## Pearson's product-moment correlation

## data:  hDANs_subtype_CC$JAM2_ExpRatio and hDANs_subtype_CC$Proportion.change
## t = 3.7328, df = 6, p-value = 0.009705
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##   0.3198006 0.9695340
## sample estimates:
##   cor 
## 0.8360628 

result <- lm(Proportion.change~JAM2_ExpRatio, data = hDANs_subtype_CC)
result
summary(result)

plot(hDANs_subtype_CC$JAM2_ExpRatio, hDANs_subtype_CC$Proportion.change, xlim = c(0, 0.28), ylim = c(-1, 3),
     xaxs = "i", yaxs = "i", xaxp = c(0, 0.28, 4), yaxp = c(-1, 3, 4))
abline(result, col = "red", lwd = 2)


# VGluT2 expression and glutamate signaling 
SNc_hDANs.combined@meta.data$SLC17A6 <- SNc_hDANs.combined@assays$RNA@counts["SLC17A6",]

SLC17A6_pos_ratio <- as.data.frame(prop.table(table(SNc_hDANs.combined$SLC17A6 > 0 , SNc_hDANs.combined$subtype.status), margin = 2))
SLC17A6_pos_ratio <- subset(SLC17A6_pos_ratio, SLC17A6_pos_ratio$Var1 == "TRUE")
SLC17A6_pos_ratio <- SLC17A6_pos_ratio[,-1]

SLC17A6_pos_ratio$status <- rep(c("HC", "PD"), 8)
SLC17A6_pos_ratio$subtype <- c(rep("hDAN_1", 2), rep("hDAN_2", 2), rep("hDAN_3", 2), rep("hDAN_4", 2), rep("hDAN_5", 2), rep("hDAN_6", 2),
                               rep("hDAN_7", 2), rep("hDAN_8", 2))

SLC17A6_pos_ratio$status <- factor(x = SLC17A6_pos_ratio$status, levels = c("HC", "PD"))
SLC17A6_pos_ratio$subtype <- factor(x = SLC17A6_pos_ratio$subtype, levels = c("hDAN_1", "hDAN_2", "hDAN_3", "hDAN_4", "hDAN_5", "hDAN_6",
                                                                              "hDAN_7", "hDAN_8"))

ggplot() + geom_bar(data = SLC17A6_pos_ratio, aes(x = subtype, y = Freq, fill = status), position = position_dodge2(padding = 0.3), stat = "identity") + 
  theme(axis.title = element_blank(), axis.text.x = element_text(size = 20), axis.ticks.x = element_blank(),
        axis.ticks.length.y = unit(.25, "cm"), axis.line = element_line(color = "black", size = 0.4), 
        axis.text.y = element_text(size = 14), panel.background = element_rect(fill = "white"), legend.position = "top", 
        legend.title = element_blank(), legend.key.height = unit(.3, "cm"), legend.key.width = unit(.8, "cm"), 
        legend.direction = "vertical", legend.spacing.x = unit(.4, "cm"), legend.text = element_text(size = 12)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,0.6, by = 0.15), limits = c(0,0.6))

hDANs_subtype_SLC17A6_averExp <- AverageExpression(SNc_hDANs.combined_HC, slot = "counts", features = "SLC17A6")$RNA
hDANs_subtype_CC$SLC17A6_averExp <- hDANs_subtype_SLC17A6_averExp[1, match(hDANs_subtype_CC$X, colnames(hDANs_subtype_SLC17A6_averExp))]

plot(hDANs_subtype_CC$SLC17A6_averExp, hDANs_subtype_CC$Proportion.change)
cor(hDANs_subtype_CC$SLC17A6_averExp, hDANs_subtype_CC$Proportion.change)
cor.test(hDANs_subtype_CC$SLC17A6_averExp, hDANs_subtype_CC$Proportion.change)
## Pearson's product-moment correlation

## data:  hDANs_subtype_CC$SLC17A6_averExp and hDANs_subtype_CC$Proportion.change
## t = 3.8235, df = 6, p-value = 0.008725
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
## 0.3377653 0.9707193
## sample estimates:
##  cor 
## 0.8420246  

result <- lm(Proportion.change~SLC17A6_averExp, data = hDANs_subtype_CC)
result
summary(result)
plot(hDANs_subtype_CC$SLC17A6_averExp, hDANs_subtype_CC$Proportion.change, xlim = c(0, 1), ylim = c(-1, 3),
     xaxs = "i", yaxs = "i", xaxp = c(0, 1, 4), yaxp = c(-1, 3, 4))
abline(result, col = "red", lwd = 2)


# hdWGCNA
theme_set(theme_cowplot())

DefaultAssay(SNc_hDANs.combined_HC) <- "integrated"

HC_DANs <- SetupForWGCNA(SNc_hDANs.combined_HC, gene_select = "variable", wgcna_name = "HC_DANs")
length(HC_DANs@misc$HC_DANs$wgcna_genes)

DefaultAssay(HC_DANs) <- "RNA"

HC_DANs <- MetacellsByGroups(seurat_obj = HC_DANs, group.by = "DAN_subtype", reduction = "pca", k = 20, min_cells = 100, max_shared = 15,
                             ident.group = "DAN_subtype")
HC_DANs <- NormalizeMetacells(HC_DANs)

HC_DANs <- SetDatExpr(HC_DANs, group_name = c("hDAN_3", "hDAN_4", "hDAN_6", "hDAN_7"), group.by = "DAN_subtype", assay = "RNA", slot = "data")
HC_DANs <- TestSoftPowers(HC_DANs, networkType = "signed")
plot_list <- PlotSoftPowers(HC_DANs)
wrap_plots(plot_list, ncol = 2)

power_table <- GetPowerTable(HC_DANs)
head(power_table)

## construct co-expression network
HC_DANs <- ConstructNetwork(HC_DANs, soft_power = 5, setDatExpr = FALSE, tom_name = "hDAN_3467")
PlotDendrogram(HC_DANs, main = "hDAN_3467 hdWGCNA Dendrogram")

HC_DANs@misc$HC_DANs$wgcna_modules %>% head
write.table(HC_DANs@misc$HC_DANs$wgcna_modules, file = "hDAN_3467_modules.txt", col.names = TRUE, sep = "\t", quote = FALSE)

table(HC_DANs@misc$HC_DANs$wgcna_modules$module)
## grey      blue       red     brown turquoise    yellow     green 
## 803       339        54       143       545        60        56 

Tom <- GetTOM(HC_DANs)

## Compute harmonized module eigengenes
HC_DANs <- ScaleData(HC_DANs, features =  VariableFeatures(HC_DANs))
HC_DANs <- ModuleEigengenes(HC_DANs, npcs = 30, reduction.use = "pca", group.by.vars = "Sample_ID")

## harmonized module eigengenes
hMEs <- GetMEs(HC_DANs)
head(hMEs)
MEs <- GetMEs(HC_DANs, harmonized = FALSE)
head(MEs)

## Compute module connectivity
HC_DANs <- ModuleConnectivity(HC_DANs, group.by = "DAN_subtype", group_name = c("hDAN_3", "hDAN_4", "hDAN_6", "hDAN_7"))
HC_DANs <- ResetModuleNames(HC_DANs, new_name = "hDAN_3467-M")

PlotKMEs(HC_DANs, ncol = 5)

modules <- GetModules(HC_DANs)

hub_df <- GetHubGenes(HC_DANs, n_hubs = 10)
head(hub_df)
write.table(hub_df, file = "hub_genes.txt", col.names = TRUE, sep = "\t", quote = FALSE)

HC_DANs <- ModuleExprScore(HC_DANs, n_genes = 25, method = "Seurat")

## Visualization
plot_list <- ModuleFeaturePlot(HC_DANs, features = "hMEs", order = TRUE)
wrap_plots(plot_list, ncol = 3) 

plot_list <- ModuleFeaturePlot(HC_DANs, features = "scores", order = "shuffle", ucell = TRUE)
wrap_plots(plot_list, ncol = 3)

MEs <- GetMEs(HC_DANs, harmonized = TRUE)
mods <- colnames(MEs); mods <- mods[mods != "grey"]

HC_DANs@meta.data <- cbind(HC_DANs@meta.data, MEs)
P <- DotPlot(HC_DANs, features = mods, group.by = "DAN_subtype")
P <- P + coord_flip() + RotatedAxis() + scale_color_gradient2(high = "red", mid = "grey95", low = "blue")
P

theme_set(theme_cowplot())
ModuleNetworkPlot(HC_DANs)

HubGeneNetworkPlot(HC_DANs, n_hubs = 5, n_other = 5, edge_prop = 0.75, mods = 'all', edge.alpha = 0.6)


# resistant DAN makers GO enrichment 
resis_hDANs.markers<-FindMarkers(SNc_hDANs.combined_HC, ident.1 = c("hDAN_3", "hDAN_4", "hDAN_6", "hDAN_7"),
                                 ident.2 = c("hDAN_1", "hDAN_2", "hDAN_5", "hDAN_8"), verbose = FALSE, only.pos = TRUE)

resis_hDANs.markers <- subset(resis_hDANs.markers, p_val_adj < 0.05)
resis_hDANs.markers <- as.character(unique(rownames(resis_hDANs.markers)))

resis_hDANs.markers.df <- bitr(resis_hDANs.markers, fromType = "SYMBOL", toType = c("ENTREZID", "ENSEMBL"), OrgDb = org.Hs.eg.db)
resis_hDANs.markers <- unique(resis_hDANs.markers.df$ENTREZID)

resis_hDANs.markers_ego <- enrichGO(gene = resis_hDANs.markers, OrgDb = org.Hs.eg.db, ont = "MF", readable = TRUE)
resis_hDANs.markers_ego <- enrichGO(gene = resis_hDANs.markers, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)
resis_hDANs.markers_ego <- enrichGO(gene = resis_hDANs.markers, OrgDb = org.Hs.eg.db, ont = "CC", readable = TRUE)

resis_hDANs.markers_ego <- clusterProfiler::simplify(resis_hDANs.markers_ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
dotplot(resis_hDANs.markers_ego, showCategory = 10, title = "resis_hDANs.markers_egoCC")

resis_hDANs.markers_KEGG <- enrichKEGG(gene = resis_hDANs.markers, keyType = "kegg", organism = "hsa", pAdjustMethod = "fdr", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
dotplot(resis_hDANs.markers_KEGG, title = "resis_hDANs.markers_KEGG")


# Glutamate receptors expression in ODCs
hODCs.list <- c(hODCs_HC, hODCs_PD)
features <- SelectIntegrationFeatures(object.list = hODCs.list)
hODCs.anchors <- FindIntegrationAnchors(object.list = hODCs.list, anchor.features = features)
hODCs.combined <- IntegrateData(anchorset = hODCs.anchors)

DefaultAssay(hODCs.combined) <- "integrated"

all.genes <- rownames(hODCs.combined)
hODCs.combined <- ScaleData(hODCs.combined, features = all.genes,verbose = FALSE)
hODCs.combined <- RunPCA(hODCs.combined, npcs = 30, verbose = FALSE)

hODCs.combined <- FindNeighbors(hODCs.combined, reduction = "pca", dims = 1:30)
hODCs.combined <- FindClusters(hODCs.combined, resolution = 0.1)
hODCs.combined <- RunUMAP(hODCs.combined, reduction = "pca", dims = 1:30)

P1 <- DimPlot(hODCs.combined, reduction = "umap", label = TRUE, repel = TRUE)
P2 <- DimPlot(hODCs.combined, reduction = "umap", group.by = "status")
CombinePlots(plots = list(P1, P2), ncol = 1)

hODCs.combined <- RenameIdents(hODCs.combined, `0` = "hODCs_1", `1` = "hODCs_2", `2` = "hODCs_3", `3` = "hODCs_4", `4` = "hODCs_5", `5` = "hODCs_6", `6` = "hODCs_7")
hODCs.combined$ODCs_subtype <- hODCs.combined@active.ident
Idents(hODCs.combined) <- "ODCs_subtype"
DimPlot(hODCs.combined, reduction = "umap", label = TRUE)

DefaultAssay(hODCs.combined) <- "RNA"
all.genes <- rownames(hODCs.combined)
hODCs.combined <- ScaleData(hODCs.combined, features = all.genes,verbose = FALSE)

hODCs.combined@meta.data$status <- factor(x = hODCs.combined@meta.data$status, levels = c("HC", "PD"))

DotPlot(hODCs.combined, features = c("GRM1", "GRM2", "GRM3", "GRM4", "GRM5", "GRM6", "GRM7", "GRM8", "GRIA1", "GRIA2", 
                                       "GRIA3", "GRIA4", "GRID1", "GRID2", "GRIK1", "GRIK2", "GRIK3", "GRIK4", "GRIK5", 
                                       "GRIN1", "GRIN2A", "GRIN2B", "GRIN2D", "GRIN3A", "GRIN3B", "GRINA"), 
        group.by = "status", scale = FALSE) + RotatedAxis()

table(hODCs.combined$status)
## HC    PD 
## 76832 80686 

table(Idents(hODCs.combined), hODCs.combined$status)
##           HC       PD
## hODCs_1 25431    20144
## hODCs_2 26614    12913
## hODCs_3 17960    13457
## hODCs_4  6698    8416
## hODCs_5     2    10705
## hODCs_6    76    8467
## hODCs_7    51    6584

prop.table(table(Idents(hODCs.combined), hODCs.combined$status), margin = 2)
##                  HC           PD
## hODCs_1 3.309949e-01 2.496592e-01
## hODCs_2 3.463921e-01 1.600402e-01
## hODCs_3 2.337568e-01 1.667823e-01
## hODCs_4 8.717722e-02 1.043056e-01
## hODCs_5 2.603082e-05 1.326748e-01
## hODCs_6 9.891712e-04 1.049377e-01
## hODCs_7 6.637859e-04 8.160028e-02




# WD ^ ^ written by Shuxuan Lyu 2024/04/23



















