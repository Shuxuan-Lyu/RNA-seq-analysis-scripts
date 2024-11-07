library(Seurat)
library(SeuratData)
library(DoubletFinder)
library(cowplot)
library(patchwork)
library(ggplot2)
library(ggsci)
library(dplyr)
library(psych)
library(tidyverse)
library(ComplexHeatmap)
library(stringr)
library(RColorBrewer)
library(clusterProfiler)
library(org.Mm.eg.db)
library(DOSE)
library(GOSemSim)
library(enrichplot)
library(COSG)
library(CellChat)
library(plot1cell)
library(ggrepel)
library(WGCNA)
library(hdWGCNA)
library(patchwork)


rm(list = ls())
options(stringsAsFactors = F)


# loading Sham and MPTP dataset
Sham.data <- Read10X(data.dir = "data/Sham.matrix/")
MPTP.data <- Read10X(data.dir = "data/MPTP.matrix/")


# Remove doublets
## Sham
Sham_mb <- CreateSeuratObject(counts = Sham.data, project = "Sham mouse midbrain", min.cells = 3)
Sham_mb@meta.data$status <- c(rep("Sham", ncol(Sham.data)))
Sham_mb
## An object of class Seurat 
## 21769 features across 10680 samples within 1 assay 
## Active assay: RNA (21769 features, 0 variable features)

Sham_mb <- NormalizeData(Sham_mb, normalization.method = "LogNormalize", scale.factor = 10000)
Sham_mb <- FindVariableFeatures(Sham_mb, selection.method = "vst", nfeatures = 3000)

all.genes <- rownames(Sham_mb)
Sham_mb <- ScaleData(Sham_mb, features = all.genes, verbose = FALSE)
Sham_mb <- RunPCA(Sham_mb, npcs = 30, verbose = FALSE)
VizDimLoadings(Sham_mb, dims = 1:20, reduction = "pca")
Sham_mb <- JackStraw(Sham_mb, num.replicate = 200)
Sham_mb <- ScoreJackStraw(Sham_mb, dims = 1:20)
JackStrawPlot(Sham_mb, dims = 1:20)
ElbowPlot(Sham_mb)

Sham_mb <- FindNeighbors(Sham_mb, reduction = "pca", dims = 1:30)
Sham_mb <- FindClusters(Sham_mb, resolution = 0.3)
Sham_mb <- RunUMAP(Sham_mb, dims = 1:30)
Sham_mb_P1 <- DimPlot(Sham_mb, reduction = "umap", label = TRUE)
Sham_mb
## An object of class Seurat 
## 21769 features across 10680 samples within 1 assay 
## Active assay: RNA (21769 features, 3000 variable features)
## 2 dimensional reductions calculated: pca, umap

sweep.res.list_Sham_mb <- paramSweep_v3(Sham_mb, PCs = 1:30, sct = FALSE)
sweep.stats_Sham_mb <- summarizeSweep(sweep.res.list_Sham_mb, GT = FALSE)
bcmvn_Sham_mb <- find.pK(sweep.stats_Sham_mb)
opt_pK_Sham_mb <- as.numeric(as.vector(bcmvn_Sham_mb$pK[which.max(bcmvn_Sham_mb$BCmetric)]))
print(opt_pK_Sham_mb)
## [1] 0.005

annotations_Sham_mb <- Sham_mb@meta.data$seurat_clusters
homotypic.prop_Sham_mb <- modelHomotypic(annotations_Sham_mb)
nExp_poi_Sham_mb <- round(0.076*nrow(Sham_mb@meta.data))
nExp_poi.adj_Sham_mb <- round(nExp_poi_Sham_mb*(1-homotypic.prop_Sham_mb))

Sham_mb <- doubletFinder_v3(Sham_mb, PCs = 1:30, pN = 0.25, pK = opt_pK_Sham_mb, nExp = nExp_poi.adj_Sham_mb, 
                            reuse.pANN = FALSE, sct = FALSE)  
Sham_mb_P2 <- DimPlot(Sham_mb, reduction = "umap", group.by = "DF.classifications_0.25_0.005_715")
CombinePlots(plots = list(Sham_mb_P1, Sham_mb_P2), ncol = 1)

table(Sham_mb@meta.data$DF.classifications_0.25_0.005_715)
## Doublet Singlet 
## 715     9965

Idents(Sham_mb) <- "DF.classifications_0.25_0.005_715"
table(Sham_mb@active.ident)

Sham_mb_RemDoub <- subset(Sham_mb, idents = "Singlet")
Sham_mb_RemDoub
## An object of class Seurat 
## 21769 features across 9965 samples within 1 assay 
## Active assay: RNA (21769 features, 3000 variable features)
## 2 dimensional reductions calculated: pca, umap

Idents(Sham_mb_RemDoub) <- "seurat_clusters"
DimPlot(Sham_mb_RemDoub, reduction = "umap", label = TRUE)

## MPTP
MPTP_mb <- CreateSeuratObject(counts = MPTP.data, project = "MPTP mouse midbrain", min.cells = 3)
MPTP_mb@meta.data$status <- c(rep("MPTP", ncol(MPTP.data)))
MPTP_mb
## An object of class Seurat 
## 21552 features across 8804 samples within 1 assay 
## Active assay: RNA (21552 features, 0 variable features)

MPTP_mb <- NormalizeData(MPTP_mb, normalization.method = "LogNormalize", scale.factor = 10000)
MPTP_mb <- FindVariableFeatures(MPTP_mb, selection.method = "vst", nfeatures = 3000)

all.genes <- rownames(MPTP_mb)
MPTP_mb <- ScaleData(MPTP_mb, features = all.genes, verbose = FALSE)
MPTP_mb <- RunPCA(MPTP_mb, npcs = 30, verbose = FALSE)
VizDimLoadings(MPTP_mb, dims = 1:20, reduction = "pca")
MPTP_mb <- JackStraw(MPTP_mb, num.replicate = 200)
MPTP_mb <- ScoreJackStraw(MPTP_mb, dims = 1:20)
JackStrawPlot(MPTP_mb, dims = 1:20)
ElbowPlot(MPTP_mb)

MPTP_mb <- FindNeighbors(MPTP_mb, reduction = "pca", dims = 1:30)
MPTP_mb <- FindClusters(MPTP_mb, resolution = 0.3)
MPTP_mb <- RunUMAP(MPTP_mb, dims = 1:30)
MPTP_mb_P1 <- DimPlot(MPTP_mb, reduction = "umap", label = TRUE)
MPTP_mb
## An object of class Seurat 
## 21552 features across 8804 samples within 1 assay 
## Active assay: RNA (21552 features, 3000 variable features)
## 2 dimensional reductions calculated: pca, umap

sweep.res.list_MPTP_mb <- paramSweep_v3(MPTP_mb, PCs = 1:30, sct = FALSE)
sweep.stats_MPTP_mb <- summarizeSweep(sweep.res.list_MPTP_mb, GT = FALSE)
bcmvn_MPTP_mb <- find.pK(sweep.stats_MPTP_mb)
opt_pK_MPTP_mb <- as.numeric(as.vector(bcmvn_MPTP_mb$pK[which.max(bcmvn_MPTP_mb$BCmetric)]))
print(opt_pK_MPTP_mb)
# [1] 0.005

annotations_MPTP_mb <- MPTP_mb@meta.data$seurat_clusters
homotypic.prop_MPTP_mb <- modelHomotypic(annotations_MPTP_mb)
nExp_poi_MPTP_mb <- round(0.069*nrow(MPTP_mb@meta.data))
nExp_poi.adj_MPTP_mb <- round(nExp_poi_MPTP_mb*(1-homotypic.prop_MPTP_mb))

MPTP_mb <- doubletFinder_v3(MPTP_mb, PCs = 1:30, pN = 0.25, pK = opt_pK_MPTP_mb, nExp = nExp_poi.adj_MPTP_mb, 
                          reuse.pANN = FALSE, sct = FALSE)  
MPTP_mb_P2 <- DimPlot(MPTP_mb, reduction = "umap", group.by = "DF.classifications_0.25_0.005_533")
CombinePlots(plots = list(MPTP_mb_P1, MPTP_mb_P2), ncol = 1)

table(MPTP_mb@meta.data$DF.classifications_0.25_0.005_533)
## Doublet Singlet 
## 533     8271 

Idents(MPTP_mb) <- "DF.classifications_0.25_0.005_533"
table(MPTP_mb@active.ident)

MPTP_mb_RemDoub <- subset(MPTP_mb, idents = "Singlet")
MPTP_mb_RemDoub
## An object of class Seurat 
## 21552 features across 8271 samples within 1 assay 
## Active assay: RNA (21552 features, 3000 variable features)
## 2 dimensional reductions calculated: pca, umap

Idents(MPTP_mb_RemDoub) <- "seurat_clusters"
DimPlot(MPTP_mb_RemDoub, reduction = "umap", label = TRUE)


# Sham clustering
Mice_Sham_mb <- CreateSeuratObject(counts = Sham_mb_RemDoub@assays$RNA@counts, project = "Sham mouse midbrain", min.cells = 3)
Mice_Sham_mb@meta.data$status <- Sham_mb_RemDoub@meta.data$status
Mice_Sham_mb
## An object of class Seurat 
## 21508 features across 9965 samples within 1 assay 
## Active assay: RNA (21508 features, 0 variable features)

Mice_Sham_mb[["percent.mt"]] <- PercentageFeatureSet(Mice_Sham_mb, pattern = "^mt-")

Mice_Sham_mb_QCP1 <- VlnPlot(Mice_Sham_mb, features = "nFeature_RNA", pt.size = 0.5) + scale_y_continuous(breaks = seq(0,10000, by = 2500), limits = c(0,10000)) + NoLegend()
Mice_Sham_mb_QCP2 <- VlnPlot(Mice_Sham_mb, features = "nCount_RNA", pt.size = 0.5) + scale_y_continuous(breaks = seq(0,60000, by = 12000), limits = c(0,60000)) + NoLegend()
Mice_Sham_mb_QCP3 <- VlnPlot(Mice_Sham_mb, features = "percent.mt", pt.size = 0.5) + scale_y_continuous(breaks = seq(0,40, by = 8), limits = c(0,40)) + NoLegend()

wrap_plots(plots = list(Mice_Sham_mb_QCP1, Mice_Sham_mb_QCP2, Mice_Sham_mb_QCP3), ncol = 3)

Mice_Sham_mb <- subset(Mice_Sham_mb, subset = nFeature_RNA > 200 & percent.mt < 3)
Mice_Sham_mb <- subset(Mice_Sham_mb, subset = nFeature_RNA < 6000)
Mice_Sham_mb
## An object of class Seurat 
## 21508 features across 9919 samples within 1 assay 
## Active assay: RNA (21508 features, 0 variable features)

Mice_Sham_mb <- NormalizeData(Mice_Sham_mb, normalization.method = "LogNormalize", scale.factor = 10000)
Mice_Sham_mb <- FindVariableFeatures(Mice_Sham_mb, selection.method = "vst", nfeatures = 3000)
Mice_Sham_mb
## An object of class Seurat 
## 21508 features across 9919 samples within 1 assay 
## Active assay: RNA (21508 features, 3000 variable features)

all.genes <- rownames(Mice_Sham_mb)
Mice_Sham_mb <- ScaleData(Mice_Sham_mb, features = all.genes,verbose = FALSE)
Mice_Sham_mb <- RunPCA(Mice_Sham_mb, npcs = 30, verbose = FALSE)
VizDimLoadings(Mice_Sham_mb, dims = 1:20, reduction = "pca")
Mice_Sham_mb <- JackStraw(Mice_Sham_mb, num.replicate = 100)
Mice_Sham_mb <- ScoreJackStraw(Mice_Sham_mb, dims = 1:20)
JackStrawPlot(Mice_Sham_mb, dims = 1:20)
ElbowPlot(Mice_Sham_mb)

Mice_Sham_mb <- FindNeighbors(Mice_Sham_mb, reduction = "pca", dims = 1:30)
Mice_Sham_mb <- FindClusters(Mice_Sham_mb, resolution = 0.3)
Mice_Sham_mb <- RunUMAP(Mice_Sham_mb, dims = 1:30)
DimPlot(Mice_Sham_mb, reduction = "umap", label = TRUE)

Mice_Sham_mb
## An object of class Seurat 
## 21508 features across 9919 samples within 1 assay 
## Active assay: RNA (21508 features, 3000 variable features)
## 2 dimensional reductions calculated: pca, umap

Mice_Sham_mb.allmarkers <- FindAllMarkers(Mice_Sham_mb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(Mice_Sham_mb.allmarkers, file = "Mice_Sham_mb.allmarkers.txt", col.names = TRUE, sep = "\t", quote = FALSE)

FeaturePlot(Mice_Sham_mb, features = c("Stmn2", "Thy1", "Eno2", "Rbfox3", "Snap25", "Syt1"), label = TRUE)
VlnPlot(Mice_Sham_mb, features = c("Stmn2", "Thy1", "Eno2", "Rbfox3", "Snap25", "Syt1"), stack = TRUE, same.y.lims = TRUE, flip = TRUE, pt.size = 0)

FeaturePlot(Mice_Sham_mb, features = c("Opalin", "Mog", "Plp1", "Mbp", "Cldn11"), label = TRUE)
VlnPlot(Mice_Sham_mb, features = c("Opalin", "Mog", "Plp1", "Mbp", "Cldn11"), stack = TRUE, same.y.lims = TRUE, flip = TRUE, pt.size = 0)

FeaturePlot(Mice_Sham_mb, features = c("Vcan", "Bcan", "Pdgfra", "Olig2"), label = TRUE)
VlnPlot(Mice_Sham_mb, features = c("Vcan", "Bcan", "Pdgfra", "Olig2"), stack = TRUE, same.y.lims = TRUE, flip = TRUE, pt.size = 0)

FeaturePlot(Mice_Sham_mb, features = c("Aqp4", "Gfap", "Slc1a3", "Slc1a2", "Slc4a4", "Ntsr2", "Gja1", "Aldh1l1"), label = TRUE)
VlnPlot(Mice_Sham_mb, features = c("Aqp4", "Gfap", "Slc1a3", "Slc1a2", "Slc4a4", "Ntsr2", "Gja1", "Aldh1l1"), stack = TRUE, same.y.lims = TRUE,
        flip = TRUE, pt.size = 0)

FeaturePlot(Mice_Sham_mb, features = c("Cx3cr1", "P2ry12", "Csf1r"), label = TRUE)
VlnPlot(Mice_Sham_mb, features = c("Cx3cr1", "P2ry12", "Csf1r"), stack = TRUE, same.y.lims = TRUE, flip = TRUE, pt.size = 0)

FeaturePlot(Mice_Sham_mb, features = c("Flt1", "Cldn5", "Ptprb", "Pecam1"), label = TRUE)
VlnPlot(Mice_Sham_mb, features = c("Flt1", "Cldn5", "Ptprb", "Pecam1"), stack = TRUE, same.y.lims = TRUE, flip = TRUE, pt.size = 0)

FeaturePlot(Mice_Sham_mb, features = c("Dcn", "Col1a1", "Kcnj8", "Abcc9"), label = TRUE)
VlnPlot(Mice_Sham_mb, features = c("Dcn", "Col1a1", "Kcnj8", "Abcc9"), stack = TRUE, same.y.lims = TRUE, flip = TRUE, pt.size = 0)

Mice_Sham_mb <- RenameIdents(Mice_Sham_mb, `0` = "Oligodendrocytes", `1` = "Neurons", `2` = "Astrocytes", `3` = "Oligodendrocytes", `4` = "Neurons",
                             `5` = "Neurons", `6` = "Endothelial cells", `7` = "Microglia", `8` = "OPCs", `9` = "Neurons", `10` = "Neurons",
                             `11` = "Pericytes", `12` = "Neurons", `13` = "Endothelial cells", `14` = "Fibroblast", `15` = "Astrocytes", `16` = "Neurons",
                             `17` = "Neurons", `18` = "Neurons")

Mice_Sham_mb$Major_celltype <- Mice_Sham_mb@active.ident

Majorcelltype_level <- c("Neurons", "Oligodendrocytes", "OPCs", "Astrocytes", "Microglia", "Endothelial cells", "Pericytes", "Fibroblast")
Mice_Sham_mb@meta.data$Major_celltype <- factor(x = Mice_Sham_mb@meta.data$Major_celltype, levels = Majorcelltype_level)
Idents(Mice_Sham_mb) <- "Major_celltype"

DimPlot(Mice_Sham_mb, reduction = "umap", label = TRUE)
VlnPlot(Mice_Sham_mb, features = c("Syt1", "Snap25", "Rbfox3", "Cldn11", "Mog", "Pdgfra", "Olig2", "Aqp4", "Ntsr2", "Cx3cr1", "P2ry12", "Flt1",
                                   "Cldn5", "Abcc9", "Kcnj8", "Dcn", "Col1a1"), stack = TRUE, same.y.lims = TRUE, flip = TRUE, pt.size = 0)

## extract Sham DA neurons data
Idents(Mice_Sham_mb) <- "seurat_clusters"

FeaturePlot(Mice_Sham_mb, features = c("Slc32a1", "Gad1", "Gad2", "Slc17a6", "Th", "Slc6a3", "Slc18a2", "Nr4a2", "Pitx3"), label = TRUE)
VlnPlot(Mice_Sham_mb, features = c("Slc32a1", "Gad1", "Gad2", "Slc17a6", "Th", "Slc6a3", "Slc18a2", "Nr4a2", "Pitx3"), stack = TRUE,
        same.y.lims = TRUE, flip = TRUE, pt.size = 1)

VlnPlot(Mice_Sham_mb, features = c("Th", "Slc6a3", "Slc18a2", "Nr4a2"), ncol = 2, pt.size = 1)

Sham_DANs_1 <- subset(Mice_Sham_mb, idents = "4")
Sham_DANs_1.data <- as.data.frame(Sham_DANs_1@assays$RNA@counts)
Sham_DANs_1.data$X <- rownames(Sham_DANs_1.data)

Sham_DANs_2 <- subset(Mice_Sham_mb, idents = c("1", "5", "9", "10", "12", "16", "17", "18"))
Sham_DANs_2 <- subset(Sham_DANs_2, subset = Th > 0, slot = 'counts')
Sham_DANs_2 <- subset(Sham_DANs_2, subset = Slc6a3 > 0, slot = 'counts')
Sham_DANs_2 <- subset(Sham_DANs_2, subset = Slc18a2 > 0, slot = 'counts')
Sham_DANs_2.data <- as.data.frame(Sham_DANs_2@assays$RNA@counts)
Sham_DANs_2.data$X <- rownames(Sham_DANs_2.data)

Sham_DANs.data <- merge(Sham_DANs_1.data, Sham_DANs_2.data, by = "X", all = TRUE)
rownames(Sham_DANs.data) <- Sham_DANs.data$X
Sham_DANs.data <- Sham_DANs.data[, -1]
Sham_DANs.data[is.na(Sham_DANs.data)] <- 0

Sham_DANs <- CreateSeuratObject(counts = Sham_DANs.data, project = "Sham mouse DANs", min.cells = 3)

Sham_DANs@meta.data$Major_celltype <- c(Sham_DANs_1@meta.data$Major_celltype, Sham_DANs_2@meta.data$Major_celltype)
Sham_DANs@meta.data$status <- c(Sham_DANs_1@meta.data$status, Sham_DANs_2@meta.data$status)

Sham_DANs
## An object of class Seurat 
## 16696 features across 652 samples within 1 assay 
## Active assay: RNA (16696 features, 0 variable features)

table(Sham_DANs$Major_celltype)
table(Sham_DANs$status)

Sham_DANs[["percent.mt"]] <- PercentageFeatureSet(Sham_DANs, pattern = "^mt-")

Sham_DANs_QCP1 <- VlnPlot(Sham_DANs, features = "nFeature_RNA", pt.size = 0.5) + scale_y_continuous(breaks = seq(0,6000, by = 1200), limits = c(0,6000)) + NoLegend()
Sham_DANs_QCP2 <- VlnPlot(Sham_DANs, features = "nCount_RNA", pt.size = 0.5) + scale_y_continuous(breaks = seq(0,30000, by = 6000), limits = c(0,30000)) + NoLegend()
Sham_DANs_QCP3 <- VlnPlot(Sham_DANs, features = "percent.mt", pt.size = 0.5) + scale_y_continuous(breaks = seq(0,3, by = 0.6), limits = c(0,3)) + NoLegend()

wrap_plots(plots = list(Sham_DANs_QCP1, Sham_DANs_QCP2, Sham_DANs_QCP3), ncol = 3)

Sham_DANs <- NormalizeData(Sham_DANs, normalization.method = "LogNormalize", scale.factor = 10000)
Sham_DANs <- FindVariableFeatures(Sham_DANs, selection.method = "vst", nfeatures = 3000)


# MPTP clustering
Mice_MPTP_mb <- CreateSeuratObject(counts = MPTP_mb_RemDoub@assays$RNA@counts, project = "MPTP mouse mb", min.cells = 3)
Mice_MPTP_mb@meta.data$status <- MPTP_mb_RemDoub@meta.data$status

Mice_MPTP_mb
## An object of class Seurat 
## 21270 features across 8271 samples within 1 assay 
## Active assay: RNA (21270 features, 0 variable features)

Mice_MPTP_mb[["percent.mt"]] <- PercentageFeatureSet(Mice_MPTP_mb, pattern = "^mt-")

Mice_MPTP_mb_QCP1 <- VlnPlot(Mice_MPTP_mb, features = "nFeature_RNA", pt.size = 0.5) + scale_y_continuous(breaks = seq(0,11000, by = 2200), limits = c(0,11000)) + NoLegend()
Mice_MPTP_mb_QCP2 <- VlnPlot(Mice_MPTP_mb, features = "nCount_RNA", pt.size = 0.5) + scale_y_continuous(breaks = seq(0,90000, by = 18000), limits = c(0,90000)) + NoLegend()
Mice_MPTP_mb_QCP3 <- VlnPlot(Mice_MPTP_mb, features = "percent.mt", pt.size = 0.5) + scale_y_continuous(breaks = seq(0,40, by = 8), limits = c(-1,40)) + NoLegend()

wrap_plots(plots = list(Mice_MPTP_mb_QCP1, Mice_MPTP_mb_QCP2, Mice_MPTP_mb_QCP3), ncol = 3)

Mice_MPTP_mb <- subset(Mice_MPTP_mb, subset = nFeature_RNA > 200 & percent.mt < 3)
Mice_MPTP_mb <- subset(Mice_MPTP_mb, subset = nFeature_RNA < 6000)
Mice_MPTP_mb
## An object of class Seurat 
## 21270 features across 8230 samples within 1 assay 
## Active assay: RNA (21270 features, 0 variable features)

Mice_MPTP_mb <- NormalizeData(Mice_MPTP_mb, normalization.method = "LogNormalize", scale.factor = 10000)
Mice_MPTP_mb <- FindVariableFeatures(Mice_MPTP_mb, selection.method = "vst", nfeatures = 3000)
Mice_MPTP_mb
## An object of class Seurat 
## 21270 features across 8230 samples within 1 assay 
## Active assay: RNA (21270 features, 3000 variable features)

all.genes <- rownames(Mice_MPTP_mb)
Mice_MPTP_mb <- ScaleData(Mice_MPTP_mb, features = all.genes,verbose = FALSE)
Mice_MPTP_mb <- RunPCA(Mice_MPTP_mb, npcs = 30, verbose = FALSE)
VizDimLoadings(Mice_MPTP_mb, dims = 1:20, reduction = "pca")
Mice_MPTP_mb <- JackStraw(Mice_MPTP_mb, num.replicate = 100)
Mice_MPTP_mb <- ScoreJackStraw(Mice_MPTP_mb, dims = 1:20)
JackStrawPlot(Mice_MPTP_mb, dims = 1:20)
ElbowPlot(Mice_MPTP_mb)

Mice_MPTP_mb <- FindNeighbors(Mice_MPTP_mb, reduction = "pca", dims = 1:30)
Mice_MPTP_mb <- FindClusters(Mice_MPTP_mb, resolution = 0.2)
Mice_MPTP_mb <- RunUMAP(Mice_MPTP_mb, dims = 1:30)
DimPlot(Mice_MPTP_mb, reduction = "umap", label = TRUE)

Mice_MPTP_mb
## An object of class Seurat 
## 21270 features across 8230 samples within 1 assay 
## Active assay: RNA (21270 features, 3000 variable features)
## 2 dimensional reductions calculated: pca, umap

Mice_MPTP_mb.allmarkers <- FindAllMarkers(Mice_MPTP_mb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(Mice_MPTP_mb.allmarkers, file = "Mice_MPTP_mb.allmarkers.txt", col.names = TRUE, sep = "\t", quote = FALSE)

FeaturePlot(Mice_MPTP_mb, features = c("Stmn2", "Thy1", "Eno2", "Rbfox3", "Snap25", "Syt1"), label = TRUE)
VlnPlot(Mice_MPTP_mb, features = c("Stmn2", "Thy1", "Eno2", "Rbfox3", "Snap25", "Syt1"), stack = TRUE, same.y.lims = TRUE, flip = TRUE, pt.size = 0)

FeaturePlot(Mice_MPTP_mb, features = c("Opalin", "Mog", "Plp1", "Mbp", "Cldn11"), label = TRUE)
VlnPlot(Mice_MPTP_mb, features = c("Opalin", "Mog", "Plp1", "Mbp", "Cldn11"), stack = TRUE, same.y.lims = TRUE, flip = TRUE, pt.size = 0)

FeaturePlot(Mice_MPTP_mb, features = c("Vcan", "Bcan", "Pdgfra", "Olig2"), label = TRUE)
VlnPlot(Mice_MPTP_mb, features = c("Vcan", "Bcan", "Pdgfra", "Olig2"), stack = TRUE, same.y.lims = TRUE, flip = TRUE, pt.size = 0)

FeaturePlot(Mice_MPTP_mb, features = c("Aqp4", "Gfap", "Slc1a3", "Slc1a2", "Slc4a4", "Ntsr2", "Gja1", "Aldh1l1"), label = TRUE)
VlnPlot(Mice_MPTP_mb, features = c("Aqp4", "Gfap", "Slc1a3", "Slc1a2", "Slc4a4", "Ntsr2", "Gja1", "Aldh1l1"), stack = TRUE, same.y.lims = TRUE,
        flip = TRUE, pt.size = 0)

FeaturePlot(Mice_MPTP_mb, features = c("Cx3cr1", "P2ry12", "Csf1r"), label = TRUE)
VlnPlot(Mice_MPTP_mb, features = c("Cx3cr1", "P2ry12", "Csf1r"), stack = TRUE, same.y.lims = TRUE, flip = TRUE, pt.size = 0)

FeaturePlot(Mice_MPTP_mb, features = c("Flt1", "Cldn5", "Ptprb", "Pecam1"), label = TRUE)
VlnPlot(Mice_MPTP_mb, features = c("Flt1", "Cldn5", "Ptprb", "Pecam1"), stack = TRUE, same.y.lims = TRUE, flip = TRUE, pt.size = 0)

FeaturePlot(Mice_MPTP_mb, features = c("Dcn", "Col1a1", "Kcnj8", "Abcc9"), label = TRUE)
VlnPlot(Mice_MPTP_mb, features = c("Dcn", "Col1a1", "Kcnj8", "Abcc9"), stack = TRUE, same.y.lims = TRUE, flip = TRUE, pt.size = 0)

Mice_MPTP_mb <- RenameIdents(Mice_MPTP_mb, `0` = "Oligodendrocytes", `1` = "Neurons", `2` = "Astrocytes", `3` = "Neurons", `4` = "Endothelial cells",
                           `5` = "Neurons", `6` = "Neurons", `7` = "Microglia", `8` = "OPCs", `9` = "Neurons", `10` = "Fibroblast",
                           `11` = "Pericytes", `12` = "Neurons", `13` = "Neurons", `14` = "Neurons", `15` = "Neurons", `16` = "Astrocytes")

Mice_MPTP_mb$Major_celltype <- Mice_MPTP_mb@active.ident
Mice_MPTP_mb@meta.data$Major_celltype <- factor(x = Mice_MPTP_mb@meta.data$Major_celltype, levels = Majorcelltype_level)
Idents(Mice_MPTP_mb) <- "Major_celltype"

DimPlot(Mice_MPTP_mb, reduction = "umap", label = TRUE)
VlnPlot(Mice_MPTP_mb, features = c("Syt1", "Snap25", "Rbfox3", "Cldn11", "Mog", "Pdgfra", "Olig2", "Aqp4", "Ntsr2", "Cx3cr1", "P2ry12", "Flt1",
                                 "Cldn5", "Abcc9", "Kcnj8", "Dcn", "Col1a1"), stack = TRUE, same.y.lims = TRUE, flip = TRUE, pt.size = 0)

## extract MPTP DA neurons data
Idents(Mice_MPTP_mb) <- "seurat_clusters"

FeaturePlot(Mice_MPTP_mb, features = c("Slc32a1", "Gad1", "Gad2", "Slc17a6", "Th", "Slc6a3", "Slc18a2", "Nr4a2", "Pitx3"), label = TRUE)
VlnPlot(Mice_MPTP_mb, features = c("Slc32a1", "Gad1", "Gad2", "Slc17a6", "Th", "Slc6a3", "Slc18a2", "Nr4a2", "Pitx3"), stack = TRUE,
        same.y.lims = TRUE, flip = TRUE, pt.size = 1)

VlnPlot(Mice_MPTP_mb, features = c("Th", "Slc6a3", "Slc18a2", "Nr4a2"), ncol = 2, pt.size = 1)

MPTP_DANs_1 <- subset(Mice_MPTP_mb, idents = "5")
MPTP_DANs_1.data <- as.data.frame(MPTP_DANs_1@assays$RNA@counts)
MPTP_DANs_1.data$X <- rownames(MPTP_DANs_1.data)

MPTP_DANs_2 <- subset(Mice_MPTP_mb, idents = c("1", "3", "6", "9", "12", "13", "14", "15"))
MPTP_DANs_2 <- subset(MPTP_DANs_2, subset = Th > 0, slot = 'counts')
MPTP_DANs_2 <- subset(MPTP_DANs_2, subset = Slc6a3 > 0, slot = 'counts')
MPTP_DANs_2 <- subset(MPTP_DANs_2, subset = Slc18a2 > 0, slot = 'counts')
MPTP_DANs_2.data <- as.data.frame(MPTP_DANs_2@assays$RNA@counts)
MPTP_DANs_2.data$X <- rownames(MPTP_DANs_2.data)

MPTP_DANs.data <- merge(MPTP_DANs_1.data, MPTP_DANs_2.data, by = "X", all = TRUE)
rownames(MPTP_DANs.data) <- MPTP_DANs.data$X
MPTP_DANs.data <- MPTP_DANs.data[, -1]
MPTP_DANs.data[is.na(MPTP_DANs.data)] <- 0

MPTP_DANs <- CreateSeuratObject(counts = MPTP_DANs.data, project = "MPTP mouse DANs", min.cells = 3)

MPTP_DANs@meta.data$Major_celltype <- c(MPTP_DANs_1@meta.data$Major_celltype, MPTP_DANs_2@meta.data$Major_celltype)
MPTP_DANs@meta.data$status <- c(MPTP_DANs_1@meta.data$status, MPTP_DANs_2@meta.data$status)

MPTP_DANs
## An object of class Seurat 
## 15410 features across 371 samples within 1 assay 
## Active assay: RNA (15410 features, 0 variable features)

table(MPTP_DANs$Major_celltype)
table(MPTP_DANs$status)

MPTP_DANs[["percent.mt"]] <- PercentageFeatureSet(MPTP_DANs, pattern = "^mt-")

MPTP_DANs_QCP1 <- VlnPlot(MPTP_DANs, features = "nFeature_RNA", pt.size = 0.5) + scale_y_continuous(breaks = seq(0,6000, by = 1200), limits = c(0,6000)) + NoLegend()
MPTP_DANs_QCP2 <- VlnPlot(MPTP_DANs, features = "nCount_RNA", pt.size = 0.5) + scale_y_continuous(breaks = seq(0,30000, by = 6000), limits = c(0,30000)) + NoLegend()
MPTP_DANs_QCP3 <- VlnPlot(MPTP_DANs, features = "percent.mt", pt.size = 0.5) + scale_y_continuous(breaks = seq(0,3, by = 0.6), limits = c(0,3)) + NoLegend()

wrap_plots(plots = list(MPTP_DANs_QCP1, MPTP_DANs_QCP2, MPTP_DANs_QCP3), ncol = 3)

MPTP_DANs <- NormalizeData(MPTP_DANs, normalization.method = "LogNormalize", scale.factor = 10000)
MPTP_DANs <- FindVariableFeatures(MPTP_DANs, selection.method = "vst", nfeatures = 3000)


# Perform Sham and MPTP dataset integration
Mice_mb.list <- c(Mice_Sham_mb, Mice_MPTP_mb)
features <- SelectIntegrationFeatures(object.list = Mice_mb.list)

Mice_mb.anchors <- FindIntegrationAnchors(object.list = Mice_mb.list, anchor.features = features)
Mice_mb.combined <- IntegrateData(anchorset = Mice_mb.anchors)

DefaultAssay(Mice_mb.combined) <- "integrated"

all.genes <- rownames(Mice_mb.combined)
Mice_mb.combined <- ScaleData(Mice_mb.combined, features = all.genes,verbose = FALSE)
Mice_mb.combined <- RunPCA(Mice_mb.combined, npcs = 30, verbose = FALSE)

VizDimLoadings(Mice_mb.combined, dims = 1:12, reduction = "pca")
Mice_mb.combined <- JackStraw(Mice_mb.combined, num.replicate = 100)
Mice_mb.combined <- ScoreJackStraw(Mice_mb.combined, dims = 1:20)
JackStrawPlot(Mice_mb.combined, dims = 1:20)
ElbowPlot(Mice_mb.combined)

Mice_mb.combined <- FindNeighbors(Mice_mb.combined, reduction = "pca", dims = 1:20)
Mice_mb.combined <- FindClusters(Mice_mb.combined, resolution = 0.2)
Mice_mb.combined <- RunUMAP(Mice_mb.combined, reduction = "pca", dims = 1:20)

DimPlot(Mice_mb.combined, reduction = "umap", label = TRUE, repel = TRUE)

Mice_mb.combined@meta.data$status <- factor(x = Mice_mb.combined@meta.data$status, levels = c("Sham", "MPTP"))
DimPlot(Mice_mb.combined, reduction = "umap", group.by = "status")

Mice_mb.combined
## An object of class Seurat 
## 24060 features across 18149 samples within 2 assays 
## Active assay: integrated (2000 features, 2000 variable features)
## 1 other assay present: RNA
## 2 dimensional reductions calculated: pca, umap

Mice_mb.combined.allmarkers <- FindAllMarkers(Mice_mb.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(Mice_mb.combined.allmarkers, file = "Mice_mb.combined.allmarkers.txt", col.names = TRUE, sep = "\t", quote = FALSE)

DefaultAssay(Mice_mb.combined) <- "RNA"
all.genes <- rownames(Mice_mb.combined)
Mice_mb.combined <- ScaleData(Mice_mb.combined, features = all.genes,verbose = FALSE)

FeaturePlot(Mice_mb.combined, features = c("Stmn2", "Thy1", "Eno2", "Rbfox3", "Snap25", "Syt1", "Slc6a3", "Th"), label = TRUE)
VlnPlot(Mice_mb.combined, features = c("Stmn2", "Thy1", "Eno2", "Rbfox3", "Snap25", "Syt1", "Slc6a3", "Th"), stack = TRUE, same.y.lims = TRUE, flip = TRUE, pt.size = 0)

FeaturePlot(Mice_mb.combined, features = c("Slc17a6", "Slc32a1", "Gad1", "Gad2", "Slc6a3", "Th", "Slc18a2", "Nr4a2"), label = TRUE)
VlnPlot(Mice_mb.combined, features = c("Slc17a6", "Slc32a1", "Gad1", "Gad2", "Slc6a3", "Th", "Slc18a2", "Nr4a2"), stack = TRUE, same.y.lims = TRUE, flip = TRUE, pt.size = 0)

FeaturePlot(Mice_mb.combined, features = c("Opalin", "Mog", "Plp1", "Mbp", "Cldn11"), label = TRUE)
VlnPlot(Mice_mb.combined, features = c("Opalin", "Mog", "Plp1", "Mbp", "Cldn11"), stack = TRUE, same.y.lims = TRUE, flip = TRUE, pt.size = 0)

FeaturePlot(Mice_mb.combined, features = c("Vcan", "Bcan", "Pdgfra", "Olig2"), label = TRUE)
VlnPlot(Mice_mb.combined, features = c("Vcan", "Bcan", "Pdgfra", "Olig2"), stack = TRUE, same.y.lims = TRUE, flip = TRUE, pt.size = 0)

FeaturePlot(Mice_mb.combined, features = c("Aqp4", "Gfap", "Slc1a3", "Slc1a2", "Slc4a4", "Ntsr2", "Gja1", "Aldh1l1"), label = TRUE)
VlnPlot(Mice_mb.combined, features = c("Aqp4", "Gfap", "Slc1a3", "Slc1a2", "Slc4a4", "Ntsr2", "Gja1", "Aldh1l1"), stack = TRUE, same.y.lims = TRUE,
        flip = TRUE, pt.size = 0)

FeaturePlot(Mice_mb.combined, features = c("Cx3cr1", "P2ry12", "Csf1r"), label = TRUE)
VlnPlot(Mice_mb.combined, features = c("Cx3cr1", "P2ry12", "Csf1r"), stack = TRUE, same.y.lims = TRUE, flip = TRUE, pt.size = 0)

FeaturePlot(Mice_mb.combined, features = c("Flt1", "Cldn5", "Ptprb", "Pecam1"), label = TRUE)
VlnPlot(Mice_mb.combined, features = c("Flt1", "Cldn5", "Ptprb", "Pecam1"), stack = TRUE, same.y.lims = TRUE, flip = TRUE, pt.size = 0)

FeaturePlot(Mice_mb.combined, features = c("Dcn", "Col1a1", "Kcnj8", "Abcc9"), label = TRUE)
VlnPlot(Mice_mb.combined, features = c("Dcn", "Col1a1", "Kcnj8", "Abcc9"), stack = TRUE, same.y.lims = TRUE, flip = TRUE, pt.size = 0)

Mice_mb.combined <- RenameIdents(Mice_mb.combined, `0` = "Oligodendrocytes", `1` = "Glu_GABA", `2` = "Astrocytes", `3` = "DAN",
                                 `4` = "Glu_GABA", `5` = "Endothelial cells", `6` = "Microglia", `7` = "Glu_GABA", `8` = "OPCs",
                                 `9` = "GABA", `10` = "GABA", `11` = "Pericytes", `12` = "Fibroblast", `13` = "Glu", `14` = "Astrocytes",
                                 `15` = "Endothelial cells", `16` = "GABA", `17` = "Glu")

Mice_mb.combined$Major_celltype <- Mice_mb.combined@active.ident

Majorcelltype_level <- c("Glu_GABA", "Glu", "GABA", "DAN", "Oligodendrocytes", "OPCs", "Astrocytes", "Microglia", "Endothelial cells",
                         "Pericytes", "Fibroblast")

Mice_mb.combined@meta.data$Major_celltype <- factor(x = Mice_mb.combined@meta.data$Major_celltype, levels = Majorcelltype_level)
Idents(Mice_mb.combined) <- "Major_celltype"

DimPlot(Mice_mb.combined, reduction = "umap", label = TRUE)
VlnPlot(Mice_mb.combined, features = c("Syt1", "Snap25", "Rbfox3", "Slc17a6", "Slc32a1", "Gad2", "Th", "Slc6a3", "Slc18a2", "Nr4a2",
                                       "Cldn11", "Mog", "Pdgfra", "Olig2", "Aqp4", "Ntsr2", "Cx3cr1", "P2ry12", "Flt1", "Cldn5", "Abcc9",
                                       "Kcnj8", "Dcn", "Col1a1"), stack = TRUE, same.y.lims = TRUE, flip = TRUE, pt.size = 0)

table(Mice_mb.combined$status)
## Sham MPTP 
## 9919 8230

table(Idents(Mice_mb.combined), Mice_mb.combined$status)
##                   Sham MPTP
## Glu_GABA          2063 2334
## Glu                164  196
## GABA               677  380
## DAN                632  387
## Oligodendrocytes  3169 2492
## OPCs               353  255
## Astrocytes        1512 1152
## Microglia          389  261
## Endothelial cells  641  469
## Pericytes          182  151
## Fibroblast         137  153

head(Mice_mb.combined, 2)

## Prepare data for integrated UMAP plotting
circ_data <- prepare_circlize_data(Mice_mb.combined, scale = 0.8)
set.seed(123456)

Major_celltype_colors <- rand_color(length(levels(Mice_mb.combined)))
status_colors <- rand_color(length(names(table(Mice_mb.combined$status))))
cluster_colors <- rand_color(length(names(table(Mice_mb.combined$seurat_clusters))))

plot_circlize(circ_data, do.label = T, pt.size = 0.2, col.use = c("#F8766D", "#DB8E00", "#AEA200", "#64B200", "#00BD5C", "#00C1A7",
                                                                  "#00BADE", "#00A6FF", "#B385FF", "#EF67EB", "#FF63B6"),
              bg.color = 'white', kde2d.n = 1000, repel = T, label.cex = 1)

add_track(circ_data, group = "status", colors = status_colors, track_num = 2)
add_track(circ_data, group = "seurat_clusters", colors = cluster_colors, track_num = 3)


# Perform DA neurons integration
Mice_DANs.list <- c(Sham_DANs, MPTP_DANs)
features <- SelectIntegrationFeatures(object.list = Mice_DANs.list)

Mice_DANs.anchors <- FindIntegrationAnchors(object.list = Mice_DANs.list, anchor.features = features)
Mice_DANs.combined <- IntegrateData(anchorset = Mice_DANs.anchors)

DefaultAssay(Mice_DANs.combined) <- "integrated"

all.genes <- rownames(Mice_DANs.combined)
Mice_DANs.combined <- ScaleData(Mice_DANs.combined, features = all.genes,verbose = FALSE)
Mice_DANs.combined <- RunPCA(Mice_DANs.combined, npcs = 30, verbose = FALSE)

VizDimLoadings(Mice_DANs.combined, dims = 1:12, reduction = "pca")
Mice_DANs.combined <- JackStraw(Mice_DANs.combined, num.replicate = 100)
Mice_DANs.combined <- ScoreJackStraw(Mice_DANs.combined, dims = 1:20)
JackStrawPlot(Mice_DANs.combined, dims = 1:20)
ElbowPlot(Mice_DANs.combined)

Mice_DANs.combined <- FindNeighbors(Mice_DANs.combined, reduction = "pca", dims = 1:20)
Mice_DANs.combined <- FindClusters(Mice_DANs.combined, resolution = 0.2)
Mice_DANs.combined <- RunUMAP(Mice_DANs.combined, reduction = "pca", dims = 1:20)

Mice_DANs.combined@meta.data$status <- factor(x = Mice_DANs.combined@meta.data$status, levels = c("Sham", "MPTP"))

P1 <- DimPlot(Mice_DANs.combined, reduction = "umap", label = TRUE, repel = TRUE)
P2 <- DimPlot(Mice_DANs.combined, reduction = "umap", group.by = "status")
wrap_plots(plots = list(P1, P2), ncol = 2)

DimPlot(Mice_DANs.combined, reduction = "umap", split.by = "status", label = TRUE)

Mice_DANs.combined <- RenameIdents(Mice_DANs.combined, `0` = "mDAN_1", `1` = "mDAN_2", `2` = "mDAN_3", `3` = "mDAN_4", `4` = "mDAN_5")
Mice_DANs.combined$DANs_subtype <- Mice_DANs.combined@active.ident
Idents(Mice_DANs.combined) <- "DANs_subtype"
DimPlot(Mice_DANs.combined, reduction = "umap", label = TRUE)

FeaturePlot(Mice_DANs.combined, features = c("Sox6", "Calb1"), label = TRUE)

table(Mice_DANs.combined$status)
## Sham MPTP
## 652  371

table(Idents(Mice_DANs.combined), Mice_DANs.combined$status)
##        Sham MPTP     
## mDAN_1  189 112
## mDAN_2  176 124
## mDAN_3  161  88
## mDAN_4   98  33 
## mDAN_5   28  14 

prop.table(table(Idents(Mice_DANs.combined), Mice_DANs.combined$status), margin = 2)
##              Sham       MPTP 
## mDAN_1 0.28987730 0.30188679
## mDAN_2 0.26993865 0.33423181 
## mDAN_3 0.24693252 0.23719677
## mDAN_4 0.15030675 0.08894879
## mDAN_5 0.04294479 0.03773585

Mice_mb.combined_meta.data <- Mice_mb.combined@meta.data
Mice_mb.combined_meta.data <- Mice_mb.combined_meta.data[, c(4,8)]
Mice_mb.combined_meta.data_sham <- subset(Mice_mb.combined_meta.data, status == "Sham")
Mice_mb.combined_meta.data_MPTP <- subset(Mice_mb.combined_meta.data, status == "MPTP")

rownames(Mice_mb.combined_meta.data_sham) <- gsub("_1","", rownames(Mice_mb.combined_meta.data_sham))
rownames(Mice_mb.combined_meta.data_MPTP) <- gsub("_2","", rownames(Mice_mb.combined_meta.data_MPTP))

Mice_DANs.combined_meta.data <- Mice_DANs.combined@meta.data
Mice_DANs.combined_meta.data <- Mice_DANs.combined_meta.data[, c(5,9)]
Mice_DANs.combined_meta.data_sham <- subset(Mice_DANs.combined_meta.data, status == "Sham")
Mice_DANs.combined_meta.data_MPTP <- subset(Mice_DANs.combined_meta.data, status == "MPTP")

Mice_mb.combined_meta.data_sham <- Mice_mb.combined_meta.data_sham[!rownames(Mice_mb.combined_meta.data_sham) %in% rownames(Mice_DANs.combined_meta.data_sham),]
Mice_mb.combined_meta.data_MPTP <- Mice_mb.combined_meta.data_MPTP[!rownames(Mice_mb.combined_meta.data_MPTP) %in% rownames(Mice_DANs.combined_meta.data_MPTP),]

table(Mice_mb.combined_meta.data_sham$Major_celltype)
## Glu_GABA     Glu        GABA         DAN         Oligodendrocytes    OPCs        Astrocytes 
## 2031         164        672          17          3169                353         1512 
## Microglia    Endothelial cells       Pericytes   Fibroblast 
## 389          641                     182         137

Mice_DANs.combined_meta.data_sham_DAN <- subset(Mice_mb.combined_meta.data_sham, Major_celltype == "DAN")
Mice_DANs.combined_meta.data_sham_DAN <- rownames(Mice_DANs.combined_meta.data_sham_DAN)

Sham_unknNeu_mat <- as.matrix(Mice_Sham_mb@assays$RNA@counts)
Sham_unknNeu_mat <- Sham_unknNeu_mat[,Mice_DANs.combined_meta.data_sham_DAN]
Sham_unknNeu_mat <- Sham_unknNeu_mat[c("Slc17a6", "Slc32a1", "Gad1", "Gad2", "Th", "Slc6a3", "Slc18a2"),]

table(Mice_mb.combined_meta.data_MPTP$Major_celltype)
## Glu_GABA     Glu        GABA         DAN         Oligodendrocytes    OPCs        Astrocytes 
## 2321         196        379          30          2492                255         1152 
## Microglia    Endothelial cells       Pericytes   Fibroblast 
## 261          469                     151         153 

Mice_DANs.combined_meta.data_MPTP_DAN <- subset(Mice_mb.combined_meta.data_MPTP, Major_celltype == "DAN")
Mice_DANs.combined_meta.data_MPTP_DAN <- rownames(Mice_DANs.combined_meta.data_MPTP_DAN)

MPTP_unknNeu_mat <- as.matrix(Mice_MPTP_mb@assays$RNA@counts)
MPTP_unknNeu_mat <- MPTP_unknNeu_mat[,Mice_DANs.combined_meta.data_MPTP_DAN]
MPTP_unknNeu_mat <- MPTP_unknNeu_mat[c("Slc17a6", "Slc32a1", "Gad1", "Gad2", "Th", "Slc6a3", "Slc18a2"),]

gplots::balloonplot(table(Mice_DANs.combined$status, Mice_DANs.combined$DANs_subtype)) 
Idents(Mice_DANs.combined) <- Mice_DANs.combined$status
table(Idents(Mice_DANs.combined))

DANs_degs = lapply(unique(Mice_DANs.combined$DANs_subtype), function(x){
  FindMarkers(Mice_DANs.combined[,Mice_DANs.combined$DANs_subtype == x],ident.1 = 'MPTP', ident.2 = 'Sham')
})

names(DANs_degs) <- unique(Mice_DANs.combined$DANs_subtype)

df2 <- data.frame()
top_n_df <- data.frame()

for (i in 1:length(DANs_degs)){
  dat <- DANs_degs[[i]]
  dat$DANs_subtype <- names(DANs_degs)[[i]]
  dat$gene <- rownames(dat)
  dat$change <- ifelse(dat$avg_log2FC > 0 & dat$p_val_adj < 0.05,'up',
                       ifelse(dat$avg_log2FC < 0 & dat$p_val_adj < 0.05,'down','nochange'))
  dat_sig <- dat[dat$change!= 'nochange',]
  up_top5 <- top_n(dat_sig[dat_sig$change == 'up',], 5, avg_log2FC)[,'gene']
  down_top5 <- top_n(dat_sig[dat_sig$change == 'down',], -5, avg_log2FC)[,'gene']
  label_gene <- c(up_top5, down_top5)
  top_n_df <- rbind(top_n_df, dat[label_gene,]) 
  df2 <- rbind(df2, dat)
}

rownames(df2) <- 1:nrow(df2)

df3 <- df2
df2 <- df2[df2$change!= 'nochange',]

p1 <- ggplot() + geom_jitter(data = df2, aes(x = DANs_subtype, y = avg_log2FC, color = change), width = 0.2, size = 2)
p1

p2 <- p1 + geom_jitter(data = top_n_df, aes(x = DANs_subtype, y = avg_log2FC, color = change), width = 0.2, size = 0.5)+
  geom_text_repel(data = top_n_df, aes(x = DANs_subtype, y = avg_log2FC, label = gene), size = 3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2

up_bar <- top_n(group_by(df2, DANs_subtype), 1, avg_log2FC)
down_bar <- top_n(group_by(df2, DANs_subtype), -1, avg_log2FC)

bar_df <- data.frame(row.names = up_bar$DANs_subtype, x = up_bar$DANs_subtype, up = up_bar$avg_log2FC + 0.3,
                     down = down_bar$avg_log2FC - 0.3)
bar_df <- bar_df[sort(unique(df2$DANs_subtype)),]

bar_df$down <- ifelse(bar_df$down >= 0, NA, bar_df$down)

ggplot() + geom_col(data = bar_df, aes(x = x, y = up), alpha = 0.2) + geom_col(data = bar_df, aes(x = x, y = down), alpha = 0.2)

p3 <- p2 + geom_col(data = bar_df, aes(x = x, y = up), alpha = 0.2) + geom_col(data = bar_df, aes(x = x, y = down), alpha = 0.2)
p3

box_df <- data.frame(row.names = up_bar$DANs_subtype, x = up_bar$DANs_subtype, y = 0)

p4 <- p3 + geom_tile(data = box_df, aes(x = x, y = y), height = 0.1, fill = c("#F8766D", "#DB8E00", "#AEA200"))
p4

write.table(DANs_degs$mDAN_1, file = "mDAN_1.MPTP_DEG.txt", col.names = TRUE, sep = "\t", quote = FALSE)
write.table(DANs_degs$mDAN_2, file = "mDAN_2.MPTP_DEG.txt", col.names = TRUE, sep = "\t", quote = FALSE)
write.table(DANs_degs$mDAN_3, file = "mDAN_3.MPTP_DEG.txt", col.names = TRUE, sep = "\t", quote = FALSE)
write.table(DANs_degs$mDAN_4, file = "mDAN_4.MPTP_DEG.txt", col.names = TRUE, sep = "\t", quote = FALSE)
write.table(DANs_degs$mDAN_5, file = "mDAN_5.MPTP_DEG.txt", col.names = TRUE, sep = "\t", quote = FALSE)

mDAN_4.MPTP_DEG_GO <- subset(DANs_degs$mDAN_4, subset = p_val_adj < 0.05)
mDAN_4.MPTP_DEG_GO <- subset(mDAN_4.MPTP_DEG_GO, subset = avg_log2FC < 0)
mDAN_4.MPTP_DEG_GO <- as.character(unique(rownames(mDAN_4.MPTP_DEG_GO)))
mDAN_4.MPTP_DEG_GO.df <- bitr(mDAN_4.MPTP_DEG_GO, fromType = "SYMBOL", toType = c("ENTREZID", "ENSEMBL"), OrgDb = org.Mm.eg.db)
mDAN_4.MPTP_DEG_GO <- unique(mDAN_4.MPTP_DEG_GO.df$ENTREZID)

mDAN_4.MPTP_DEG_ego <- enrichGO(gene = mDAN_4.MPTP_DEG_GO, OrgDb = org.Mm.eg.db, ont = "MF", readable = TRUE)
mDAN_4.MPTP_DEG_ego <- enrichGO(gene = mDAN_4.MPTP_DEG_GO, OrgDb = org.Mm.eg.db, ont = "BP", readable = TRUE)
mDAN_4.MPTP_DEG_ego <- enrichGO(gene = mDAN_4.MPTP_DEG_GO, OrgDb = org.Mm.eg.db, ont = "CC", readable = TRUE)

mDAN_4.MPTP_DEG_ego  <- clusterProfiler::simplify(mDAN_4.MPTP_DEG_ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
dotplot(mDAN_4.MPTP_DEG_ego, showCategory = 10, title = "mDAN_4.MPTP_DEG_ego_CC")
goplot(mDAN_4.MPTP_DEG_ego)

mDAN_4.MPTP_DEG_kegg <- enrichKEGG(gene = mDAN_4.MPTP_DEG_GO, keyType = "kegg", organism = 'mmu', pAdjustMethod = "fdr",
                                pvalueCutoff = 0.5, qvalueCutoff = 0.5)
dotplot(mDAN_4.MPTP_DEG_kegg, title = "mDAN_4.MPTP_DEG_kegg")

## subset sham DANs
Mice_DANs.sham <- subset(Mice_DANs.combined, idents = "Sham")
Mice_DANs.MPTP <- subset(Mice_DANs.combined, idents = "MPTP")
Idents(Mice_DANs.sham) <- "DANs_subtype"
Idents(Mice_DANs.MPTP) <- "DANs_subtype"

DefaultAssay(Mice_DANs.sham) <- "RNA"
all.genes <- rownames(Mice_DANs.sham)
Mice_DANs.sham <- ScaleData(Mice_DANs.sham, features = all.genes, verbose = FALSE)

FeaturePlot(Mice_DANs.sham, features = c("Sox6", "Calb1"), label = TRUE)

DANs_subtype.markers <- FindAllMarkers(Mice_DANs.sham, only.pos = TRUE)
write.table(DANs_subtype.markers, file = "DANs_subtype.markers.txt", col.names = TRUE, sep = "\t", quote = FALSE)

DANs_subtype.markers_cosg <- cosg(Mice_DANs.sham, groups = "all", assay = "RNA", slot = "data", mu = 1, remove_lowly_expressed = T,
                                  expressed_pct = 0.5, n_genes_user = 100)
write.table(DANs_subtype.markers_cosg$names, file = "DANs_subtype.markers_cosg.txt", col.names = TRUE, sep = "\t", quote = FALSE)

VlnPlot(Mice_DANs.sham, features = c("Sox6", "Calb1", "Sorcs2", "Lmo3", "Rph3a", "Slc44a5", "Kcnab1", "Sorcs3", "Tmsb4x", "Tuba1a", "Fstl4", "Megf11"), 
        stack = TRUE, same.y.lims = TRUE, flip = TRUE, pt.size = 0)


# hdWGCNA
theme_set(theme_cowplot())

Mice_DANs.sham <- FindVariableFeatures(Mice_DANs.sham, selection.method = "vst", nfeatures = 2000)

MB_DANs <- SetupForWGCNA(Mice_DANs.sham, gene_select = "variable", wgcna_name = "MB_DANs")
length(MB_DANs@misc$MB_DANs$wgcna_genes)

MB_DANs <- MetacellsByGroups(seurat_obj = MB_DANs, group.by = "DANs_subtype", reduction = "pca", k = 15, min_cells = 20, max_shared = 5, ident.group = "DANs_subtype")
MB_DANs <- NormalizeMetacells(MB_DANs)

MB_DANs <- SetDatExpr(MB_DANs, group_name = "mDAN_4", group.by = "DANs_subtype", assay = "RNA", slot = "data")
MB_DANs <- TestSoftPowers(MB_DANs, networkType = "signed")
plot_list <- PlotSoftPowers(MB_DANs)
wrap_plots(plot_list, ncol = 2)

power_table <- GetPowerTable(MB_DANs)
head(power_table)

## construct co-expression network
MB_DANs <- ConstructNetwork(MB_DANs, soft_power = 10, setDatExpr = FALSE, tom_name = "mDAN_4")
PlotDendrogram(MB_DANs, main = "mDAN_4 hdWGCNA Dendrogram")

MB_DANs@misc$MB_DANs$wgcna_modules %>% head
write.table(MB_DANs@misc$MB_DANs$wgcna_modules, file = "mDAN_4_modules.txt", col.names = TRUE, sep = "\t", quote = FALSE)

table(MB_DANs@misc$MB_DANs$wgcna_modules$module)
## yellow  brown   turquoise   red     black    blue    pink    green    grey    magenta    purple  greenyellow 
## 154     191     237         115     99       192     87      117      81      85         66      51 

Tom <- GetTOM(MB_DANs)

## Compute harmonized module eigengenes
MB_DANs <- ScaleData(MB_DANs, features =  VariableFeatures(MB_DANs))
MB_DANs <- ModuleEigengenes(MB_DANs, npcs = 20, reduction.use = "pca")

## harmonized module eigengenes
hMEs <- GetMEs(MB_DANs)
head(hMEs)
MEs <- GetMEs(MB_DANs, harmonized = FALSE)
head(MEs)

## Compute module connectivity
MB_DANs <- ModuleConnectivity(MB_DANs, group.by = "DANs_subtype", group_name = "mDAN_4")
MB_DANs <- ResetModuleNames(MB_DANs, new_name = "mDAN_4-M")

PlotKMEs(MB_DANs, ncol = 5)

modules <- GetModules(MB_DANs)

hub_df <- GetHubGenes(MB_DANs, n_hubs = 10)
head(hub_df)
write.table(hub_df, file = "hub_genes.txt", col.names = TRUE, sep = "\t", quote = FALSE)

MB_DANs <- ModuleExprScore(MB_DANs, n_genes = 25, method = "Seurat")

## Visualization
plot_list <- ModuleFeaturePlot(MB_DANs, features = "hMEs", order = TRUE)
wrap_plots(plot_list, ncol = 4) 

plot_list <- ModuleFeaturePlot(MB_DANs, features = "scores", order = "shuffle", ucell = TRUE)
wrap_plots(plot_list, ncol = 4)

MEs <- GetMEs(MB_DANs, harmonized = TRUE)
mods <- colnames(MEs); mods <- mods[mods != "grey"]

MB_DANs@meta.data <- cbind(MB_DANs@meta.data, MEs)
P <- DotPlot(MB_DANs, features = mods, group.by = "DANs_subtype")
P <- P + coord_flip() + RotatedAxis() + scale_color_gradient2(high = "red", mid = "grey95", low = "blue")
P

theme_set(theme_cowplot())
ModuleNetworkPlot(MB_DANs)

HubGeneNetworkPlot(MB_DANs, n_hubs = 5, n_other = 5, edge_prop = 0.75, mods = 'all', edge.alpha = 0.6)

## Enrichment analysis
mDAN_4_M4.df <- subset(modules, module == "mDAN_4-M4")
mDAN_4_M4.df <- bitr(rownames(mDAN_4_M4.df), fromType = "SYMBOL", toType = c("ENTREZID", "ENSEMBL"), OrgDb = org.Mm.eg.db)

mDAN_4_M4_ID <- unique(mDAN_4_M4.df$ENTREZID)

mDAN_4_M4_ego <- enrichGO(gene = mDAN_4_M4_ID, OrgDb = org.Mm.eg.db, ont = "BP", readable = TRUE)
mDAN_4_M4_ego <- clusterProfiler::simplify(mDAN_4_M4_ego, cutoff = 0.7, by = "p.adjust", select_fun = min)

dotplot(mDAN_4_M4_ego, title = "mDAN_4_M4_ego_BP")
goplot(mDAN_4_M4_ego)

mDAN_4_M4_kegg <- enrichKEGG(gene = mDAN_4_M4_ID, keyType = "kegg", organism = 'mmu', pAdjustMethod = "fdr",
                             pvalueCutoff = 0.5, qvalueCutoff = 0.5)
dotplot(mDAN_4_M4_kegg, title = "mDAN_4_M4_KEGG")


# Sham DANs_subtype and ODCs CellChat
Idents(Mice_Sham_mb) <- "Major_celltype"
Sham_ODCs <- subset(Mice_Sham_mb, idents = "Oligodendrocytes")
Sham_mDANs_ODCs_CC <- merge(Sham_ODCs, y = Mice_DANs.sham, add.cell.ids = c("ODCs", "mDANs"), project = "Sham_mDANs_ODCs_CC",
                   merge.data = TRUE)

Sham_mDANs_ODCs_CC$cell_type <- Sham_mDANs_ODCs_CC@active.ident
data.input <- GetAssayData(Sham_mDANs_ODCs_CC, assay = "RNA", slot = "data")
identity <- subset(Sham_mDANs_ODCs_CC@meta.data, select = "cell_type")

Sham_mDANs_ODCs_CC <- createCellChat(object = data.input, meta = identity, group.by = "cell_type")
levels(Sham_mDANs_ODCs_CC@idents)
groupSize <- as.numeric(table(Sham_mDANs_ODCs_CC@idents))

CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB
Sham_mDANs_ODCs_CC@DB <- CellChatDB.use

Sham_mDANs_ODCs_CC <- subsetData(Sham_mDANs_ODCs_CC)
Sham_mDANs_ODCs_CC <- identifyOverExpressedGenes(Sham_mDANs_ODCs_CC)
Sham_mDANs_ODCs_CC <- identifyOverExpressedInteractions(Sham_mDANs_ODCs_CC)
Sham_mDANs_ODCs_CC <- projectData(Sham_mDANs_ODCs_CC, PPI.mouse)

Sham_mDANs_ODCs_CC <- computeCommunProb(Sham_mDANs_ODCs_CC, raw.use = TRUE)
Sham_mDANs_ODCs_CC <- filterCommunication(Sham_mDANs_ODCs_CC, min.cells = 10)

Sham_mDANs_ODCs_CC <- computeCommunProbPathway(Sham_mDANs_ODCs_CC)
Sham_mDANs_ODCs_CC <- aggregateNet(Sham_mDANs_ODCs_CC)

par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(Sham_mDANs_ODCs_CC@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge = F, title.name = "Number of interactions")

netVisual_circle(Sham_mDANs_ODCs_CC@net$weight, vertex.weight = groupSize, weight.scale = T,
                 label.edge = F, title.name = "Interaction strength")

mat_c <- Sham_mDANs_ODCs_CC@net$count
mat_w <- Sham_mDANs_ODCs_CC@net$weight

par(mfrow = c(2, 3), xpd = TRUE)

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

write.table(mat_c, file = "Sham_mDANs_ODCs_CC_count.txt", col.names = TRUE, sep = "\t", quote = FALSE)
write.table(mat_w, file = "Sham_mDANs_ODCs_CC_weight.txt", col.names = TRUE, sep = "\t", quote = FALSE)

mDANs_subtype_CC <- read.csv(file = "SP table 6 mDANs subtype and ODCs CellChat.csv", header = TRUE)

plot(mDANs_subtype_CC$ODCs.outgoing.count, mDANs_subtype_CC$Proportion.change)
cor(mDANs_subtype_CC$ODCs.outgoing.count, mDANs_subtype_CC$Proportion.change)
cor.test(mDANs_subtype_CC$ODCs.outgoing.count, mDANs_subtype_CC$Proportion.change)
## Pearson's product-moment correlation

## data:  mDANs_subtype_CC$ODCs.outgoing.count and mDANs_subtype_CC$Proportion.change
## t = 3.8403, df = 3, p-value = 0.03114
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##   0.1497088 0.9942298
## sample estimates:
##   cor 
## 0.911572 

result <- lm(Proportion.change~ODCs.outgoing.count, data = mDANs_subtype_CC)
result
summary(result)
plot(mDANs_subtype_CC$ODCs.outgoing.count, mDANs_subtype_CC$Proportion.change, xlim = c(0, 40), ylim = c(-0.8, 0),
     xaxs = "i", yaxs = "i", xaxp = c(0, 40, 4), yaxp = c(-0.8, 0, 4))
abline(result, col = "red", lwd = 2)

plot(mDANs_subtype_CC$ODCs.incoming.count, mDANs_subtype_CC$Proportion.change)
cor(mDANs_subtype_CC$ODCs.incoming.count, mDANs_subtype_CC$Proportion.change)
cor.test(mDANs_subtype_CC$ODCs.incoming.count, mDANs_subtype_CC$Proportion.change)
## Pearson's product-moment correlation

## data:  mDANs_subtype_CC$ODCs.incoming.count and mDANs_subtype_CC$Proportion.change
## t = 6.3707, df = 3, p-value = 0.007828
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##   0.5564160 0.9977724
## sample estimates:
##   cor 
## 0.9649712 

result <- lm(Proportion.change~ODCs.incoming.count, data = mDANs_subtype_CC)
result
summary(result)
plot(mDANs_subtype_CC$ODCs.incoming.count, mDANs_subtype_CC$Proportion.change, xlim = c(0, 40), ylim = c(-0.8, 0),
     xaxs = "i", yaxs = "i", xaxp = c(0, 40, 4), yaxp = c(-0.8, 0, 4))
abline(result, col = "red", lwd = 2)

plot(mDANs_subtype_CC$ODCs.outgoing.weight, mDANs_subtype_CC$Proportion.change)
cor(mDANs_subtype_CC$ODCs.outgoing.weight, mDANs_subtype_CC$Proportion.change)
cor.test(mDANs_subtype_CC$ODCs.outgoing.weight, mDANs_subtype_CC$Proportion.change)
## Pearson's product-moment correlation

## data:  mDANs_subtype_CC$ODCs.outgoing.weight and mDANs_subtype_CC$Proportion.change
## t = 3.6219, df = 3, p-value = 0.0362
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##   0.09742916 0.99358527
## sample estimates:
##   cor 
## 0.9021483  

result <- lm(Proportion.change~ODCs.outgoing.weight, data = mDANs_subtype_CC)
result
summary(result)
plot(mDANs_subtype_CC$ODCs.outgoing.weight, mDANs_subtype_CC$Proportion.change, xlim = c(0, 2), ylim = c(-0.8, 0),
     xaxs = "i", yaxs = "i", xaxp = c(0, 2, 4), yaxp = c(-0.8, 0, 4))
abline(result, col = "red", lwd = 2)

plot(mDANs_subtype_CC$ODCs.incoming.weight, mDANs_subtype_CC$Proportion.change)
cor(mDANs_subtype_CC$ODCs.incoming.weight, mDANs_subtype_CC$Proportion.change)
cor.test(mDANs_subtype_CC$ODCs.incoming.weight, mDANs_subtype_CC$Proportion.change)
## Pearson's product-moment correlation

## data:  mDANs_subtype_CC$ODCs.incoming.weight and mDANs_subtype_CC$Proportion.change
## t = 3.7177, df = 3, p-value = 0.03386
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##   0.1207658 0.9938803
## sample estimates:
##   cor 
## 0.9064519 

result <- lm(Proportion.change~ODCs.incoming.weight, data = mDANs_subtype_CC)
result
summary(result)
plot(mDANs_subtype_CC$ODCs.incoming.weight, mDANs_subtype_CC$Proportion.change, xlim = c(0, 2), ylim = c(-0.8, 0),
     xaxs = "i", yaxs = "i", xaxp = c(0, 2, 4), yaxp = c(-0.8, 0, 4))
abline(result, col = "red", lwd = 2)

Sham_mDANs_ODCs_CC <- netAnalysis_computeCentrality(Sham_mDANs_ODCs_CC, slot.name = "netP")

write.table(Sham_mDANs_ODCs_CC@netP$pathways, file = "pathways.txt", col.names = TRUE, sep = "\t", quote = FALSE)

pathways.show <- "LAMININ"
pathways.show <- "JAM"

netVisual_aggregate(Sham_mDANs_ODCs_CC, signaling = pathways.show, layout = "chord")
netVisual_heatmap(Sham_mDANs_ODCs_CC, signaling = pathways.show, color.heatmap = "Reds")
netAnalysis_contribution(Sham_mDANs_ODCs_CC, signaling = pathways.show)

netVisual_bubble(Sham_mDANs_ODCs_CC, sources.use = 6, targets.use = c(1,2,3,4,5), remove.isolate = FALSE)
netVisual_bubble(Sham_mDANs_ODCs_CC, sources.use = c(1,2,3,4,5), targets.use = 6, remove.isolate = FALSE)

netVisual_chord_gene(Sham_mDANs_ODCs_CC, sources.use = 6, targets.use = c(1,2,3,4,5), signaling = pathways.show, legend.pos.x = 8)
netVisual_chord_gene(Sham_mDANs_ODCs_CC, sources.use = c(1,2,3,4,5), targets.use = 6, signaling = pathways.show, legend.pos.x = 8)

netVisual_chord_gene(Sham_mDANs_ODCs_CC, sources.use = 6, targets.use = c(1,2,3,4,5), slot.name = "netP", legend.pos.x = 10)
netVisual_chord_gene(Sham_mDANs_ODCs_CC, sources.use = c(1,2,3,4,5), targets.use = 6, slot.name = "netP", legend.pos.x = 10)

FeaturePlot(Mice_DANs.sham, features = "Jam3", label = TRUE)


# MPTP DANs_subtype and ODCs CellChat
Idents(Mice_MPTP_mb) <- "Major_celltype"
MPTP_ODCs <- subset(Mice_MPTP_mb, idents = "Oligodendrocytes")
MPTP_mDANs_ODCs_CC <- merge(MPTP_ODCs, y = Mice_DANs.MPTP, add.cell.ids = c("ODCs", "mDANs"), project = "MPTP_mDANs_ODCs_CC",
                      merge.data = TRUE)

MPTP_mDANs_ODCs_CC$cell_type <- MPTP_mDANs_ODCs_CC@active.ident
data.input <- GetAssayData(MPTP_mDANs_ODCs_CC, assay = "RNA", slot = "data")
identity <- subset(MPTP_mDANs_ODCs_CC@meta.data, select = "cell_type")

MPTP_mDANs_ODCs_CC <- createCellChat(object = data.input, meta = identity, group.by = "cell_type")
levels(MPTP_mDANs_ODCs_CC@idents)
groupSize <- as.numeric(table(MPTP_mDANs_ODCs_CC@idents))

MPTP_mDANs_ODCs_CC@DB <- CellChatDB.use

MPTP_mDANs_ODCs_CC <- subsetData(MPTP_mDANs_ODCs_CC)
MPTP_mDANs_ODCs_CC <- identifyOverExpressedGenes(MPTP_mDANs_ODCs_CC)
MPTP_mDANs_ODCs_CC <- identifyOverExpressedInteractions(MPTP_mDANs_ODCs_CC)
MPTP_mDANs_ODCs_CC <- projectData(MPTP_mDANs_ODCs_CC, PPI.mouse)

MPTP_mDANs_ODCs_CC <- computeCommunProb(MPTP_mDANs_ODCs_CC, raw.use = TRUE)
MPTP_mDANs_ODCs_CC <- filterCommunication(MPTP_mDANs_ODCs_CC, min.cells = 10)

MPTP_mDANs_ODCs_CC <- computeCommunProbPathway(MPTP_mDANs_ODCs_CC)
MPTP_mDANs_ODCs_CC <- aggregateNet(MPTP_mDANs_ODCs_CC)

MPTP_mDANs_ODCs_CC <- netAnalysis_computeCentrality(MPTP_mDANs_ODCs_CC, slot.name = "netP")


# CellChat integration
CellChat.list <- list(Sham = Sham_mDANs_ODCs_CC, MPTP = MPTP_mDANs_ODCs_CC)
I_CellChat <- mergeCellChat(CellChat.list, add.names = names(CellChat.list), cell.prefix = TRUE)

I_CellChat
## An object of class CellChat created from a merged object with multiple datasets 
## 828 signaling genes.
## 6694 cells. 
## CellChat analysis of single cell RNA-seq data!

gg1 <- compareInteractions(I_CellChat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(I_CellChat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(I_CellChat, weight.scale = T)
netVisual_diffInteraction(I_CellChat, weight.scale = T, measure = "weight")

gg1 <- netVisual_heatmap(I_CellChat)
gg2 <- netVisual_heatmap(I_CellChat, measure = "weight")
gg1 + gg2

weight.max <- getMaxWeight(CellChat.list, attribute = c("idents","count"))
weight.max <- getMaxWeight(CellChat.list, attribute = c("idents","weight"))

par(mfrow = c(1,2), xpd = TRUE)
for (i in 1:length(CellChat.list)) {
  netVisual_circle(CellChat.list[[i]]@net$weight, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Interaction strength - ", names(CellChat.list)[i]))
}

num.link <- sapply(CellChat.list, function(x) {rowSums(x@net$count) + cODCsums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) 

P1 <- netAnalysis_signalingRole_scatter(CellChat.list$Sham, title = names(CellChat.list$Sham), weight.MinMax = weight.MinMax) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(6,16, by = 2), limits = c(6,16)) + 
  scale_x_continuous(expand = c(0, 0), breaks = seq(6,16, by = 2), limits = c(6,16)) + ggtitle("Sham")
P2 <- netAnalysis_signalingRole_scatter(CellChat.list$MPTP, title = names(CellChat.list$MPTP), weight.MinMax = weight.MinMax) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(6,16, by = 2), limits = c(6,16)) + 
  scale_x_continuous(expand = c(0, 0), breaks = seq(6,16, by = 2), limits = c(6,16)) + ggtitle("MPTP")
P1 + P2

rankNet(I_CellChat, mode = "comparison", stacked = T, do.stat = TRUE)

i = 1
pathway.union <- union(CellChat.list[[i]]@netP$pathways, CellChat.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(CellChat.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(CellChat.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(CellChat.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(CellChat.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(CellChat.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(CellChat.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(CellChat.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(CellChat.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(CellChat.list[[i]], pattern = "all", signaling = pathway.union, title = names(CellChat.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(CellChat.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(CellChat.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

netVisual_bubble(I_CellChat, sources.use = 6, targets.use = c(1:5),  comparison = c(1, 2), angle.x = 45)
netVisual_bubble(I_CellChat, sources.use = c(1:5), targets.use = 6,  comparison = c(1, 2), angle.x = 45)

gg1 <- netVisual_bubble(I_CellChat, sources.use = 6, targets.use = c(1:5),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in MPTP", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(I_CellChat, sources.use = 6, targets.use = c(1:5),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in MPTP", angle.x = 45, remove.isolate = T)
gg1 + gg2

gg1 <- netVisual_bubble(I_CellChat, sources.use = c(1:5), targets.use = 6,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in MPTP", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(I_CellChat, sources.use = c(1:5), targets.use = 6,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in MPTP", angle.x = 45, remove.isolate = T)
gg1 + gg2

## define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "MPTP"
features.name = pos.dataset

I_CellChat <- identifyOverExpressedGenes(I_CellChat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
net <- netMappingDEG(I_CellChat, features.name = features.name)

net.up <- subsetCommunication(I_CellChat, net = net, datasets = "MPTP",ligand.logFC = 0.2, receptor.logFC = NULL)
net.down <- subsetCommunication(I_CellChat, net = net, datasets = "Sham",ligand.logFC = -0.1, receptor.logFC = -0.1)

gene.up <- extractGeneSubsetFromPair(net.up, I_CellChat)
gene.down <- extractGeneSubsetFromPair(net.down, I_CellChat)

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(I_CellChat, pairLR.use = pairLR.use.up, sources.use = 6, targets.use = c(1:5), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(CellChat.list)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(I_CellChat, pairLR.use = pairLR.use.down, sources.use = 6, targets.use = c(1:5), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(CellChat.list)[2]))
gg1 + gg2

gg1 <- netVisual_bubble(I_CellChat, pairLR.use = pairLR.use.up, sources.use = c(1:5), targets.use = 6, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(CellChat.list)[2]))
gg2 <- netVisual_bubble(I_CellChat, pairLR.use = pairLR.use.down, sources.use = c(1:5), targets.use = 6, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(CellChat.list)[2]))
gg1 + gg2

par(mfrow = c(1,2), xpd = TRUE)
netVisual_chord_gene(CellChat.list[[2]], sources.use = c(1:5), targets.use = 6, slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(CellChat.list)[2]))
netVisual_chord_gene(CellChat.list[[1]], sources.use = 6, targets.use = c(1:5), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(CellChat.list)[2]))

pathways.show <- "JAM" 
weight.max <- getMaxWeight(CellChat.list, slot.name = "netP", attribute = pathways.show)

par(mfrow = c(1,2), xpd = TRUE)
for (i in 1:length(CellChat.list)) {
  netVisual_aggregate(CellChat.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(CellChat.list)[i]))
}

ht <- list()
for (i in 1:length(CellChat.list)) {
  ht[[i]] <- netVisual_heatmap(CellChat.list[[i]], signaling = pathways.show, color.heatmap = "Reds", title.name = paste(pathways.show, "signaling ",names(CellChat.list)[i]))
}

ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

for (i in 1:length(CellChat.list)) {
  netVisual_aggregate(CellChat.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(CellChat.list)[i]))
}

for (i in 1:length(CellChat.list)) {
  netVisual_chord_gene(CellChat.list[[i]], sources.use = 6, targets.use = c(4,5), lab.cex = 0.5, title.name = paste0("Signaling from ODC - ", names(CellChat.list)[i]))
}

for (i in 1:length(CellChat.list)) {
  netVisual_chord_gene(CellChat.list[[i]], sources.use = c(1:5), targets.use = 6,  title.name = paste0("Signaling received by ODC ", names(CellChat.list)[i]), legend.pos.x = 10)
}

for (i in 1:length(CellChat.list)) {
  netVisual_chord_gene(CellChat.list[[i]], sources.use = c(1:5), targets.use = 6,slot.name = "netP", title.name = paste0("Signaling pathways sending from mDANs - ", names(CellChat.list)[i]), legend.pos.x = 10)
}

I_CellChat@meta$datasets = factor(I_CellChat@meta$datasets, levels = c("Sham", "MPTP"))
plotGeneExpression(I_CellChat, signaling = "JAM", split.by = "datasets", colors.ggplot = T)


# Jam3 expression
Mice_DANs.combined@meta.data$Jam3 <- Mice_DANs.combined@assays$RNA@counts["Jam3",]
Idents(Mice_DANs.combined) <- "DANs_subtype"
Mice_DANs.combined <- ScaleData(Mice_DANs.combined, features = rownames(Mice_DANs.combined), verbose = FALSE)

VlnPlot(Mice_DANs.combined, features = "Jam3", split.by  = "status", assay = "RNA")

Mice_DANs.combined$subtype.status <- paste(Idents(Mice_DANs.combined), Mice_DANs.combined$status, sep = "_")

Jam3_pos_ratio <- as.data.frame(prop.table(table(Mice_DANs.combined$Jam3 > 0 , Mice_DANs.combined$subtype.status), margin = 2))
Jam3_pos_ratio <- subset(Jam3_pos_ratio, Jam3_pos_ratio$Var1 == "TRUE")
Jam3_pos_ratio <- Jam3_pos_ratio[,-1]
Jam3_pos_ratio$status <- rep(c("MPTP", "Sham"), 5)
Jam3_pos_ratio$DANs_subtype <- c("mDAN_1", "mDAN_1", "mDAN_2", "mDAN_2", "mDAN_3", "mDAN_3", "mDAN_4", "mDAN_4", "mDAN_5", "mDAN_5")

Jam3_pos_ratio$status <- factor(x = Jam3_pos_ratio$status, levels = c("Sham", "MPTP"))
Jam3_pos_ratio$DANs_subtype <- factor(x = Jam3_pos_ratio$DANs_subtype, levels = c("mDAN_1", "mDAN_2", "mDAN_3", "mDAN_4", "mDAN_5"))

ggplot() + geom_bar(data = Jam3_pos_ratio, aes(x = DANs_subtype, y = Freq, fill = status), position = position_dodge2(padding = 0.3), stat = "identity") + 
  theme(axis.title = element_blank(), axis.text.x = element_text(size = 20), axis.ticks.x = element_blank(),
        axis.ticks.length.y = unit(.25, "cm"), axis.line = element_line(color = "black", size = 0.4), 
        axis.text.y = element_text(size = 14), panel.background = element_rect(fill = "white"), legend.position = "top", 
        legend.title = element_blank(), legend.key.height = unit(.3, "cm"), legend.key.width = unit(.8, "cm"), 
        legend.direction = "vertical", legend.spacing.x = unit(.4, "cm"), legend.text = element_text(size = 12)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,0.8, by = 0.16), limits = c(0,0.8))

mDANs_subtype_Jam3_averExp <- AverageExpression(Mice_DANs.sham, slot = "counts", features = "Jam3")$RNA
mDANs_subtype_CC$Jam3_averExp <- mDANs_subtype_Jam3_averExp[1, match(mDANs_subtype_CC$X, colnames(mDANs_subtype_Jam3_averExp))]

plot(mDANs_subtype_CC$Jam3_averExp, mDANs_subtype_CC$Proportion.change)
cor(mDANs_subtype_CC$Jam3_averExp, mDANs_subtype_CC$Proportion.change)
cor.test(mDANs_subtype_CC$Jam3_averExp, mDANs_subtype_CC$Proportion.change)
## Pearson's product-moment correlation

## data:  mDANs_subtype_CC$Jam3_averExp and mDANs_subtype_CC$Proportion.change
## t = 2.035, df = 3, p-value = 0.1347
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##   -0.3679781  0.9832060
## sample estimates:
##   cor 
## 0.7615191 

result <- lm(Proportion.change~Jam3_averExp, data = mDANs_subtype_CC)
result
summary(result)

plot(mDANs_subtype_CC$Jam3_averExp, mDANs_subtype_CC$Proportion.change, xlim = c(0, 1.4), ylim = c(-0.8, 0),
     xaxs = "i", yaxs = "i", xaxp = c(0, 1.4, 4), yaxp = c(-0.8, 0, 4))
abline(result, col = "red", lwd = 2)

Jam3_pos_ratio <- subset(Jam3_pos_ratio, Jam3_pos_ratio$status == "Sham")
mDANs_subtype_CC$Jam3_ExpRatio <- Jam3_pos_ratio[match(mDANs_subtype_CC$X, Jam3_pos_ratio$DANs_subtype), 2]

plot(mDANs_subtype_CC$Jam3_ExpRatio, mDANs_subtype_CC$Proportion.change)
cor(mDANs_subtype_CC$Jam3_ExpRatio, mDANs_subtype_CC$Proportion.change)
cor.test(mDANs_subtype_CC$Jam3_ExpRatio, mDANs_subtype_CC$Proportion.change)
## Pearson's product-moment correlation

## data:  mDANs_subtype_CC$Jam3_ExpRatio and mDANs_subtype_CC$Proportion.change
## t = 2.328, df = 3, p-value = 0.1023
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##   -0.2736966  0.9863715
## sample estimates:
##   cor 
## 0.8023051 

result <- lm(Proportion.change~Jam3_ExpRatio, data = mDANs_subtype_CC)
result
summary(result)

plot(mDANs_subtype_CC$Jam3_ExpRatio, mDANs_subtype_CC$Proportion.change, xlim = c(0, 0.8), ylim = c(-0.8, 0),
     xaxs = "i", yaxs = "i", xaxp = c(0, 0.8, 4), yaxp = c(-0.8, 0, 4))
abline(result, col = "red", lwd = 2)


# VGluT2 expression and glutamate signaling 
Mice_DANs.combined@meta.data$Slc17a6 <- Mice_DANs.combined@assays$RNA@counts["Slc17a6",]

Slc17a6_pos_ratio <- as.data.frame(prop.table(table(Mice_DANs.combined$Slc17a6 >0 , Mice_DANs.combined$subtype.status), margin = 2))
Slc17a6_pos_ratio <- subset(Slc17a6_pos_ratio, Slc17a6_pos_ratio$Var1 == "TRUE")
Slc17a6_pos_ratio <- Slc17a6_pos_ratio[,-1]
Slc17a6_pos_ratio$status <- rep(c("MPTP", "Sham"), 5)
Slc17a6_pos_ratio$subtype <- c("mDAN_1", "mDAN_1", "mDAN_2", "mDAN_2", "mDAN_3", "mDAN_3", "mDAN_4", "mDAN_4", "mDAN_5", "mDAN_5")

Slc17a6_pos_ratio$status <- factor(x = Slc17a6_pos_ratio$status, levels = c("Sham", "MPTP"))
Slc17a6_pos_ratio$subtype <- factor(x = Slc17a6_pos_ratio$subtype, levels = c("mDAN_1", "mDAN_2", "mDAN_3", "mDAN_4", "mDAN_5"))

ggplot() + geom_bar(data = Slc17a6_pos_ratio, aes(x = subtype, y = Freq, fill = status), position = position_dodge2(padding = 0.3), stat = "identity") + 
  theme(axis.title = element_blank(), axis.text.x = element_text(size = 20), axis.ticks.x = element_blank(),
        axis.ticks.length.y = unit(.25, "cm"), axis.line = element_line(color = "black", size = 0.4), 
        axis.text.y = element_text(size = 14), panel.background = element_rect(fill = "white"), legend.position = "top", 
        legend.title = element_blank(), legend.key.height = unit(.3, "cm"), legend.key.width = unit(.8, "cm"), 
        legend.direction = "vertical", legend.spacing.x = unit(.4, "cm"), legend.text = element_text(size = 12)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,0.6, by = 0.15), limits = c(0,0.6))

mDANs_subtype_Slc17a6_averExp <- AverageExpression(Mice_DANs.sham, slot = "counts", features = "Slc17a6")$RNA
mDANs_subtype_CC$Slc17a6_averExp <- mDANs_subtype_Slc17a6_averExp[1, match(mDANs_subtype_CC$X, colnames(mDANs_subtype_Slc17a6_averExp))]

plot(mDANs_subtype_CC$Slc17a6_averExp, mDANs_subtype_CC$Proportion.change)
cor(mDANs_subtype_CC$Slc17a6_averExp, mDANs_subtype_CC$Proportion.change)
cor.test(mDANs_subtype_CC$Slc17a6_averExp, mDANs_subtype_CC$Proportion.change)
## Pearson's product-moment correlation

## data:  mDANs_subtype_CC$Slc17a6_averExp and mDANs_subtype_CC$Proportion.change
## t = -0.21743, df = 3, p-value = 0.8418
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##   -0.9071355  0.8512568
## sample estimates:
##   cor 
## -0.1245543 

result <- lm(Proportion.change~Slc17a6_averExp, data = mDANs_subtype_CC)
result
summary(result)

plot(mDANs_subtype_CC$Slc17a6_averExp, mDANs_subtype_CC$Proportion.change, xlim = c(0, 1), ylim = c(-0.8, 0),
     xaxs = "i", yaxs = "i", xaxp = c(0, 1, 4), yaxp = c(-0.8, 0, 4))
abline(result, col = "red", lwd = 2)

## hdWGCNA
theme_set(theme_cowplot())

resist_DANs <- SetupForWGCNA(Mice_DANs.sham, gene_select = "variable", wgcna_name = "resist_DANs")
length(resist_DANs@misc$resist_DANs$wgcna_genes)

resist_DANs <- MetacellsByGroups(seurat_obj = resist_DANs, group.by = "DANs_subtype", reduction = "pca", k = 15, min_cells = 20, max_shared = 5, ident.group = "DANs_subtype")
resist_DANs <- NormalizeMetacells(resist_DANs)

resist_DANs <- SetDatExpr(resist_DANs, group_name = c("mDAN_1", "mDAN_2"), group.by = "DANs_subtype", assay = "RNA", slot = "data")
resist_DANs <- TestSoftPowers(resist_DANs, networkType = "signed")
plot_list <- PlotSoftPowers(resist_DANs)
wrap_plots(plot_list, ncol = 2)

power_table <- GetPowerTable(resist_DANs)
head(power_table)

## construct co-expression network
resist_DANs <- ConstructNetwork(resist_DANs, soft_power = 7, setDatExpr = FALSE, tom_name = "mDAN_12")
PlotDendrogram(resist_DANs, main = "mDAN_12 hdWGCNA Dendrogram")

resist_DANs@misc$resist_DANs$wgcna_modules %>% head
write.table(resist_DANs@misc$resist_DANs$wgcna_modules, file = "wgcna_modules.txt", col.names = TRUE, sep = "\t", quote = FALSE)

table(resist_DANs@misc$resist_DANs$wgcna_modules$module)
## pink greenyellow     magenta      yellow        blue       green       black      salmon       brown   turquoise         red      purple 
## 98          85          90         115         137         113         100          54         131         201         106          87 
## grey         tan 
## 576          70  

Tom <- GetTOM(resist_DANs)

## Compute harmonized module eigengenes
resist_DANs <- ScaleData(resist_DANs, features =  VariableFeatures(resist_DANs))
resist_DANs <- ModuleEigengenes(resist_DANs, npcs = 20, reduction.use = "pca")

## harmonized module eigengenes
hMEs <- GetMEs(resist_DANs)
head(hMEs)
MEs <- GetMEs(resist_DANs, harmonized = FALSE)
head(MEs)

## Compute module connectivity
resist_DANs <- ModuleConnectivity(resist_DANs, group.by = "DANs_subtype", group_name = c("mDAN_1", "mDAN_2"))
resist_DANs <- ResetModuleNames(resist_DANs, new_name = "mDAN_12-M")

PlotKMEs(resist_DANs, ncol = 5)

modules <- GetModules(resist_DANs)

hub_df <- GetHubGenes(resist_DANs, n_hubs = 10)
head(hub_df)
write.table(hub_df, file = "hub_genes.txt", col.names = TRUE, sep = "\t", quote = FALSE)

resist_DANs <- ModuleExprScore(resist_DANs, n_genes = 25, method = "Seurat")

## Visualization
plot_list <- ModuleFeaturePlot(resist_DANs, features = "hMEs", order = TRUE)
wrap_plots(plot_list, ncol = 4) 

plot_list <- ModuleFeaturePlot(resist_DANs, features = "scores", order = "shuffle", ucell = TRUE)
wrap_plots(plot_list, ncol = 4)

MEs <- GetMEs(resist_DANs, harmonized = TRUE)
mods <- colnames(MEs); mods <- mods[mods != "grey"]

resist_DANs@meta.data <- cbind(resist_DANs@meta.data, MEs)
P <- DotPlot(resist_DANs, features = mods, group.by = "DANs_subtype")
P <- P + coord_flip() + RotatedAxis() + scale_color_gradient2(high = "red", mid = "grey95", low = "blue")
P

theme_set(theme_cowplot())
ModuleNetworkPlot(resist_DANs)

HubGeneNetworkPlot(resist_DANs, n_hubs = 5, n_other = 5, edge_prop = 0.75, mods = 'all', edge.alpha = 0.6)

## Enrichment analysis
M3.df <- subset(modules, module == "mDAN_12-M3")
M3.df <- bitr(rownames(M3.df), fromType = "SYMBOL", toType = c("ENTREZID", "ENSEMBL"), OrgDb = org.Mm.eg.db)
M3_ID <- unique(M3.df$ENTREZID)

M3_GO_ego <- enrichGO(gene = M3_ID, OrgDb = org.Mm.eg.db, ont = "CC", readable = TRUE)
M3_GO_ego <- clusterProfiler::simplify(M3_GO_ego, cutoff = 0.7, by = "p.adjust", select_fun = min)

dotplot(M3_GO_ego, title = "M3_GO_ego_CC")
goplot(M3_GO_ego)

M3_kegg <- enrichKEGG(gene = M3_ID, keyType = "kegg", organism = 'mmu', pAdjustMethod = "fdr", pvalueCutoff = 0.5, qvalueCutoff = 0.5)
dotplot(M3_kegg, title = "M3_KEGG")


# Glutamate receptors expression in ODCs
Mice_ODCs.list <- c(Sham_ODCs, MPTP_ODCs)
features <- SelectIntegrationFeatures(object.list = Mice_ODCs.list)
Mice_ODCs.anchors <- FindIntegrationAnchors(object.list = Mice_ODCs.list, anchor.features = features)
Mice_ODCs.combined <- IntegrateData(anchorset = Mice_ODCs.anchors)

DefaultAssay(Mice_ODCs.combined) <- "integrated"

all.genes <- rownames(Mice_ODCs.combined)
Mice_ODCs.combined <- ScaleData(Mice_ODCs.combined, features = all.genes,verbose = FALSE)
Mice_ODCs.combined <- RunPCA(Mice_ODCs.combined, npcs = 30, verbose = FALSE)

Mice_ODCs.combined <- FindNeighbors(Mice_ODCs.combined, reduction = "pca", dims = 1:30)
Mice_ODCs.combined <- FindClusters(Mice_ODCs.combined, resolution = 0.3)
Mice_ODCs.combined <- RunUMAP(Mice_ODCs.combined, reduction = "pca", dims = 1:30)

P1 <- DimPlot(Mice_ODCs.combined, reduction = "umap", label = TRUE, repel = TRUE)
P2 <- DimPlot(Mice_ODCs.combined, reduction = "umap", group.by = "status")
wrap_plots(plots = list(P1, P2), ncol = 1)

Mice_ODCs.combined <- RenameIdents(Mice_ODCs.combined, `0` = "mODCs_1", `1` = "mODCs_2", `2` = "mODCs_3")
Mice_ODCs.combined$ODCs_subtype <- Mice_ODCs.combined@active.ident
Idents(Mice_ODCs.combined) <- "ODCs_subtype"
DimPlot(Mice_ODCs.combined, reduction = "umap", label = TRUE)

DefaultAssay(Mice_ODCs.combined) <- "RNA"
all.genes <- rownames(Mice_ODCs.combined)
Mice_ODCs.combined <- ScaleData(Mice_ODCs.combined, features = all.genes,verbose = FALSE)

Mice_ODCs.combined@meta.data$status <- factor(x = Mice_ODCs.combined@meta.data$status, levels = c("Sham", "MPTP"))

DotPlot(Mice_ODCs.combined, features = c("Grm1", "Grm2", "Grm3", "Grm4", "Grm5", "Grm6", "Grm7", "Grm8", "Gria1", "Gria2", 
                                        "Gria3", "Gria4", "Grid1", "Grid2", "Grik1", "Grik2", "Grik3", "Grik4", "Grik5", 
                                        "Grin1", "Grin2A", "Grin2b", "Grin2d", "Grin3a", "Grin3b", "Grina"),
        group.by = "status", scale = FALSE) + RotatedAxis()

table(Mice_ODCs.combined$status)
## Sham MPTP 
## 3177 2494

table(Idents(Mice_ODCs.combined), Mice_ODCs.combined$status)
##         Sham  MPTP
## mODCs_1  2065  1580
## mODCs_2  751   606
## mODCs_3  361   308

prop.table(table(Idents(Mice_ODCs.combined), Mice_ODCs.combined$status), margin = 2)
##              Sham       MPTP
## mODCs_1  0.6499843  0.6335204
## mODCs_2  0.2363865  0.2429832
## mODCs_3  0.1136292  0.1234964




# WD ^ ^ written by Shuxuan Lyu 2023/12/05





