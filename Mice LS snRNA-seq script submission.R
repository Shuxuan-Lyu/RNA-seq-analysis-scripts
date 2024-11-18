library(DoubletFinder)
library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ComplexHeatmap)
library(stringr)
library(RColorBrewer)
library(clusterProfiler)
library(org.Mm.eg.db)
library(plot1cell)
library(cowplot)
library(WGCNA)
library(hdWGCNA)
library(patchwork)
library(igraph)
library(enrichR)
library(GeneOverlap)
library(ggsci)
library(psych)
library(DOSE)
library(GOSemSim)
library(enrichplot)


rm(list = ls())
options(stringsAsFactors = F)


# Male_dataset
male_1.data <- Read10X(data.dir = "data/male/c57-02.matrix/")
male_2.data <- Read10X(data.dir = "data/male/c57-03.matrix/")

## Remove doublets
## male_1
male_1 <- CreateSeuratObject(counts = male_1.data, project = "male_1_LS", min.cells = 3)
male_1@meta.data$sample_ID <- c(rep("male_1", ncol(male_1.data)))
male_1
## An object of class Seurat 
## 22224 features across 9050 samples within 1 assay 
## Active assay: RNA (22224 features, 0 variable features)

male_1 <- NormalizeData(male_1, normalization.method = "LogNormalize", scale.factor = 10000)
male_1 <- FindVariableFeatures(male_1, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(male_1)
male_1 <- ScaleData(male_1, features = all.genes, verbose = FALSE)
male_1 <- RunPCA(male_1, npcs = 30, verbose = FALSE)

male_1 <- FindNeighbors(male_1, reduction = "pca", dims = 1:20)
male_1 <- FindClusters(male_1, resolution = 0.3)
male_1 <- RunUMAP(male_1, dims = 1:20)
male_1_p1 <- DimPlot(male_1, reduction = "umap", label = TRUE)
male_1
## An object of class Seurat 
## 22224 features across 9050 samples within 1 assay 
## Active assay: RNA (22224 features, 2000 variable features)
## 2 dimensional reductions calculated: pca, umap

sweep.res.list_male_1 <- paramSweep_v3(male_1, PCs = 1:20, sct = FALSE)
sweep.stats_male_1 <- summarizeSweep(sweep.res.list_male_1, GT = FALSE)
bcmvn_male_1 <- find.pK(sweep.stats_male_1)
opt_pK_male_1 <- as.numeric(as.vector(bcmvn_male_1$pK[which.max(bcmvn_male_1$BCmetric)]))
print(opt_pK_male_1)
## [1] 0.23

annotations_male_1 <- male_1@meta.data$seurat_clusters
homotypic.prop_male_1 <- modelHomotypic(annotations_male_1)
nExp_poi_male_1 <- round(0.069*nrow(male_1@meta.data))
nExp_poi.adj_male_1 <- round(nExp_poi_male_1*(1-homotypic.prop_male_1))

male_1 <- doubletFinder_v3(male_1, PCs = 1:20, pN = 0.25, pK = opt_pK_male_1, nExp = nExp_poi.adj_male_1, 
                           reuse.pANN = FALSE, sct = FALSE)  
male_1_p2 <- DimPlot(male_1, reduction = "umap", group.by = "DF.classifications_0.25_0.23_570")
CombinePlots(plots = list(male_1_p1, male_1_p2), ncol = 1)

table(male_1@meta.data$DF.classifications_0.25_0.23_570)
## Doublet Singlet 
## 570     8480

Idents(male_1) <- "DF.classifications_0.25_0.23_570"
table(male_1@active.ident)

male_1_RemDoub <- subset(male_1, idents = "Singlet")
male_1_RemDoub
## An object of class Seurat 
## 22224 features across 8480 samples within 1 assay 
## Active assay: RNA (22224 features, 2000 variable features)
## 2 dimensional reductions calculated: pca, umap

Idents(male_1_RemDoub) <- "seurat_clusters"
DimPlot(male_1_RemDoub, reduction = "umap", label = TRUE)

## male_2
male_2 <- CreateSeuratObject(counts = male_2.data, project = "male_2_LS", min.cells = 3)
male_2@meta.data$sample_ID <- c(rep("male_2", ncol(male_2.data)))
male_2
## An object of class Seurat 
## 21490 features across 8204 samples within 1 assay 
## Active assay: RNA (21490 features, 0 variable features)

male_2 <- NormalizeData(male_2, normalization.method = "LogNormalize", scale.factor = 10000)
male_2 <- FindVariableFeatures(male_2, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(male_2)
male_2 <- ScaleData(male_2, features = all.genes, verbose = FALSE)
male_2 <- RunPCA(male_2, npcs = 30, verbose = FALSE)

male_2 <- FindNeighbors(male_2, reduction = "pca", dims = 1:20)
male_2 <- FindClusters(male_2, resolution = 0.3)
male_2 <- RunUMAP(male_2, dims = 1:20)
male_2_p1 <- DimPlot(male_2, reduction = "umap", label = TRUE)
male_2
## An object of class Seurat 
## 21490 features across 8204 samples within 1 assay 
## Active assay: RNA (21490 features, 2000 variable features)
## 2 dimensional reductions calculated: pca, umap

sweep.res.list_male_2 <- paramSweep_v3(male_2, PCs = 1:20, sct = FALSE)
sweep.stats_male_2 <- summarizeSweep(sweep.res.list_male_2, GT = FALSE)
bcmvn_male_2 <- find.pK(sweep.stats_male_2)
opt_pK_male_2 <- as.numeric(as.vector(bcmvn_male_2$pK[which.max(bcmvn_male_2$BCmetric)]))
print(opt_pK_male_2)
# [1] 0.1

annotations_male_2 <- male_2@meta.data$seurat_clusters
homotypic.prop_male_2 <- modelHomotypic(annotations_male_2)
nExp_poi_male_2 <- round(0.061*nrow(male_2@meta.data))
nExp_poi.adj_male_2 <- round(nExp_poi_male_2*(1-homotypic.prop_male_2))

male_2 <- doubletFinder_v3(male_2, PCs = 1:20, pN = 0.25, pK = opt_pK_male_2, nExp = nExp_poi.adj_male_2, 
                           reuse.pANN = FALSE, sct = FALSE)  
male_2_p2 <- DimPlot(male_2, reduction = "umap", group.by = "DF.classifications_0.25_0.1_449")
CombinePlots(plots = list(male_2_p1, male_2_p2), ncol = 1)

table(male_2@meta.data$DF.classifications_0.25_0.1_449)
## Doublet Singlet 
## 449     7755

Idents(male_2) <- "DF.classifications_0.25_0.1_449"
table(male_2@active.ident)

male_2_RemDoub <- subset(male_2, idents = "Singlet")
male_2_RemDoub
## An object of class Seurat 
## 21490 features across 7755 samples within 1 assay 
## Active assay: RNA (21490 features, 2000 variable features)
## 2 dimensional reductions calculated: pca, umap

Idents(male_2_RemDoub) <- "seurat_clusters"
DimPlot(male_2_RemDoub, reduction = "umap", label = TRUE)


## Initialize the Seurat object with the raw (non-normalized data)
male_1_RemDoub.data <- as.data.frame(male_1_RemDoub@assays$RNA@counts)
male_1_RemDoub.data$X <- rownames(male_1_RemDoub.data)
male_2_RemDoub.data <- as.data.frame(male_2_RemDoub@assays$RNA@counts)
male_2_RemDoub.data$X <- rownames(male_2_RemDoub.data)

LS_male.data <- merge(male_1_RemDoub.data, male_2_RemDoub.data, by = "X", all = TRUE)
rownames(LS_male.data) <- LS_male.data$X
LS_male.data <- LS_male.data[, -1]
LS_male.data[is.na(LS_male.data)] <- 0

LS_male <- CreateSeuratObject(counts = LS_male.data, project = "LS_male", min.cells = 3)
LS_male@meta.data$sample_ID <- c(rep("male_1", ncol(male_1_RemDoub)), rep("male_2", ncol(male_2_RemDoub)))

LS_male
## An object of class Seurat
## 22295 features across 16235 samples within 1 assay 
## Active assay: RNA (22295 features, 0 variable features)

table(LS_male$sample_ID)
## male_1 male_2 
## 8480   7755 

## perform QC and treat data for PCA
LS_male[["percent.mt"]] <- PercentageFeatureSet(LS_male, pattern = "^mt-")
VlnPlot(LS_male, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(LS_male, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, split.by = "sample_ID", 
        pt.size = .1)
LS_male <- subset(LS_male, subset = nFeature_RNA > 200 & percent.mt < 5)
LS_male <- subset(LS_male, subset = nFeature_RNA < 6000)

LS_male <- NormalizeData(LS_male, normalization.method = "LogNormalize", scale.factor = 10000)
LS_male <- FindVariableFeatures(LS_male, selection.method = "vst", nfeatures = 3000)
LS_male
## An object of class Seurat 
## 22295 features across 16184 samples within 1 assay 
## Active assay: RNA (22295 features, 3000 variable features)



# Female_dataset
Female_1.data <- Read10X(data.dir = "data/female/LJL_0928.matrix/")
Female_2.data <- Read10X(data.dir = "data/female/LJL_1012.matrix/")

## Remove doublets
## Female_1
Female_1 <- CreateSeuratObject(counts = Female_1.data, project = "Female_1_LS", min.cells = 3)
Female_1@meta.data$sample_ID <- c(rep("Female_1", ncol(Female_1.data)))
Female_1
## An object of class Seurat 
## 22385 features across 7759 samples within 1 assay 
## Active assay: RNA (22385 features, 0 variable features)

Female_1 <- NormalizeData(Female_1, normalization.method = "LogNormalize", scale.factor = 10000)
Female_1 <- FindVariableFeatures(Female_1, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(Female_1)
Female_1 <- ScaleData(Female_1, features = all.genes, verbose = FALSE)
Female_1 <- RunPCA(Female_1, npcs = 30, verbose = FALSE)

Female_1 <- FindNeighbors(Female_1, reduction = "pca", dims = 1:20)
Female_1 <- FindClusters(Female_1, resolution = 0.3)
Female_1 <- RunUMAP(Female_1, dims = 1:20)
Female_1_p1 <- DimPlot(Female_1, reduction = "umap", label = TRUE)
Female_1
## An object of class Seurat 
## 22385 features across 7759 samples within 1 assay 
## Active assay: RNA (22385 features, 2000 variable features)
## 2 dimensional reductions calculated: pca, umap

sweep.res.list_Female_1 <- paramSweep_v3(Female_1, PCs = 1:20, sct = FALSE)
sweep.stats_Female_1 <- summarizeSweep(sweep.res.list_Female_1, GT = FALSE)
bcmvn_Female_1 <- find.pK(sweep.stats_Female_1)
opt_pK_Female_1 <- as.numeric(as.vector(bcmvn_Female_1$pK[which.max(bcmvn_Female_1$BCmetric)]))
print(opt_pK_Female_1)
## [1] 0.2

annotations_Female_1 <- Female_1@meta.data$seurat_clusters
homotypic.prop_Female_1 <- modelHomotypic(annotations_Female_1)
nExp_poi_Female_1 <- round(0.061*nrow(Female_1@meta.data))
nExp_poi.adj_Female_1 <- round(nExp_poi_Female_1*(1-homotypic.prop_Female_1))

Female_1 <- doubletFinder_v3(Female_1, PCs = 1:20, pN = 0.25, pK = opt_pK_Female_1, nExp = nExp_poi.adj_Female_1, 
                             reuse.pANN = FALSE, sct = FALSE)  
Female_1_p2 <- DimPlot(Female_1, reduction = "umap", group.by = "DF.classifications_0.25_0.2_431")
CombinePlots(plots = list(Female_1_p1, Female_1_p2), ncol = 1)

table(Female_1@meta.data$DF.classifications_0.25_0.2_431)
## Doublet Singlet 
## 431     7328

Idents(Female_1) <- "DF.classifications_0.25_0.2_431"
table(Female_1@active.ident)

Female_1_RemDoub <- subset(Female_1, idents = "Singlet")
Female_1_RemDoub
## An object of class Seurat 
## 22385 features across 7328 samples within 1 assay 
## Active assay: RNA (22385 features, 2000 variable features)
## 2 dimensional reductions calculated: pca, umap

Idents(Female_1_RemDoub) <- "seurat_clusters"
DimPlot(Female_1_RemDoub, reduction = "umap", label = TRUE)

## Female_2
Female_2 <- CreateSeuratObject(counts = Female_2.data, project = "Female_2_LS", min.cells = 3)
Female_2@meta.data$sample_ID <- c(rep("Female_2", ncol(Female_2.data)))
Female_2
## An object of class Seurat 
## 22066 features across 8191 samples within 1 assay 
## Active assay: RNA (22066 features, 0 variable features)

Female_2 <- NormalizeData(Female_2, normalization.method = "LogNormalize", scale.factor = 10000)
Female_2 <- FindVariableFeatures(Female_2, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(Female_2)
Female_2 <- ScaleData(Female_2, features = all.genes, verbose = FALSE)
Female_2 <- RunPCA(Female_2, npcs = 30, verbose = FALSE)

Female_2 <- FindNeighbors(Female_2, reduction = "pca", dims = 1:20)
Female_2 <- FindClusters(Female_2, resolution = 0.3)
Female_2 <- RunUMAP(Female_2, dims = 1:20)
Female_2_p1 <- DimPlot(Female_2, reduction = "umap", label = TRUE)
Female_2
## An object of class Seurat 
## 22066 features across 8191 samples within 1 assay 
## Active assay: RNA (22066 features, 2000 variable features)
## 2 dimensional reductions calculated: pca, umap

sweep.res.list_Female_2 <- paramSweep_v3(Female_2, PCs = 1:20, sct = FALSE)
sweep.stats_Female_2 <- summarizeSweep(sweep.res.list_Female_2, GT = FALSE)
bcmvn_Female_2 <- find.pK(sweep.stats_Female_2)
opt_pK_Female_2 <- as.numeric(as.vector(bcmvn_Female_2$pK[which.max(bcmvn_Female_2$BCmetric)]))
print(opt_pK_Female_2)
# [1] 0.21

annotations_Female_2 <- Female_2@meta.data$seurat_clusters
homotypic.prop_Female_2 <- modelHomotypic(annotations_Female_2)
nExp_poi_Female_2 <- round(0.061*nrow(Female_2@meta.data))
nExp_poi.adj_Female_2 <- round(nExp_poi_Female_2*(1-homotypic.prop_Female_2))

Female_2 <- doubletFinder_v3(Female_2, PCs = 1:20, pN = 0.25, pK = opt_pK_Female_2, nExp = nExp_poi.adj_Female_2, 
                             reuse.pANN = FALSE, sct = FALSE)  
Female_2_p2 <- DimPlot(Female_2, reduction = "umap", group.by = "DF.classifications_0.25_0.21_453")
CombinePlots(plots = list(Female_2_p1, Female_2_p2), ncol = 1)

table(Female_2@meta.data$DF.classifications_0.25_0.21_453)
## Doublet Singlet 
## 453     7738

Idents(Female_2) <- "DF.classifications_0.25_0.21_453"
table(Female_2@active.ident)

Female_2_RemDoub <- subset(Female_2, idents = "Singlet")
Female_2_RemDoub
## An object of class Seurat 
## 22066 features across 7738 samples within 1 assay 
## Active assay: RNA (22066 features, 2000 variable features)
## 2 dimensional reductions calculated: pca, umap

Idents(Female_2_RemDoub) <- "seurat_clusters"
DimPlot(Female_2_RemDoub, reduction = "umap", label = TRUE)


## Initialize the Seurat object with the raw (non-normalized data)
Female_1_RemDoub.data <- as.data.frame(Female_1_RemDoub@assays$RNA@counts)
Female_1_RemDoub.data$X <- rownames(Female_1_RemDoub.data)
Female_2_RemDoub.data <- as.data.frame(Female_2_RemDoub@assays$RNA@counts)
Female_2_RemDoub.data$X <- rownames(Female_2_RemDoub.data)

LS_Female.data <- merge(Female_1_RemDoub.data, Female_2_RemDoub.data, by = "X", all = TRUE)
rownames(LS_Female.data) <- LS_Female.data$X
LS_Female.data <- LS_Female.data[, -1]
LS_Female.data[is.na(LS_Female.data)] <- 0

LS_Female <- CreateSeuratObject(counts = LS_Female.data, project = "LS_Female", min.cells = 3)
LS_Female@meta.data$sample_ID <- c(rep("Female_1", ncol(Female_1_RemDoub)), rep("Female_2", ncol(Female_2_RemDoub)))

LS_Female
## An object of class Seurat 
## 22637 features across 15066 samples within 1 assay 
## Active assay: RNA (22637 features, 0 variable features)

table(LS_Female$sample_ID)
## Female_1  Female_2 
## 7328      7738  

## perform QC and treat data for PCA
LS_Female[["percent.mt"]] <- PercentageFeatureSet(LS_Female, pattern = "^mt-")
VlnPlot(LS_Female, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(LS_Female, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, split.by = "sample_ID", 
        pt.size = .1)
LS_Female <- subset(LS_Female, subset = nFeature_RNA > 200 & percent.mt < 5)
LS_Female <- subset(LS_Female, subset = nFeature_RNA < 6000)

LS_Female <- NormalizeData(LS_Female, normalization.method = "LogNormalize", scale.factor = 10000)
LS_Female <- FindVariableFeatures(LS_Female, selection.method = "vst", nfeatures = 3000)
LS_Female
## An object of class Seurat 
## 22637 features across 14983 samples within 1 assay 
## Active assay: RNA (22637 features, 3000 variable features)



# Perform integration
LS.list = c(LS_male, LS_Female)
features <- SelectIntegrationFeatures(object.list = LS.list, nfeatures = 4000)

LS.anchors <- FindIntegrationAnchors(object.list = LS.list, anchor.features = features)
LS.combined <- IntegrateData(anchorset = LS.anchors)

DefaultAssay(LS.combined) <- "integrated"

all.genes <- rownames(LS.combined)
LS.combined <- ScaleData(LS.combined, features = all.genes, verbose = FALSE, vars.to.regress = "percent.mt")
LS.combined <- RunPCA(LS.combined, npcs = 50, verbose = FALSE)
LS.combined
## An object of class Seurat 
## 26608 features across 31167 samples within 2 assays 
## Active assay: integrated (3516 features, 3516 variable features)
## 1 other assay present: RNA
## 1 dimensional reduction calculated: pca

table(LS.combined$sample_ID)
## Female_1  Female_2   male_1   male_2 
## 7266      7717       8472     7712

LS.combined$gender <- LS.combined@active.ident
table(LS.combined$gender)
## LS_Female   LS_male 
## 14983       16184

LS.combined$sample_ID <- factor(x = LS.combined$sample_ID, levels = c("male_1", "male_2", "Female_1", "Female_2"))
LS.combined$gender <- factor(x = LS.combined$gender, levels = c("LS_male", "LS_Female"))

VizDimLoadings(LS.combined, dims = 1:2, reduction = "pca")
DimPlot(object = LS.combined, reduction = "pca", pt.size = .1, group.by = "sample_ID")
VlnPlot(object = LS.combined, features = "PC_1", group.by = "sample_ID", pt.size = .1)

plot1 <- FeatureScatter(LS.combined, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "sample_ID")
plot2 <- FeatureScatter(LS.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "sample_ID")
CombinePlots(plots = list(plot1, plot2))

## Run Harmony
LS.combined <- RunHarmony(LS.combined, group.by.vars = "sample_ID", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(LS.combined, 'harmony')
harmony_embeddings[1:5, 1:5]

DimPlot(object = LS.combined, reduction = "harmony", pt.size = .1, group.by = "sample_ID")
VlnPlot(object = LS.combined, features = "harmony_1", group.by = "sample_ID", pt.size = .1)

LS.combined
## An object of class Seurat 
## 26608 features across 31167 samples within 2 assays 
## Active assay: integrated (3516 features, 3516 variable features)
## 1 other assay present: RNA
## 2 dimensional reductions calculated: pca, harmony

## U-MAP and Clustering with harmony
LS.combined <- FindNeighbors(LS.combined, reduction = "harmony", dims = 1:30)
LS.combined <- FindClusters(LS.combined, resolution = 0.2)
LS.combined <- RunUMAP(LS.combined, reduction = "harmony", dims = 1:30)

DimPlot(LS.combined, reduction = "umap", label = TRUE)
DimPlot(LS.combined, reduction = "umap", group.by = "gender")
DimPlot(LS.combined, reduction = "umap", group.by = "sample_ID")

## find markers for every cluster compared to all remaining cells, report only the positive ones
DefaultAssay(LS.combined) <- "RNA"
all.genes <- rownames(LS.combined)
LS.combined <- ScaleData(LS.combined, features = all.genes,verbose = FALSE, vars.to.regress = "percent.mt")

LS.combined.allmarkers <- FindAllMarkers(LS.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(LS.combined.allmarkers, file = "LS.AllMarkers.txt", col.names = TRUE, sep = "\t", quote = FALSE)

## annotate each cluster based on canonical marker genes
FeaturePlot(LS.combined, features = c("Aqp4", "Gfap", "Slc1a3", "Slc1a2", "Slc4a4", "Ntsr2"), label = TRUE)
VlnPlot(LS.combined, features = c("Aqp4", "Gfap", "Slc1a3", "Slc1a2", "Slc4a4", "Ntsr2"), ncol = 1, pt.size = 0.1)

FeaturePlot(LS.combined, features = c("Cx3cr1", "P2ry12", "Csf1r"), label = TRUE)
VlnPlot(LS.combined, features = c("Cx3cr1", "P2ry12", "Csf1r"), ncol = 1, pt.size = 0.1)

FeaturePlot(LS.combined, features = c("Opalin", "Mog", "Plp1", "Mbp"), label = TRUE)
VlnPlot(LS.combined, features = c("Opalin", "Mog", "Plp1", "Mbp"), ncol = 1, pt.size = 0.1)

FeaturePlot(LS.combined, features = c("Vcan", "Bcan", "Pdgfra", "Olig2"), label = TRUE)
VlnPlot(LS.combined, features = c("Vcan", "Bcan", "Pdgfra", "Olig2"), ncol = 1, pt.size = 0.1)

FeaturePlot(LS.combined, features = c("Flt1", "Cldn5", "Ptprb"), label = TRUE)
VlnPlot(LS.combined, features = c("Flt1", "Cldn5", "Ptprb"), ncol = 1, pt.size = 0.1)

FeaturePlot(LS.combined, features = c("Hdc", "Fam216b", "Foxj1", "Dnah11", "Dnah12", "Spag16", "Ccdc153"), label = TRUE)
VlnPlot(LS.combined, features = c("Hdc", "Fam216b", "Foxj1", "Dnah11", "Dnah12", "Spag16", "Ccdc153"), ncol = 1, pt.size = 0.1)

FeaturePlot(LS.combined, features = c("Slc17a6", "Slc32a1", "Gad1", "Gad2", "Eno2", "Rbfox3", "Snap25", "Syt1"), label = TRUE)
VlnPlot(LS.combined, features = c("Slc32a1", "Gad1", "Gad2", "Eno2", "Rbfox3", "Snap25", "Syt1"), ncol = 1, pt.size = 0.1)

VlnPlot(LS.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), stack = TRUE, flip = TRUE)
table(LS.combined$seurat_clusters)

LS.combined <- RenameIdents(LS.combined, `0` = "Astrocytes", `1` = "GABA Neu", `2` = "GABA Neu", `3` = "GABA Neu", `4` = "GABA Neu",
                            `5` = "GABA Neu", `6` = "GABA Neu", `7` = "Ependymal cells", `8` = "GABA Neu", `9` = "diff.Neu",
                            `10` = "Oligodendrocytes", `11` = "GABA Neu", `12` = "OPCs", `13` = "Microglia", `14` = "Immature Neu",
                            `15` = "GABA Neu", `16` = "Endothelial cells", `17` = "GABA Neu", `18` = "OPCs", `19` = "Oligodendrocytes")

LS.combined$Major_celltype <- LS.combined@active.ident
Idents(LS.combined) <- "Major_celltype"

Majorcelltype_level <- c("GABA Neu", "diff.Neu", "Immature Neu", "Astrocytes", "Oligodendrocytes", "OPCs", "Microglia", "Ependymal cells",
                         "Endothelial cells")
LS.combined$Major_celltype <- factor(x = LS.combined$Major_celltype, levels = Majorcelltype_level)
Idents(LS.combined) <- "Major_celltype"
DimPlot(LS.combined, reduction = "umap", label = TRUE)

VlnPlot(LS.combined, features = c("Gad1", "Gad2", "Rbfox3", "Lars2", "Adarb2", "Ntsr2", "Mag", "Pdgfra", "P2ry12", "Dnah11", "Flt1"),
        stack = TRUE, same.y.lims = TRUE, flip = TRUE)

## Prepare data for integrated UMAP plotting
circ_data <- prepare_circlize_data(LS.combined, scale = 0.8)
set.seed(123456)

Major_celltype_colors <- rand_color(length(levels(LS.combined)))
gender_colors <- rand_color(length(names(table(LS.combined$gender))))
cluster_colors <- rand_color(length(names(table(LS.combined$seurat_clusters))))

plot_circlize(circ_data, do.label = T, pt.size = 0.2, col.use = Major_celltype_colors, bg.color = 'white', kde2d.n = 1000,
              repel = T, label.cex = 1)

add_track(circ_data, group = "gender", colors = gender_colors, track_num = 2)
add_track(circ_data, group = "seurat_clusters", colors = cluster_colors, track_num = 3)


# LS all Neurons re-clustering
LS_Neurons.data <- subset(LS.combined, idents = "GABA Neu")
LS_Neurons <- CreateSeuratObject(counts = LS_Neurons.data@assays$RNA@counts, project = "LS_Neurons", min.cells = 3)
LS_Neurons@meta.data$sample_ID <- LS_Neurons.data$sample_ID
LS_Neurons@meta.data$gender <- LS_Neurons.data$gender
LS_Neurons
## An object of class Seurat 
## 22450 features across 18715 samples within 1 assay 
## Active assay: RNA (22450 features, 0 variable features)

table(LS_Neurons$sample_ID)
## male_1   male_2   Female_1   Female_2 
## 5237     4535     4353       4590

table(LS_Neurons$gender)
## LS_male   LS_Female 
## 9772      8943 

LS_Neurons[["percent.mt"]] <- PercentageFeatureSet(LS_Neurons, pattern = "^mt-")
VlnPlot(LS_Neurons, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(LS_Neurons, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, split.by = "sample_ID",pt.size = .1)

LS_Neurons <- subset(LS_Neurons, subset = nFeature_RNA > 1000 & percent.mt < 2)
LS_Neurons
## An object of class Seurat 
## 22450 features across 18693 samples within 1 assay 
## Active assay: RNA (22450 features, 0 variable features)

## remove mt- genes
LS_Neurons <- LS_Neurons[-c(grep("^mt-", rownames(LS_Neurons))),]
grep("^mt-", rownames(LS_Neurons))
LS_Neurons
## An object of class Seurat 
## 22437 features across 18693 samples within 1 assay 
## Active assay: RNA (22437 features, 0 variable features)

LS_Neurons[["percent.mt"]] <- PercentageFeatureSet(LS_Neurons, pattern = "^mt-")
VlnPlot(LS_Neurons, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

Idents(LS_Neurons) <- "gender"
LS_Neurons_male <- subset(LS_Neurons, idents = "LS_male")
LS_Neurons_Female <- subset(LS_Neurons, idents = "LS_Female")

LS_Neurons_male <- NormalizeData(LS_Neurons_male, normalization.method = "LogNormalize", scale.factor = 10000)
LS_Neurons_male <- FindVariableFeatures(LS_Neurons_male, selection.method = "vst", nfeatures = 3000)

LS_Neurons_Female <- NormalizeData(LS_Neurons_Female, normalization.method = "LogNormalize", scale.factor = 10000)
LS_Neurons_Female <- FindVariableFeatures(LS_Neurons_Female, selection.method = "vst", nfeatures = 3000)

LS_Neurons.list = c(LS_Neurons_male, LS_Neurons_Female)
features <- SelectIntegrationFeatures(object.list = LS_Neurons.list, nfeatures = 3000)

LS_Neurons.anchors <- FindIntegrationAnchors(object.list = LS_Neurons.list, anchor.features = features)
LS_Neurons.combined <- IntegrateData(anchorset = LS_Neurons.anchors)

DefaultAssay(LS_Neurons.combined) <- "integrated"
LS_Neurons.combined@meta.data$gender <- factor(x = LS_Neurons.combined@meta.data$gender, levels = c("LS_male", "LS_Female"))

all.genes <- rownames(LS_Neurons.combined)
LS_Neurons.combined <- ScaleData(LS_Neurons.combined, features = all.genes, verbose = FALSE)
LS_Neurons.combined <- RunPCA(LS_Neurons.combined, npcs = 50, verbose = FALSE)
LS_Neurons.combined
## An object of class Seurat 
## 25437 features across 18693 samples within 2 assays 
## Active assay: integrated (3000 features, 3000 variable features)
## 1 other assay present: RNA
## 1 dimensional reduction calculated: pca

table(LS_Neurons.combined$sample_ID)
## Female_1 Female_2   male_1   male_2 
## 4352     4584     5231     4526 

table(LS_Neurons.combined$gender)
## LS_male   LS_Female 
## 9757      8936

VizDimLoadings(LS_Neurons.combined, dims = 1:2, reduction = "pca")
DimPlot(object = LS_Neurons.combined, reduction = "pca", pt.size = .1, group.by = "sample_ID")
VlnPlot(object = LS_Neurons.combined, features = "PC_1", group.by = "sample_ID", pt.size = .1)

plot1 <- FeatureScatter(LS_Neurons.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(LS_Neurons.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

## Run Harmony
LS_Neurons.combined <- RunHarmony(LS_Neurons.combined, group.by.vars = "sample_ID", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(LS_Neurons.combined, 'harmony')
harmony_embeddings[1:5, 1:5]

DimPlot(object = LS_Neurons.combined, reduction = "harmony", pt.size = .1, group.by = "sample_ID")
VlnPlot(object = LS_Neurons.combined, features = "harmony_1", group.by = "sample_ID", pt.size = .1)

LS_Neurons.combined
## An object of class Seurat 
## 25437 features across 18693 samples within 2 assays 
## Active assay: integrated (3000 features, 3000 variable features)
## 1 other assay present: RNA
## 2 dimensional reductions calculated: pca, harmony

## U-MAP and Clustering with harmony
LS_Neurons.combined <- FindNeighbors(LS_Neurons.combined, reduction = "harmony", dims = 1:20)
LS_Neurons.combined <- FindClusters(LS_Neurons.combined, resolution = 0.2)
LS_Neurons.combined <- RunUMAP(LS_Neurons.combined, reduction = "harmony", dims = 1:20)

DimPlot(LS_Neurons.combined, reduction = "umap", label = TRUE)
DimPlot(LS_Neurons.combined, reduction = "umap", group.by = "gender")

## find markers for every neuronal subtype
LS_Neu.markers <- FindAllMarkers(LS_Neurons.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(LS_Neu.markers, file = "LS_Neu.AllMarkers.txt", col.names = TRUE, sep = "\t", quote = FALSE)

LS_Neurons.combined <- RenameIdents(LS_Neurons.combined, `0` = "N1", `1` = "N2", `2` = "N3", `3` = "N4", `4` = "N5", `5` = "N6", `6` = "N7",
                                    `7` = "N8", `8` = "N9", `9` = "N10", `10` = "N11", `11` = "N12", `12` = "N13", `13` = "N14", `14` = "N15")

LS_Neurons.combined$Neu_subtype <- LS_Neurons.combined@active.ident
DimPlot(LS_Neurons.combined, reduction = "umap", label = TRUE)

table(LS_Neurons.combined$Neu_subtype)
## N1   N2   N3   N4   N5   N6   N7   N8    N9   N10  N11  N12  N13   N14  N15 
## 2440 2257 2158 2107 2078 1913 1779 1455  788  516  452  343  261   98   48 

DefaultAssay(LS_Neurons.combined) <- "RNA"
all.genes <- rownames(LS_Neurons.combined)
LS_Neurons.combined <- ScaleData(LS_Neurons.combined, features = all.genes,verbose = FALSE)

top5NeuMarkers <- LS_Neu.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5Neufeatures <- unique(top5NeuMarkers$gene)

DoHeatmap(LS_Neurons.combined, features = top5Neufeatures, size = 2, draw.lines = FALSE, angle = 45, slot = "scale.data", hjust = 0.2) +
  theme(axis.text.y = element_text(size = 4)) + scale_fill_gradientn(colors = colorRampPalette(brewer.pal(n = 3, name = "RdBu"))(100))

VlnPlot(LS_Neurons.combined, features = c("Scn1a", "Scn2a", "Scn3a", "Scn4a", "rna_Scn5a", "Scn7a", "Scn8a", "Scn9a",
                                          "Scn10a", "Scn11a"), 
        stack = TRUE, same.y.lims = TRUE, flip = TRUE, fill.by = "ident") + NoLegend()

## Scn5a and coexpression
FeaturePlot(LS_Neurons.combined, feature = "Scn5a",label = TRUE)
VlnPlot(LS_Neurons.combined, feature = "Scn5a", split.by = "gender")

complex_vlnplot_single(LS_Neurons.combined, feature = "Scn5a", groups = "gender", alpha = 1)
complex_vlnplot_multiple(LS_Neurons.combined, features = c("Scn5a", "Nts", "Crhr1", "Adora2a", "Oxtr", "Avpr1a", "Ntrk2"),
                         group = "gender", add.dot = TRUE, alpha = 0.01, font.size = 10, pt.size = 0.01)

FeaturePlot(LS_Neurons.combined, features = c("Scn5a", "Sst"), blend = TRUE, blend.threshold = 0)
FeaturePlot(LS_Neurons.combined, features = c("Scn5a", "Crhr2"), blend = TRUE, blend.threshold = 0)
FeaturePlot(LS_Neurons.combined, features = c("Scn5a", "Drd3"), blend = TRUE, blend.threshold = 0)
FeaturePlot(LS_Neurons.combined, features = c("Scn5a", "Glp1r"), blend = TRUE, blend.threshold = 0)
FeaturePlot(LS_Neurons.combined, features = c("Scn5a", "Nts"), blend = TRUE, blend.threshold = 0)
FeaturePlot(LS_Neurons.combined, features = c("Scn5a", "Crhr1"), blend = TRUE, blend.threshold = 0)
FeaturePlot(LS_Neurons.combined, features = c("Scn5a", "Adora2a"), blend = TRUE, blend.threshold = 0)
FeaturePlot(LS_Neurons.combined, features = c("Scn5a", "Oxtr"), blend = TRUE, blend.threshold = 0)
FeaturePlot(LS_Neurons.combined, features = c("Scn5a", "Avpr1a"), blend = TRUE, blend.threshold = 0)
FeaturePlot(LS_Neurons.combined, features = c("Scn5a", "Ntrk2"), blend = TRUE, blend.threshold = 0)

LS_Neurons.combined@meta.data$Scn5a <- LS_Neurons.combined@assays$RNA@counts["Scn5a",]
LS_Neurons.combined@meta.data$Sst <- LS_Neurons.combined@assays$RNA@counts["Sst",]
LS_Neurons.combined@meta.data$Crhr2 <- LS_Neurons.combined@assays$RNA@counts["Crhr2",]
LS_Neurons.combined@meta.data$Drd3 <- LS_Neurons.combined@assays$RNA@counts["Drd3",]
LS_Neurons.combined@meta.data$Glp1r <- LS_Neurons.combined@assays$RNA@counts["Glp1r",]
LS_Neurons.combined@meta.data$Nts <- LS_Neurons.combined@assays$RNA@counts["Nts",]
LS_Neurons.combined@meta.data$Crhr1 <- LS_Neurons.combined@assays$RNA@counts["Crhr1",]
LS_Neurons.combined@meta.data$Adora2a <- LS_Neurons.combined@assays$RNA@counts["Adora2a",]
LS_Neurons.combined@meta.data$Oxtr <- LS_Neurons.combined@assays$RNA@counts["Oxtr",]
LS_Neurons.combined@meta.data$Avpr1a <- LS_Neurons.combined@assays$RNA@counts["Avpr1a",]
LS_Neurons.combined@meta.data$Ntrk2 <- LS_Neurons.combined@assays$RNA@counts["Ntrk2",]

write.table(LS_Neurons.combined@meta.data, file = "meta.data.txt", col.names = TRUE, sep = "\t", quote = FALSE)


# hdWGCNA
theme_set(theme_cowplot())

DefaultAssay(LS_Neurons.combined) <- "integrated"
All_Neurons <- SetupForWGCNA(LS_Neurons.combined, gene_select = "variable", wgcna_name = "LS_Neurons")
length(All_Neurons@misc$LS_Neurons$wgcna_genes)

DefaultAssay(All_Neurons) <- "RNA"
All_Neurons <- MetacellsByGroups(seurat_obj = All_Neurons, group.by = "Neu_subtype", reduction = "harmony", k = 20, max_shared = 5, ident.group = "Neu_subtype")
All_Neurons <- NormalizeMetacells(All_Neurons)

All_Neurons <- SetDatExpr(All_Neurons, group_name = "N9", group.by = "Neu_subtype", assay = "RNA", slot = "data")
All_Neurons <- TestSoftPowers(All_Neurons, networkType = "signed")
plot_list <- PlotSoftPowers(All_Neurons)
wrap_plots(plot_list, ncol = 2)

power_table <- GetPowerTable(All_Neurons)
head(power_table)

## construct co-expression network
All_Neurons <- ConstructNetwork(All_Neurons, soft_power = 6, setDatExpr = FALSE, tom_name = "N9")
PlotDendrogram(All_Neurons, main = "N9 hdWGCNA Dendrogram")

All_Neurons@misc$LS_Neurons$wgcna_modules %>% head
write.table(All_Neurons@misc$LS_Neurons$wgcna_modules, file = "wgcna_modules.txt", col.names = TRUE, sep = "\t", quote = FALSE)

table(All_Neurons@misc$LS_Neurons$wgcna_modules$module)
## blue      grey turquoise     brown       red   magenta     green    yellow    purple      pink     black 
## 158      1068       567       149       119        52       119       122        51        60        73

Tom <- GetTOM(All_Neurons)

## Compute harmonized module eigengenes
All_Neurons <- ScaleData(All_Neurons, features = VariableFeatures(All_Neurons))
All_Neurons <- ModuleEigengenes(All_Neurons, group.by.vars = "sample_ID", reduction.use = "pca")

## harmonized module eigengenes
hMEs <- GetMEs(All_Neurons)
head(hMEs)
MEs <- GetMEs(All_Neurons, harmonized = FALSE)
head(MEs)

## Compute module connectivity
All_Neurons <- ModuleConnectivity(All_Neurons, group.by = "Neu_subtype", group_name = "N9")
All_Neurons <- ResetModuleNames(All_Neurons, new_name = "N9-M")

PlotKMEs(All_Neurons, ncol = 5)

modules <- GetModules(All_Neurons)

hub_df <- GetHubGenes(All_Neurons, n_hubs = 10)
head(hub_df)
write.table(hub_df, file = "hub_genes.txt", col.names = TRUE, sep = "\t", quote = FALSE)

All_Neurons <- ModuleExprScore(All_Neurons, n_genes = 25, method = "Seurat")

## Visualization
plot_list <- ModuleFeaturePlot(All_Neurons, features = "hMEs", order = TRUE)
wrap_plots(plot_list, ncol = 4) 

plot_list <- ModuleFeaturePlot(All_Neurons, features = "scores", order = "shuffle", ucell = TRUE)
wrap_plots(plot_list, ncol = 4)

MEs <- GetMEs(All_Neurons, harmonized = TRUE)
mods <- colnames(MEs); mods <- mods[mods != "grey"]

All_Neurons@meta.data <- cbind(All_Neurons@meta.data, MEs)
P <- DotPlot(All_Neurons, features = mods, group.by = "Neu_subtype")
P <- P + coord_flip() + RotatedAxis() + scale_color_gradient2(high = "red", mid = "grey95", low = "blue")
P

theme_set(theme_cowplot())
ModuleNetworkPlot(All_Neurons)

HubGeneNetworkPlot(All_Neurons, n_hubs = 5, n_other = 5, edge_prop = 0.75, mods = 'all', edge.alpha = 0.6)

## Enrichment analysis
set.seed(12345)
dbs <- c('GO_Biological_Process_2023','KEGG_2019_Mouse','Jensen_DISEASES')

All_Neurons <- RunEnrichr(All_Neurons, dbs = dbs, max_genes = Inf)
enrich_df <- GetEnrichrTable(All_Neurons)
write.table(enrich_df, file = "enrich_df.txt", col.names = TRUE, sep = "\t", quote = FALSE)

N9DO_enrichment <- read.csv(file = "N9 DO enrichment.csv")
N9DO_enrichment$Term <- factor(x = N9DO_enrichment$Term, levels = c("Mental depression", "Alcohol dependence", "Schizophrenia",
                                                                    "Cognitive disorder", "Nicotine dependence", "Bipolar disorder",
                                                                    "Major depressive disorder", "Opiate dependence", "Cannabis dependence",
                                                                    "Generalized anxiety disorder", "Heroin dependence", "Anxiety disorder",
                                                                    "Cocaine dependence", "Disease of mental health"))

ggplot() + geom_bar(data = N9DO_enrichment, aes(x = Term, y = Combined.Score, fill = Term),
                    position = position_dodge2(padding = 0.3), stat = "identity") + 
  theme(axis.title = element_blank(), axis.text.x = element_text(size = 12), axis.ticks.x = element_blank(),
        axis.ticks.length.y = unit(.25, "cm"), axis.line = element_line(color = "black", size = 0.4), 
        axis.text.y = element_text(size = 14), panel.background = element_rect(fill = "white"), legend.position = "top", 
        legend.title = element_blank(), legend.key.height = unit(.3, "cm"), legend.key.width = unit(.8, "cm"), 
        legend.direction = "vertical", legend.spacing.x = unit(.4, "cm"), legend.text = element_text(size = 12)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(-0.5,2, by = 0.5), limits = c(-0.5,2)) +  NoLegend() + RotatedAxis()




### Shuxuan Lyu 2024/03/22


