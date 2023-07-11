library(Seurat)
library(tidyverse)

setwd("C:/Users/Alice/OneDrive - King's College London/SOX2scRNAseq/All")

P15Sox2.raw.data <- read.csv("AliceCountTableRaw.csv", header = TRUE, row.names = 1)
P15Sox2 <- CreateSeuratObject(counts = P15Sox2.raw.data, project = "P15Sox2", min.cells = 5, min.features = 200)
P15Sox2[["percent.mt"]] <- PercentageFeatureSet(object = P15Sox2, pattern = "^mt-")
VlnPlot(object = P15Sox2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
P15Sox2 <- subset(P15Sox2, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & percent.mt < 20)
FeatureScatter(P15Sox2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
saveRDS(P15Sox2, file = "P15Sox2NoProcess.rds")

P15Sox2a <- P15Sox2
P15Sox2a <- NormalizeData(object = P15Sox2a, normalization.method = "LogNormalize", scale.factor = 10000)
P15Sox2a <- FindVariableFeatures(object = P15Sox2a, selection.method = "vst", nfeatures = 2000)
P15Sox2atop10 <- head(x = VariableFeatures(object = P15Sox2a), 10)
P15Sox2a.genes <- rownames(P15Sox2a)
P15Sox2a <- ScaleData(object = P15Sox2a, features = P15Sox2a.genes)
P15Sox2a <- RunPCA(object = P15Sox2a, features = VariableFeatures(object = P15Sox2a))
DimPlot(object = P15Sox2a, reduction = "pca")
P15Sox2a <- JackStraw(object = P15Sox2a, num.replicate = 100)
P15Sox2a <- ScoreJackStraw(object = P15Sox2a, dims = 1:20)
JackStrawPlot(object = P15Sox2a, dims = 1:20)
ElbowPlot(object = P15Sox2a)
P15Sox2a <- FindNeighbors(object = P15Sox2a, dims = 1:8)
saveRDS(P15Sox2a, file = "P15Sox2aWORes.rds")

#tried res 0.2 to 0.8, picked 0.4
P15Sox2_Res04 <- FindClusters(object = P15Sox2a, resolution = 0.4)
P15Sox2_Res04 <- RunUMAP(object = P15Sox2_Res04, dims = 1:8)
DimPlot(object = P15Sox2_Res04, reduction = "umap")
ggsave("P15Sox2_Res04_UMAP.tiff", plot = last_plot(), device = NULL)
P15Sox2_Res04.markers <- FindAllMarkers(object = P15Sox2_Res04, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- P15Sox2_Res04.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(object = P15Sox2_Res04, features = top10$gene) + NoLegend()
ggsave("P15Sox2_Res04_Heatmap.tiff", plot = last_plot(), device = NULL)
VlnPlot(P15Sox2_Res04, features = "Sox2")
ggsave("P15Sox2_Res04_VlnSox2.tiff", plot = last_plot(), device = NULL)
FeaturePlot(P15Sox2_Res04, "Sox2")
ggsave("P15Sox2_Res04_FeatureSox2.tiff", plot = last_plot(), device = NULL)

#subsetting for Sox2 expression
raw.data.Sox2pos <- as.matrix(GetAssayData(P15Sox2_Res04, slot = "counts")[, (WhichCells(P15Sox2_Res04, expression = (Sox2 > 0)))])
write.csv(raw.data.Sox2pos, file = "Sox2only.csv")

#analysis of the subsetted matrix 
P15Sox2.raw.data <- read.csv("Sox2only.csv", header = TRUE, row.names = 1)

P15Sox2 <- CreateSeuratObject(counts = P15Sox2.raw.data, project = "P15Sox2", min.cells = 5, min.features = 200)
P15Sox2[["percent.mt"]] <- PercentageFeatureSet(object = P15Sox2, pattern = "^mt-")
VlnPlot(object = P15Sox2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

saveRDS(P15Sox2, file = "P15Sox2onlyNoProcess.rds")

P15Sox2a <- P15Sox2
P15Sox2a <- NormalizeData(object = P15Sox2a, normalization.method = "LogNormalize", scale.factor = 10000)
P15Sox2a <- FindVariableFeatures(object = P15Sox2a, selection.method = "vst", nfeatures = 2000)
P15Sox2atop10 <- head(x = VariableFeatures(object = P15Sox2a), 10)
plot1 <- VariableFeaturePlot (P15Sox2a)
plot2 <- LabelPoints (plot = plot1, points = P15Sox2atop10, repel = TRUE)
plot1 + plot2
P15Sox2a.genes <- rownames(P15Sox2a)
P15Sox2a <- ScaleData(object = P15Sox2a, features = P15Sox2a.genes)
P15Sox2a <- RunPCA(object = P15Sox2a, features = VariableFeatures(object = P15Sox2a))
DimPlot(object = P15Sox2a, reduction = "pca")
P15Sox2a <- JackStraw(object = P15Sox2a, num.replicate = 100)
P15Sox2a <- ScoreJackStraw(object = P15Sox2a, dims = 1:20)
JackStrawPlot(object = P15Sox2a, dims = 1:20)
ElbowPlot(object = P15Sox2a)
P15Sox2a <- FindNeighbors(object = P15Sox2a, dims = 1:13)
saveRDS(P15Sox2a, file = "P15Sox2aonlyWORes.rds")

#tried res 0.8 to 0.2, picked 0.4
P15Sox2_Res04 <- FindClusters(object = P15Sox2a, resolution = 0.4)
P15Sox2_Res04 <- RunUMAP(object = P15Sox2_Res04, dims = 1:13)
DimPlot(object = P15Sox2_Res04, reduction = "umap")
ggsave("P15Sox2only_Res04_UMAP.tiff", plot = last_plot(), device = NULL)
P15Sox2_Res04.markers <- FindAllMarkers(object = P15Sox2_Res04, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
top10_Res04 <- P15Sox2_Res04.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(object = P15Sox2_Res04, features = top10_Res04$gene) + NoLegend()
ggsave("P15Sox2_Res04only_heatmap.tiff", plot = last_plot(), device = NULL)
saveRDS(P15Sox2_Res04, "P15Sox2_Res04.rds")
write.csv(top10_Res04, "P15Sox2_REs04_top10DEGs.csv")
top100_Res04 <- P15Sox2_Res04.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
write.csv(top100_Res04, "P15Sox2_REs04_top100DEGs.csv")

FeaturePlot (P15Sox2_Res04, "Sox2", order = T, cols = c("lightgrey", "red"), min.cutoff = 0, pt.size = 0.5)
ggsave("P15Sox2_Res04_FeaturePlot_Sox2.tiff", plot = last_plot(), device = NULL, width = 2.5, height = 2.25)
