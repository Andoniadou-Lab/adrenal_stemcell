library(Seurat)
library(tidyverse)

#whole dataset
setwd("C:/Users/Alice/OneDrive - King's College London/Whole medulla P15 seq/WholeNew20210305")

P15whole.raw.data <- read.csv("AliceWholeCountTableRaw.csv", header = TRUE, row.names = 1)
P15whole <- CreateSeuratObject(counts = P15whole.raw.data, project = "P15whole", min.cells = 5, min.features = 200)
P15whole[["percent.mt"]] <- PercentageFeatureSet(object = P15whole, pattern = "^mt-")
FeatureScatter(P15whole, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(P15whole, feature1 = "nCount_RNA", feature2 = "percent.mt")
VlnPlot(object = P15whole, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
P15whole <- subset(P15whole, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 20)
FeatureScatter(P15whole, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(P15whole, feature1 = "nCount_RNA", feature2 = "percent.mt")
VlnPlot(object = P15whole, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
saveRDS(P15whole, file = "P15medullaNoProcess.rds")

P15whole1 <- P15whole
P15whole1 <- NormalizeData(object = P15whole1, normalization.method = "LogNormalize", scale.factor = 10000)
P15whole1 <- FindVariableFeatures(object = P15whole1, selection.method = "vst", nfeatures = 2000)
P15whole1top10 <- head(x = VariableFeatures(object = P15whole1), 10)
plot1 <- VariableFeaturePlot (P15whole1)
plot2 <- LabelPoints (plot = plot1, points = P15whole1top10, repel = TRUE)
plot1 + plot2

#eliminating red blood cells
VlnPlot (P15whole1, "Hbb-bs")
VlnPlot (P15whole1, "Hbb-bt")
P15wholeNBcells <- WhichCells(P15whole1, expression = `Hbb-bt` < 3)
P15wholeNB <- subset (P15whole1, cells = P15wholeNBcells)
P15wholeNBcells1 <- WhichCells (P15wholeNB, expression = `Hbb-bs` < 3)
P15wholeNB1 <- subset (P15wholeNB, cells = P15wholeNBcells1)
VlnPlot (P15wholeNB1, "Hbb-bs")
VlnPlot (P15wholeNB1, "Hbb-bt")

#continuing pre-processing
P15whole2 <-P15wholeNB1
P15whole2 <- FindVariableFeatures(object = P15whole2, selection.method = "vst", nfeatures = 2000)
P15whole2top10 <- head(x = VariableFeatures(object = P15whole2), 10)
plot1 <- VariableFeaturePlot (P15whole2)
plot2 <- LabelPoints (plot = plot1, points = P15whole2top10, repel = TRUE)
plot1 + plot2
P15whole2.genes <- rownames(P15whole2)
P15whole2 <- ScaleData(object = P15whole2, features = P15whole2.genes)
P15whole2 <- RunPCA(object = P15whole2, features = VariableFeatures(object = P15whole2))
DimPlot(object = P15whole2, reduction = "pca")
P15whole2 <- JackStraw(object = P15whole2, num.replicate = 100)
P15whole2 <- ScoreJackStraw(object = P15whole2, dims = 1:20)
JackStrawPlot(object = P15whole2, dims = 1:10)
ElbowPlot(object = P15whole2)
P15whole2 <- FindNeighbors(object = P15whole2, dims = 1:15)
VizDimLoadings(P15whole2, dims = 1:2, reduction = "pca")
saveRDS(P15whole2, file = "P15medullaWORes.rds")

#tried res0.2 to 0.8, picked 0.4 based on annotation
P15whole_Res04 <- FindClusters(object = P15whole2, resolution = 0.4)
P15whole_Res04 <- RunUMAP(object = P15whole_Res04, dims = 1:10)
DimPlot(object = P15whole_Res04, reduction = "umap", pt.size = 0.5, order = T, label = T, label.size = 3) + NoLegend()
ggsave("P15whole_Res04_SCInt_UMAP_label.tiff", plot = last_plot(), device = NULL, width = 5, height = 5)
ggsave("P15whole_Res04_SCInt_UMAP_label_medium.tiff", plot = last_plot(), device = NULL, width = 3.5, height = 3.5)
P15Whole_Res04.markers <- FindAllMarkers(object = P15whole_Res04, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
top10_Res04 <- P15Whole_Res04.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10_Res04 <-read.csv (file = "P15medulla_Res04Top10DEGs.csv")
DoHeatmap(object = P15whole_Res04, features = top10_Res04$gene, hjust = T, angle = 0, size = 4) + NoLegend ()
ggsave("P15whole_Res04_heatmap_new_small.tiff", plot = last_plot(), device = NULL, width = 4, height = 3)
ggsave("P15whole_Res04_heatmap_new.tiff", plot = last_plot(), device = NULL, width = 12, height = 9)
saveRDS(P15whole_Res04, file = "P15medulla_Res04.rds")
write.csv(top10_Res04, file = "P15medulla_Res04Top10DEGs.csv")
top100_Res04 <- P15Whole_Res04.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
write.csv(top100_Res04, file = "P15medulla_Res04Top100DEGs.csv")
top50_Res04 <- P15Whole_Res04.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(top50_Res04, file = "P15medulla_Res04Top50DEGs.csv")

#markers used for annotation
VlnPlot(P15whole_Res04, features = c(
  "Cxcl2", "C1qc", "Emcn", "Cd34","Star", "Nr5a1","Chga", "Th", "Pnmt","Sox10", "S100b", "Gfap", "Mki67", "Top2a"), pt.size = 0, same.y.lims = T, stack = T, flip = T) + NoLegend ()  
ggsave("P15whole_Res04_VlnPlot_Markers_stacked_2.tiff", plot = last_plot(), device = NULL, width = 5, height = 6)

#subsetting for medulla only clusters
P15MedSubNB1.raw <- as.matrix(GetAssayData(P15whole_Res04, slot = "counts")[, WhichCells(P15whole_Res04, ident = c('2', '5', '7', '8', '11', '12', '14'), slot = "counts")])
write.csv(P15MedSubNB1.raw, file = "P15MedSubNB1raw.csv")

#process starting from count table
P15MedSub1.raw.data <- read.csv("P15MedSubNB1raw.csv", header = TRUE, row.names = 1)
P15MedSub1 <- CreateSeuratObject(counts = P15MedSub1.raw.data, project = "P15MedSub1", min.cells = 5, min.features = 200)
P15MedSub1[["percent.mt"]] <- PercentageFeatureSet(object = P15MedSub1, pattern = "^mt-")
VlnPlot(object = P15MedSub1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#subsetted enough
P15MedSub1 <- subset(P15MedSub1, subset = nFeature_RNA > 500 & nFeature_RNA < 5500 & percent.mt < 20)
#did not run the above
FeatureScatter(P15MedSub1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
saveRDS(P15MedSub1, file = "P15MedSub1NoProcess.rds")

P15MedSub1a <- P15MedSub1
P15MedSub1a <- NormalizeData(object = P15MedSub1a, normalization.method = "LogNormalize", scale.factor = 10000)
P15MedSub1a <- FindVariableFeatures(object = P15MedSub1a, selection.method = "vst", nfeatures = 2000)
P15MedSub1atop10 <- head(x = VariableFeatures(object = P15MedSub1a), 10)
plot1 <- VariableFeaturePlot (P15MedSub1a)
plot2 <- LabelPoints (plot = plot1, points = P15MedSub1atop10, repel = TRUE)
plot1 + plot2
P15MedSub1a.genes <- rownames(P15MedSub1a)
P15MedSub1a <- ScaleData(object = P15MedSub1a, features = P15MedSub1a.genes)
P15MedSub1a <- RunPCA(object = P15MedSub1a, features = VariableFeatures(object = P15MedSub1a))
DimPlot(object = P15MedSub1a, reduction = "pca")
P15MedSub1a <- JackStraw(object = P15MedSub1a, num.replicate = 100)
P15MedSub1a <- ScoreJackStraw(object = P15MedSub1a, dims = 1:20)
JackStrawPlot(object = P15MedSub1a, dims = 1:20)
ElbowPlot(object = P15MedSub1a)
P15MedSub1a <- FindNeighbors(object = P15MedSub1a, dims = 1:8)
saveRDS(P15MedSub1a, file = "P15MedSub1WORes.rds")
VizDimLoadings(P15MedSub1a, dims = 1:2, reduction = "pca")

#tried res0.2 to 0.8, picked 0.4 based on annotation
P15MedSub1_Res04 <- FindClusters(object = P15MedSub1a, resolution = 0.4)
P15MedSub1_Res04 <- RunUMAP(object = P15MedSub1_Res04, dims = 1:8)
DimPlot(object = P15MedSub1_Res04, reduction = "umap")
ggsave("P15MedSub1_Res04_SCInt_UMAP.tiff", plot = last_plot(), device = NULL)
P15MedSub1_Res04.markers <- FindAllMarkers(object = P15MedSub1_Res04, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
top10_Res04 <- P15MedSub1_Res04.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(object = P15MedSub1_Res04, features = top10_Res04$gene) + NoLegend()
ggsave("P15MedSub1_Res04_heatmap.tiff", plot = last_plot(), device = NULL)

saveRDS(P15MedSub1_Res04, file = "P15MedSub1_Res04.rds")
write.csv(top10_Res04, file = "P15MedSub1_Res04Top10DEGs.csv")

#markers used for annotation
VlnPlot(P15MedSub1_Res04, features = c(
  "Cxcl2", "C1qc", "Emcn", "Cd34","Star", "Nr5a1","Chga", "Th", "Pnmt","Sox10", "S100b", "Gfap", "Mki67", "Top2a"), pt.size = 0, same.y.lims = T)  
#needs further subsetting

#subsetting for medulla only clusters
P15MedSubNB2.raw <- as.matrix(GetAssayData(P15MedSub1_Res04, slot = "counts")[, WhichCells(P15MedSub1_Res04, ident = c('0', '1', '2', '3', '4', '5', '9'), slot = "counts")])
write.csv(P15MedSubNB2.raw, file = "P15MedSub2raw.csv")

#final "medulla only" dataset

P15MedSub2.raw.data <- read.csv("P15MedSub2raw.csv", header = TRUE, row.names = 1)

P15MedSub2 <- CreateSeuratObject(counts = P15MedSub2.raw.data, project = "P15MedSub2", min.cells = 5, min.features = 200)
P15MedSub2[["percent.mt"]] <- PercentageFeatureSet(object = P15MedSub2, pattern = "^mt-")
VlnPlot(object = P15MedSub2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#subsetted enough
P15MedSub2 <- subset(P15MedSub2, subset = nFeature_RNA > 500 & nFeature_RNA < 5500 & percent.mt < 20)
#did not run the above
FeatureScatter(P15MedSub2, feature1 = "nCount_RNA", feature2 = "percent.mt")
saveRDS(P15MedSub2, file = "P15MedSub2NoProcess.rds")

P15MedSub2a <- P15MedSub2
P15MedSub2a <- NormalizeData(object = P15MedSub2a, normalization.method = "LogNormalize", scale.factor = 10000)
P15MedSub2a <- FindVariableFeatures(object = P15MedSub2a, selection.method = "vst", nfeatures = 2000)
P15MedSub2atop10 <- head(x = VariableFeatures(object = P15MedSub2a), 10)
plot1 <- VariableFeaturePlot (P15MedSub2a)
plot2 <- LabelPoints (plot = plot1, points = P15MedSub2atop10, repel = TRUE)
plot1 + plot2
P15MedSub2a.genes <- rownames(P15MedSub2a)
P15MedSub2a <- ScaleData(object = P15MedSub2a, features = P15MedSub2a.genes)
P15MedSub2a <- RunPCA(object = P15MedSub2a, features = VariableFeatures(object = P15MedSub2a))
DimPlot(object = P15MedSub2a, reduction = "pca")
P15MedSub2a <- JackStraw(object = P15MedSub2a, num.replicate = 100)
P15MedSub2a <- ScoreJackStraw(object = P15MedSub2a, dims = 1:20)
JackStrawPlot(object = P15MedSub2a, dims = 1:20)
ElbowPlot(object = P15MedSub2a)
P15MedSub2a <- FindNeighbors(object = P15MedSub2a, dims = 1:8)
VizDimLoadings(P15MedSub2a, dims = 1:2, reduction = "pca")
saveRDS(P15MedSub2a, file = "P15MedSub2WORes.rds")
P15MedSub2a <- ("P15MedSub2WORes.rds")

#tried res0.2 to 0.8, picked 0.4 based on annotation
P15MedSub2_Res04 <- FindClusters(object = P15MedSub2a, resolution = 0.4)
P15MedSub2_Res04 <- RunUMAP(object = P15MedSub2_Res04, dims = 1:8)
DimPlot(object = P15MedSub2_Res04, reduction = "umap", pt.size = 0.5) + NoLegend ()
ggsave("P15MedSub2_Res04_UMAP_small.tiff", plot = last_plot(), device = NULL, width = 3.5, height = 3)
P15MedSub2_Res04.markers <- FindAllMarkers(object = P15MedSub2_Res04, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
top10_Res04 <- P15MedSub2_Res04.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(object = P15MedSub2_Res04, features = top10_Res04$gene) + NoLegend()
ggsave("P15MedSub2_Res04_heatmap_legend.tiff", plot = last_plot(), device = NULL, width = 8, height = 10)
P15MedSub2_Res04.markers <- FindAllMarkers(object = P15MedSub2_Res04, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
top5_Res04 <- P15MedSub2_Res04.markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_log2FC)
DoHeatmap(object = P15MedSub2_Res04, features = top5_Res04$gene) + NoLegend()
ggsave("P15MedSub2_Res04_heatmap_top5.tiff", plot = last_plot(), device = NULL, width = 5, height = 4)

saveRDS(P15MedSub2_Res04, file = "P15MedSub2_Res04.rds")
write.csv(top10_Res04, file = "P15MedSub2_Res04Top10DEGs.csv")
top50_Res04 <- P15MedSub2_Res04.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.csv(top50_Res04, file = "P15MedSub2_Res04Top50DEGs.csv")
top100_Res04 <- P15MedSub2_Res04.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
write.csv(top100_Res04, file = "P15MedSub2_Res04Top100DEGs.csv")

#markers for "whole adrenal" annotation
p <- DotPlot(P15MedSub2_Res04, features = c(
  "Cxcl2", "C1qc", "Emcn", "Cd34","Star", "Nr5a1","Chga", "Th", "Pnmt","Sox10", "S100b", "Gfap", "Mki67", "Top2a"), cols = "RdYlBu", dot.scale = 12, scale.by = "size")
p + theme(axis.text.x = element_text(angle = 90, face ="italic")) 
ggsave("P15MedSub2_Res04_DotPlotMarkers.tiff", plot = last_plot(), device = NULL)

#specific markers for adrenomedullary cells annotation
VlnPlot(P15MedSub2_Res04, features = c(
  "Chga", "Th", "Ddc", "Dbh", "Pnmt","Sox10", "S100b", "Gfap", "Mki67", "Top2a"), pt.size = 0, same.y.lims = T)  
ggsave("P15MedSub2_Res04_VlnPlotMarkersMedulla.tiff", plot = last_plot(), device = NULL, width = 7, height = 5)

FeaturePlot(P15MedSub2_Res04, features = c(
  "Chga", "Th", "Ddc", "Dbh", "Pnmt","Sox10", "S100b", "Gfap", "Mki67", "Top2a"), pt.size = 0.5, order = T, min.cutoff = 0, ncol = 5, cols = c("lightgrey", "red"), keep.scale = "all") 
ggsave("P15MedSub2_Res04_FeaturPlotMarkersMedulla.tiff", plot = last_plot(), device = NULL, width = 12.5, height = 4.5)

#naming clusters based on adrenomedulalry annotation
new.cluster.ids <- c(
  "0-AdrChr", 
  "1-AdrChr", 
  "2-NorAdrChr", 
  "3-T-NorAdrChr", 
  "4-Sustentacular", 
  "5-AdrChr", 
  "6-T-AdrChr",
  "7-Cycling AdrChr"
  
)
names(new.cluster.ids) <- levels(P15MedSub2_Res04)
P15MedSub2_Res04_Named <- RenameIdents(P15MedSub2_Res04, new.cluster.ids) 
saveRDS(P15MedSub2_Res04_Named, file = "P15MedSub2_Res04_Named.rds")

#ordering clusters
myLevels <- c(
  "4-Sustentacular",
  "3-T-NorAdrChr",
  "2-NorAdrChr",
  "6-T-AdrChr",
  "0-AdrChr", 
  "1-AdrChr", 
  "5-AdrChr", 
  "7-Cycling AdrChr") 
factor(Idents(P15MedSub2_Res04_Named), levels= myLevels) 
Idents(P15MedSub2_Res04_Named) <- factor(Idents(P15MedSub2_Res04_Named), levels= myLevels)
P15MedSub2_Res04_Named_ordered <- P15MedSub2_Res04_Named
saveRDS(P15MedSub2_Res04_Named_ordered, file = "P15MedSub2_Res04_Named_ordered.rds")

#"canonical" sustentacular cell markers
FeaturePlot (P15MedSub2_Res04, features = c("S100b", "Gfap", "Sox10", "Vim", "Nes"), order = T, pt.size = 0.5,  cols = c("lightgrey", "red"), keep.scale = "all", min.cutoff = 0, ncol = 5)
ggsave("P15MedSub2_Res04_FeaturPlotKnownSustMarkers.tiff", plot = last_plot(), device = NULL, width = 12.5, height = 2.25)

#new sustentacular cell markers
FeaturePlot (P15MedSub2_Res04, features = c("Scn7a", "Plp1", "Col28a1", "Kcna1", "Lgi4", "Cdh19", "Serpine2", "Atp1a2", "Apod", "Lpar1", "Sostdc1", "Fabp7", "Gfra3", "Scrg1", "Col18a1", "Sfrp1", "Postn", "Kcna2", "Htra1"), order = T, pt.size = 0.5,  cols = c("lightgrey", "red"), keep.scale = "all", min.cutoff = 0, ncol = 4)
ggsave("P15MedSub2_Res04_FeaturPlotNewSustMarkers.tiff", plot = last_plot(), device = NULL, width = 10, height = 11.25)

#old + new sustentacular cell markers
FeaturePlot(P15MedSub2_Res04, features = c("S100b", "Gfap", "Sox10", "Plp1", "Sox2"), order = T, pt.size = .5,  cols = c("lightgrey", "red"), keep.scale = "all", min.cutoff = 0, ncol = 2)
ggsave("P15MedSub2_Res04_FeaturPlotSustSox2.tiff", plot = last_plot(), device = NULL, width = 4, height = 5.4)

#noradrenaline cell markers
FeaturePlot(P15MedSub2_Res04, features = c("Cox8b", "Lix1", "Penk" ), order = T, pt.size = .5,  cols = c("lightgrey", "red"), keep.scale = "all", min.cutoff = 0, ncol = 3)
ggsave("P15MedSub2_Res04_FeaturPlotNorMarkers.tiff", plot = last_plot(), device = NULL, width = 7.5, height = 2.25)

#showing overlap of marker genes
FeaturePlot(P15MedSub2_Res04, features = c("Sox10", "Th"), order = T, pt.size = 0.4,  cols = c("grey28","deeppink3", "greenyellow"), blend = T, blend.threshold = 0.1)
ggsave("P15MedSub2_Res04_FeaturePlotSox10Th.tiff", plot = last_plot(), device = NULL, width = 9, height = 3)
FeaturePlot(P15MedSub2_Res04, features = c("Sox10", "Penk"), order = T, pt.size = 0.4,  cols = c("grey28","deeppink3", "greenyellow"), blend = T, blend.threshold = 0.1)
ggsave("P15MedSub2_Res04_FeaturePlotSox10Penk.tiff", plot = last_plot(), device = NULL, width = 9, height = 3)
FeaturePlot(P15MedSub2_Res04, features = c("Sox10", "Pnmt"), order = T, pt.size = 0.4,  cols = c("grey28","deeppink3", "greenyellow"), blend = T, blend.threshold = 0.1)
ggsave("P15MedSub2_Res04_FeaturePlotSox10Pnmt.tiff", plot = last_plot(), device = NULL, width = 9, height = 3)
FeaturePlot(P15MedSub2_Res04, features = c("Sox2", "Th"), order = T, pt.size = 0.4,  cols = c("grey28","deeppink3", "greenyellow"), blend = T, blend.threshold = 0.1)
ggsave("P15MedSub2_Res04_FeaturePlotSox2Th.tiff", plot = last_plot(), device = NULL, width = 9, height = 3)
FeaturePlot(P15MedSub2_Res04, features = c("Sox2", "Penk"), order = T, pt.size = 0.4,  cols = c("grey28","deeppink3", "greenyellow"), blend = T, blend.threshold = 0.1)
ggsave("P15MedSub2_Res04_FeaturePlotSox2Penk.tiff", plot = last_plot(), device = NULL, width = 9, height = 3)
FeaturePlot(P15MedSub2_Res04, features = c("Sox2", "Pnmt"), order = T, pt.size = 0.4,  cols = c("grey28","deeppink3", "greenyellow"), blend = T, blend.threshold = 0.1)
ggsave("P15MedSub2_Res04_FeaturePlotSox2Pnmt.tiff", plot = last_plot(), device = NULL, width = 9, height = 3)

#cell cycle analysis
setwd("C:/Users/Alice/OneDrive - King's College London/Whole medulla P15 seq/Cell cycle")

s.genes <- read.table(file = "Sgenelist.txt", header = TRUE)
s.genes <- as.character(s.genes[,1])
g2m.genes <- read.table(file = "G2Mgenelist.txt", header = TRUE)
g2m.genes <- as.character(g2m.genes[,1])

setwd("C:/Users/Alice/OneDrive - King's College London/Whole medulla P15 seq/WholeNew20210305/MedSub1/MedSub2")

P15MedSub2_Res04_CC <- CellCycleScoring(P15MedSub2_Res04, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(P15MedSub2_Res04_CC[[]])

P15MedSub2_Res04_CC <- RunPCA(P15MedSub2_Res04_CC, features = c(s.genes, g2m.genes))
saveRDS(P15MedSub2_Res04_CC, file = "P15MedSub2_Res04_CC.rds")

DimPlot(P15MedSub2_Res04_CC, pt.size = 1, order = T, split.by = "Phase")

S <- WhichCells(P15MedSub2_Res04_CC, idents = "S")
G1 <- WhichCells(P15MedSub2_Res04_CC, idents = "G1")
G2M <- WhichCells(P15MedSub2_Res04_CC, idents = "G2M")

DimPlot(P15MedSub2_Res04_CC, pt.size = 1, order = T, cols= c("deeppink3", "deepskyblue4", "darkslategray3"))
ggsave("P15MedSub2_Res04_CC_1.tiff", plot = last_plot(), device = NULL, width = 5, height = 4.5)
DimPlot(P15MedSub2_Res04_CC, pt.size = 2, order = T)
ggsave("P15MedSub2only_Res04_CC_2.tiff", plot = last_plot(), device = NULL)
