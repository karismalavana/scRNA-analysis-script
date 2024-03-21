#script that performs differential gene analysis on scRNA-seq data
#data: 	Single cell RNA-seq of primary human glioblastomas
#data source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57872

#load libraries
library(Seurat)
library(tidyverse)
library(stringr)

#load dataset
raw.data <- read.table('/Users/karismalavana/Desktop/PolygenceProject/GSE57872_GBM_data_matrix.txt')
raw.data[1:5,1:5]

#create seurat obj with gene expression counts matrix
seurat.obj <- CreateSeuratObject(counts = raw.data, project = "Glioblastomas", min.cells = 3, min.features = 200)
s = apply(raw.data, 1, min)
normdata = raw.data - s
normdata[1:5,1:5]
mat = as.matrix(normdata)
SetAssayData(seurat.obj, assay ='RNA', slot ='scale.data', mat)


#1. Quality Control  ---------------
# %MT, poor quality cells have higher % mitochondrial genes bc it indicates cell dying
#poor quality cells has extremely low # of genes or high # of genes (due to doublets or multiple cells sequenced together)
seurat.obj[['percent.mt']]<- PercentageFeatureSet(seurat.obj, pattern = "^MT-")
View(seurat.obj@meta.data)

VlnPlot(seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol = 3)
FeatureScatter(seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')

#2.Filtering (based on meta data columns)  ---------------
seurat.obj <- subset(seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
seurat.obj #5948 features across 183 samples within 1 assay 


#3. Normalization --------------- NOT DONE BC DATA ALREADY NORMALIZED
#seurat.obj <- NormalizeData(seurat.obj)
#str(seurat.obj)

#4. Identify Highly Variable Features  ---------------

#Identify 10 most variable features
VariableFeatures(seurat.obj,assay= 'RNA') <- row.names(raw.data)
topTen <- head(VariableFeatures(seurat.obj),10)

#plot variable features
plotVar <- VariableFeaturePlot(seurat.obj)
plotVar
FindVariableFeatures(seurat.obj)
LabelPoints(plot = plotVar, points = topTen, repel = TRUE)
#Cannot find points provided

#5. Scaling (remove unwanted sources of variation)  ---------------
seurat.obj <- ScaleData(seurat.obj, features =rownames(seurat.obj))

#6. Perform Linear Dimensionality Reduction (PCA)  ---------------
seurat.obj <- RunPCA(seurat.obj, features = VariableFeatures(obj = seurat.obj))

#visualize PCA results
print(seurat.obj[['pca']], dims = 1:5, nfeatures = 5)
DimHeatmap(seurat.obj, dims = 1, cells = 500, balanced = TRUE)

#determine dimensionality
ElbowPlot(seurat.obj)

#7. Clustering ---------------
#cluster cells with similar expression patterns
seurat.obj <- FindNeighbors(seurat.obj, dims = 1:15)

#resolution
seurat.obj <- FindClusters(seurat.obj, resolution = c(0.1, 0.3, 0.5, 0.7, 1,2))

#UMAP
seurat.obj <- RunUMAP(seurat.obj, dims = 1:15)

View(seurat.obj@meta.data)

DimPlot(seurat.obj, group.by = "RNA_snn_res.0.5", label = TRUE)

#identity of clusters
Idents(seurat.obj) <- 'RNA_snn_res.0.5'
markers <- FindAllMarkers(seurat.obj)

# find markers for every cluster compared to all remaining cells, report only the positive ones
comp_markers <- FindAllMarkers(seurat.obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

FeaturePlot(seurat.obj, features = c('EGFR'), min.cutoff = 'q10')|FeaturePlot(seurat.obj, features = c('FIGNL1'), min.cutoff = 'q10')
FeaturePlot(seurat.obj, features = c('S100A6'), min.cutoff = 'q10')|FeaturePlot(seurat.obj, features = c('C3'), min.cutoff = 'q10')
FeaturePlot(seurat.obj, features = c('SPP1'), min.cutoff = 'q10')|FeaturePlot(seurat.obj, features = c('CHI3L1'), min.cutoff = 'q10')

pdf(file='geneperident.pdf', height = 8, width = 8)
FeaturePlot(seurat.obj, features = c('FIGNL1'), split.by = 'orig.ident')
dev.off()

DimPlot(seurat.obj, group.by = "orig.ident", label = TRUE)

#cluster 0 marker
FeaturePlot(seurat.obj, features = c('EGFR'), min.cutoff = 'q10')|FeaturePlot(seurat.obj, features = c('LHFPL3'), min.cutoff = 'q10', cols= c("lightgrey","red"))

#cluster 1 marker
FeaturePlot(seurat.obj, features = c('SPP1'), min.cutoff = 'q10', cols= c("lightgrey","red"))|FeaturePlot(seurat.obj, features = c('CHI3L1'), min.cutoff = 'q10')

#cluster 2 marker
FeaturePlot(seurat.obj, features = c('TF'), min.cutoff = 'q10', cols= c("lightgrey","red"))|FeaturePlot(seurat.obj, features = c('CCDC170'), min.cutoff = 'q10')

#cluster 3 marker
FeaturePlot(seurat.obj, features = c('FIGNL1'), min.cutoff = 'q10', cols= c("lightgrey","red"))|FeaturePlot(seurat.obj, features = c('ATP1A2'), min.cutoff = 'q10')

#cluster 4 marker
FeaturePlot(seurat.obj, features = c('C3'), min.cutoff = 'q10',  cols= c("lightgrey","red"))|FeaturePlot(seurat.obj, features = c('SOD2'), min.cutoff = 'q10')

#cluster 5 marker
FeaturePlot(seurat.obj, features = c('SLC38A1'), min.cutoff = 'q10', cols= c("lightgrey","red"))|FeaturePlot(seurat.obj, features = c('MYL12A'), min.cutoff = 'q10')

#cluster 6 marker
FeaturePlot(seurat.obj, features = c('UBB'), min.cutoff = 'q10', cols= c("lightgrey","red"))|FeaturePlot(seurat.obj, features = c('S100A6'), min.cutoff = 'q10')

#patient information

install.packages("openxlsx")
library(openxlsx)

write.xlsx(seurat.obj@meta.data, '/Users/karismalavana/Desktop/PolygenceProject', asTable = TRUE, overwrite = TRUE)
