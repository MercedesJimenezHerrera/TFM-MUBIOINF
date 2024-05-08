# Cargar el archivo .RDS
library(dplyr)
library(org.Mm.eg.db)
bonemarrow <- readRDS(file = "C:/Users/merce/Documents/Máster en Bioinformática/TFM/be8db92e-e52c-430c-b3d1-c27b506f8c2e.rds")
str(bonemarrow)
slotNames(bonemarrow)
View(bonemarrow@meta.data)
# Extraer la matriz de conteo de expresión del objeto Seurat 
### En este código, pbmc.data1@assays$RNA@counts accede a la matriz de conteo de expresión dentro del objeto Seurat pbmc.data1. Luego, esta matriz de conteo se pasa como argumento counts a CreateSeuratObject
counts_matrix <- bonemarrow@assays$RNA@counts
###Intento de agregar los genesymbol###
#ensembl_ids <- rownames(counts_matrix)
#gene_symbols <- select(org.Mm.eg.db, keys=ensembl_ids, columns="SYMBOL", keytype="ENSEMBL")

# Ensure the number of gene symbols in count matrix is same as list to add
#nrow(counts_matrix) == length(gene_symbols$SYMBOL)

# If the above is FALSE, you will need to modify gene_symbols$SYMBOL to correct
# Something using match might work.
#gene_symbols <- gene_symbols[match(rownames(counts_matrix), gene_symbols$ENSEMBL), ]
#rownames(counts_matrix) <- gene_symbols$SYMBOL
# Identificar y eliminar filas duplicadas
duplicated_rows <- apply(counts_matrix, 1, anyDuplicated)
counts_matrix_unique <- counts_matrix[!duplicated_rows == 0, ]
# Agregar nombres de fila únicos a la matriz de recuento
#rownames(counts_matrix_unique) <- make.unique(rownames(counts_matrix_unique))

# Crear objeto Seurat con la matriz de recuento sin filas duplicadas
#bonecounts <- CreateSeuratObject(counts_matrix_unique, project = "seuratM", 
                                 assay = "RNA", min.cells = 0, min.features = 0, 
                                 names.field = 1, names.delim = "_", 
                                 meta.data = bonemarrow@meta.data)
View(bonecounts@meta.data)

head(counts_matrix)
View(counts_matrix)
# Crear un nuevo objeto Seurat utilizando la matriz de conteo extraída
bonecounts <- CreateSeuratObject(counts = counts_matrix, project = "seuratM")
bonecounts

# Control de calidad (QCONTROL)

View(bonecounts@meta.data)
#1. PERCENT-MT
#calcular el porcentaje de mt
# Calcular el porcentaje de expresión génica mitocondrial y asignarlo al objeto Seurat 

#bonecounts[["percent.mt"]] <- PercentageFeatureSet(bonecounts, pattern = "^Mt-") #Paso que no ha podido realizarse ya que no se han incluido los genes mitocondriales
#View(bonecounts@meta.data)

# Visualize QC metrics as a violin plot
VlnPlot(bonecounts, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
library(ggplot2)
FeatureScatter(bonecounts, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method ='lm')
range(bonecounts$nFeature_RNA)
#bonecounts<- subset(bonecounts, subset = nFeature_RNA > 200 & nFeature_RNA < 2499)
bonecounts
#FeatureScatter(bonecounts, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

##Nos faltaria otro comparando o counts o features con percent.mt
#plot1<- FeatureScatter(bonecounts, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot1+plot2
#the example below, we visualize QC metrics, and use these to filter cells.

#We filter cells that have unique feature counts over 2,500 or less than 200
#We filter cells that have >5% mitochondrial counts


bonecounts <- subset(bonecounts, subset = nFeature_RNA > 200 & nFeature_RNA < 2000)

#normalize data
bonecounts<- NormalizeData(bonecounts, normalization.method = "LogNormalize", scale.factor = 10000)
str(bonecounts)

##Identificación de características altamente variables (selección de características)
bonecounts <- FindVariableFeatures(bonecountsn, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(bonecounts), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(bonecounts)
LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1

all.genes <- rownames(bonecounts)
bonecounts<- ScaleData(bonecounts, features = all.genes)

#linear dimensionality reduction
str(bonecounts)
#Para los primeros componentes principales, Seurat genera una lista de genes con las cargas más positivas y negativas, que representan módulos de genes que exhiben correlación (o anticorrelación) entre células individuales en el conjunto de datos.
bonecounts<- RunPCA(bonecounts, features = VariableFeatures(object = bonecounts))
# Examine and visualize PCA results a few different ways
print(bonecounts[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(bonecounts, dims = 1:2, reduction = "pca")
DimPlot(bonecounts, reduction = "pca") + NoLegend()
DimHeatmap(bonecounts, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(bonecounts, dims = 1:15, cells = 500, balanced = TRUE)
#Para superar el gran ruido técnico en cualquier característica única para los datos de scRNA-seq, Seurat agrupa las células en función de sus puntuaciones de PCA, y cada PC representa esencialmente una "metacaracterística" que combina información en un conjunto de características correlacionadas. Por lo tanto, los componentes principales superiores representan una compresión sólida del conjunto de datos
ElbowPlot(bonecounts)
#Agrupar las células (clustering)
bonecounts<- FindNeighbors(bonecounts, dims = 1:15)
bonecounts<- FindClusters(bonecounts, resolution = c(0,0.1,0.3,0.5,0.7,1))
View(bonecounts@meta.data)
DimPlot(bonecounts, group.by = "RNA_snn_res.0.1", label = TRUE)
#Setting identity of clusters
#Look at cluster IDs of the first 5 cells
Idents(bonecounts)
Idents(bonecounts)<-"RNA_snn_res.0.1"
Idents(bonecounts)
#Ejecutar reducción dimensional no lineal (UMAP/tSNE)
bonecounts <- RunUMAP(bonecounts, dims = 1:15)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(bonecounts, reduction = "umap")

#Encontrar características expresadas diferencialmente (biomarcadores de grupo)
# find all markers of cluster 2

# Check available cluster identities
Idents(pbmc)

# Adjust ident.1 parameter based on available identities
FindAllMarkers(bonecounts, logfc.threshold = 0.25, min.pct = 0.1, only.pos = TRUE, slot = 'counts')
cluster2.markers <- FindMarkers(bonecounts, ident.1 = "2")

#cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n = 5)
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(bonecounts, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
every.markers <- FindAllMarkers(bonecounts, only.pos = TRUE)
every.markers %>%
group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
cluster0.markers <- FindMarkers(bonecounts, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(bonecounts, features = c("ENSMUSG00000024659", "ENSMUSG00000056071"))
# you can plot raw counts as well
VlnPlot(bonecounts, features = c("ENSMUSG00000024659", "ENSMUSG00000056071"), layer = "counts", log = TRUE)
FeaturePlot(bonecounts, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP","CD8A"))
every.markers %>%
group_by(cluster) %>%
dplyr::filter(avg_log2FC > ) %>%
slice_head(n = 10) %>%
ungroup() -> top10
DoHeatmap(bonecounts, features = top10$gene) + NoLegend()

#Asignación de identidad de tipo celular a grupos
#new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
#"NK", "DC", "Platelet")
#names(new.cluster.ids) <- levels(pbmc)
#pbmc <- RenameIdents(pbmc, new.cluster.ids)
#DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
#library(ggplot2)
#plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
#theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
#ggsave(filename = "../output/images/pbmc3k_umap.jpg", height = 7, width = 12, plot = plot, quality = 50)
