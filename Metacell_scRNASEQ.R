##Build Metacells with Supercell and Analyze Metacells with Seurat
install.packages("supercells", repos = "https://nowosad.r-universe.dev")
library(supercells)
library(Seurat)
library(tidyverse)
library(ggplot2)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SingleCellExperiment")

if (!requireNamespace("remotes")) install.packages("remotes")
  remotes::install_github("GfellerLab/SuperCell")
  remotes::install_github("mariiabilous/SuperCellBM")
    
library(SingleCellExperiment)
library(SuperCell)
library(SuperCellBM)

### 1. Load data
bonemarrowM <- readRDS(file = "C:/Users/merce/Documents/Máster en Bioinformática/TFM/be8db92e-e52c-430c-b3d1-c27b506f8c2e.rds")
str(bonemarrow)
slotNames(bonemarrow)
View(bonemarrow@meta.data)
DimPlot(bonemarrowM, group.by = "cell_type", label = T, raster = FALSE)

#Get gene expression matrix
GE<- bonemarrowM@assays$RNA@counts
dim(GE)

### 2. Build Metacells with Supercell
#Set the graining level gamma =10~50 and k.knn to simplify single-cell data
gamma<- 20 # graining level (nivel de granulado)
k.knn<- 5 #numero de vecinos cercanos para construir la red knn
#Building metacells from gene expresion matrix (GE)
SC<- SCimplify(GE, k.knn = k.knn, gamma = gamma, n.var.genes = 2000)
supercell_plot(SC$graph.supercells, #network
               main = paste("Metacell network, gamma = ", gamma),
               seed = 1)
#Average expression within each metacell
SC.GE<- supercell_GE(GE, SC$membership)
dim(SC.GE)
row.names(SC.GE)
colnames(SC.GE)
colnames(SC.GE)<- paste0("Metacell", 1:ncol(SC.GE))
colnames(SC.GE)
### 3. Analyze metacells with Seurat
bonemarrow_Metacell <- CreateSeuratObject(counts = SC.GE, project = "bonemarrowMet")
bonemarrow_Metacell<-NormalizeData(bonemarrow_Metacell)
bonemarrow_Metacell<- FindVariableFeatures(bonemarrow_Metacell)
bonemarrow_Metacell<- ScaleData(bonemarrow_Metacell)
bonemarrow_Metacell<- RunPCA(bonemarrow_Metacell, ndims = 50)
ElbowPlot(bonemarrow_Metacell, ndims = 50)
bonemarrow_Metacell<- RunUMAP(bonemarrow_Metacell, dims = 1:30)
bonemarrow_Metacell<- FindNeighbors(bonemarrow_Metacell)
bonemarrow_Metacell<- FindClusters(bonemarrow_Metacell, resolution = 0.2)
DimPlot(bonemarrow_Metacell, label = T)
view(bonemarrow_Metacell@meta.data)
#FEATUREPLOTS cd14+ mono
#FeaturePlot(bonemarrow_Metacell, features = c(), label = T)
#T CELLS
#FeaturePlot(bonemarrow_Metacell, features = c(), label = T)
#B CELLS
#FeaturePlot(bonemarrow_Metacell, features = c(), label = T)

#Rename clusters
#bonemarrow_Metacell<- RenameIdents((bonemarrow_Metacell,'1','2','3',...)
#DimPlot(bonemarrow_Metacell, label=T)
#P1<- DimPlot(bonemarrow_Metacell, group.by = "cell_type", label = T)
#P2<- DimPlot(bonemarrow_Metacell, label =T)
#P1+P2
