#Vikas Bansal 10 Oct 2020
#Calculate Pseudotime for FOUNDIN scRNA

set.seed(786)
setwd("/data/vikas/FOUNDIN_scRNA/")
.libPaths( c( "/home/vikas/Rlib2/", .libPaths()) )

library(monocle3)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)

pathto.outPlots <- "/data/vikas/FOUNDIN_scRNA/OutputPlots/iPSCsDopaALL/"
pathto.outData <- "/data/vikas/FOUNDIN_scRNA/OutputData/iPSCsDopaALL/"
outName <- "iPSCsDopaALL_withPseudotime"


iPSCsDopa.integrated <- readRDS("/data/vikas/FOUNDIN_scRNA/OutputData/iPSCsDopaALL/iPSCsDopaALL_integratedAfterBroadCellType.RDS")
DefaultAssay(iPSCsDopa.integrated) <- "integrated"


cds <- as.cell_data_set(iPSCsDopa.integrated)
cds <- cluster_cells(cds)


integrated.sub <- subset(as.Seurat(cds), monocle3_partitions == 1)
cds <- as.cell_data_set(integrated.sub)
cds <- learn_graph(cds)

NeuroepithelialCellsAsBase <- (integrated.sub@meta.data[(integrated.sub@meta.data$CellType == "NE"),"WholeID"])

cds <- order_cells(cds, root_cells = NeuroepithelialCellsAsBase)

saveRDS(cds,file = paste0(pathto.outData,outName,"_Monocle_cds.RDS"))


integrated.sub2 <- as.Seurat(cds, assay = "integrated")
iPSCsDopa.integrated@meta.data <- integrated.sub2@meta.data



DefaultAssay(iPSCsDopa.integrated) <- "RNA"
saveRDS(iPSCsDopa.integrated,file = paste0(pathto.outData,outName,"_integratedAfterMonocleCellType.RDS"))
