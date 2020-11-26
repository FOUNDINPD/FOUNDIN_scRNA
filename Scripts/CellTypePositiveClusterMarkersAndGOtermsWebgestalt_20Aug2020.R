#Vikas Bansal 20 Aug 2020
#Foundin project. Cell type name assign and marker identification

set.seed(786)
setwd("/data/vikas/FOUNDIN_scRNA/")

library("WebGestaltR")
library(scclusteval)

iPSCsDopa.integrated <- readRDS("OutputData/iPSCsDopaALL/iPSCsDopaALL_iPSCsDopa.integrated.RDS")


DefaultAssay(iPSCsDopa.integrated) <- "RNA"

pathto.outData <- "/data/vikas/FOUNDIN_scRNA/OutputData/iPSCsDopaALL/ClusterPositiveMarkersCellTypesRes0.2/"


iPSCsDopa.integrated <- RenameIdents(object = iPSCsDopa.integrated, '0' = 'eProg1', '1' = 'lProg2', '2' = 'iDA2', '3' = 'iDA1', '4' = 'iDA4', '5' = 'eProg2', '6' = 'PFPP', '7' = 'lProg1', '8' = 'NE', '9' = 'Ependymal', '10' = 'iDA3')

iPSCsDopa.markers <- FindAllMarkers(iPSCsDopa.integrated, only.pos = TRUE)


saveRDS(iPSCsDopa.markers,file = paste0(pathto.outData,outName,"_iPSCsDopa.PositiveMarkers.RDS"))

All_markers_sig <-iPSCsDopa.markers[(iPSCsDopa.markers$p_val_adj<0.05),]



Markers_eachCluster <- (split(All_markers_sig,f=as.factor(All_markers_sig$cluster)))



for (i in 1:length(Markers_eachCluster)) {
  write.table(Markers_eachCluster[i], file=paste0(pathto.outData,names(Markers_eachCluster)[i], "markers.txt"), sep="\t", row.names=F, quote=F)

  jaTmp <- WebGestaltR(enrichMethod = "ORA", organism = "hsapiens",
                       enrichDatabase = "geneontology_Biological_Process_noRedundant", interestGene = rownames(as.data.frame(Markers_eachCluster[i])), interestGeneType = "genesymbol",
                       referenceSet = "genome_protein-coding", outputDirectory = pathto.outData,
                       projectName = paste0(names(Markers_eachCluster)[i], "markers"))
}


#Broad cell types
iPSCsDopa.integrated <- RenameIdents(object = iPSCsDopa.integrated, 'eProg1' = "Early neuron Progenitor", 'eProg2' = "Early neuron Progenitor",'lProg1' = "Late neuron Progenitor", 'lProg2' = "Late neuron Progenitor", 'iDA2' = "Dopaminergic Neurons", 'iDA3' = "Dopaminergic Neurons", 'iDA1' = "Dopaminergic Neurons", 'iDA4' = "Immature Dopaminergic Neurons",  'PFPP' = "Proliferating Floor Plate Progenitors",  'NE' = "Neuroepithelial-like Cells", 'Ependymal' = "Ependymal-like Cells")
