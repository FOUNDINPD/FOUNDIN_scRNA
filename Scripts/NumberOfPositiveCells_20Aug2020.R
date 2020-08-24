#Vikas Bansal Aug 2020
#Calculate the number/percentage of positive cells for a gene (at least one count).


set.seed(786)
setwd("/data/vikas/FOUNDIN_scRNA/")

.libPaths( c( "/home/vikas/Rlib/", .libPaths()) )
library(Seurat)
library(gplots)
library(reshape2)

pathto.outTable <- "/data/vikas/FOUNDIN_scRNA/OutputTable/iPSCsDopaALL/"


outName <- "SinglePositive"



#iPSCsDopa.integrated <- readRDS("OutputData/iPSCsDopaALL/iPSCsDopaALL_integratedAfterCellType.RDS")
DefaultAssay(iPSCsDopa.integrated) <- "RNA"

GeneList <- c("TH", "MAP2")


#Can read any seurat object used in the manuscript
iPSCsDopa.integratedv2 <- iPSCsDopa.integrated 
#iPSCsDopa.integratedv2@meta.data$NewIDsComb <- paste0(iPSCsDopa.integratedv2$SampleID,"_", iPSCsDopa.integratedv2$CellType)
#iPSCsDopa.integratedv2@meta.data$NewIDsComb <- paste0(iPSCsDopa.integratedv2$CellType)
iPSCsDopa.integratedv2@meta.data$NewIDsComb <- paste0(iPSCsDopa.integratedv2$SampleID)

GeneName <- GeneList



for (i in 1:length(GeneName)){
  possibleError <- tryCatch(
    Take_cells1 <- GetAssayData(object = iPSCsDopa.integratedv2, slot = "data")[c(GeneName[i],GeneName[i]),],
    error=function(e) e
  )
  
  #error handling
  if(!inherits(possibleError, "error")){
    
    
    
    Take_cells1 <- GetAssayData(object = iPSCsDopa.integratedv2, slot = "data")[GeneName[i],]
    Take_cells1[Take_cells1 > 0]=1
    
    Take_cells2 <- GetAssayData(object = iPSCsDopa.integratedv2, slot = "data")[GeneName[i],]
    Take_cells2[Take_cells2 > 0]=1
    Take_cells <- Take_cells1 + Take_cells2 
    
    Take_PoscellNames <- (names(which(Take_cells>1)))
  }else(Take_PoscellNames <- list())
  
  
  
  Total_cells <- (data.frame(table(iPSCsDopa.integratedv2@meta.data[,"NewIDsComb"])))
  colnames(Total_cells)[2] <- "TotalNumberofCells"
  
  
  if(length(Take_PoscellNames)>0){Positive_Cells <- (data.frame(table(iPSCsDopa.integratedv2@meta.data[Take_PoscellNames,"NewIDsComb"])))
  colnames(Positive_Cells)[2] <- "PositiveNumberofCells"}else{
    Positive_Cells <- (data.frame(table(iPSCsDopa.integratedv2@meta.data[,"NewIDsComb"])))
    colnames(Positive_Cells)[2] <- "PositiveNumberofCells"
    Positive_Cells[,2] <- 0
  }
  
  Merged_positiveCells <- merge(Positive_Cells,Total_cells, by="Var1",all=T)
  Merged_positiveCells[is.na(Merged_positiveCells)] = 0
  
  Merged_positiveCells$PercentageExpressed <- (Merged_positiveCells$PositiveNumberofCells/Merged_positiveCells$TotalNumberofCells)*100
  colnames(Merged_positiveCells)[1] <- "NewIDsComb"
  colnames(Merged_positiveCells)[c(2,4)] <- paste0(colnames(Merged_positiveCells)[c(2,4)],"_",GeneName[i])
  
  
  if(i==1){Merged_positiveCellsV2 <- Merged_positiveCells}else{
    Merged_positiveCellsV2 <- merge(Merged_positiveCellsV2,Merged_positiveCells)
  } 
  
}


#write.table(Merged_positiveCellsV2, file=paste(pathto.outTable,"Percentage_",outName,"_SamplesInEachCluster.txt"), sep="\t", quote=F, row.names = F)
#write.table(Merged_positiveCellsV2, file=paste(pathto.outTable,"Percentage_",outName,"_OverallCluster.txt"), sep="\t", quote=F, row.names = F)

write.table(Merged_positiveCellsV2, file=paste(pathto.outTable,"Percentage_",outName,"_OverallSamples.txt"), sep="\t", quote=F, row.names = F)




