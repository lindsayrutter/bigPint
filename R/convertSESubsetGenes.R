convertSESubsetGenes <- function(dataSE, geneList){
    
    dfFormat <- as.data.frame(rowData(dataSE))
    keepDF <- which(dfFormat$ID %in% geneList)
    dataMetrics <- dfFormat[keepDF, ]
    
    oldData <- as.data.frame(assay(dataSE))
    keepDataCol <- which(rownames(oldData) %in% geneList)
    keepOldData <- oldData[keepDataCol, ]
     
    data <- DelayedArray(keepOldData)
    returnDFPair <- SummarizedExperiment(assays = data, rowData = dataMetrics)
    return(returnDFPair)
}