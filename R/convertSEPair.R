convertSEPair <- function(dataSE, group1, group2){
    
    dfFormat <- as.data.frame(rowData(dataSE))
    dataMetrics <- lapply(split.default(dfFormat[-1], sub("\\..*", "",names(dfFormat[-1]))), function(x) cbind(dfFormat[1], setNames(x, sub(".*\\.", "", names(x)))))
    
    oldData <- as.data.frame(assay(dataSE))
    keepDataCol <- which(sapply(colnames(oldData), function(x) strsplit(x, "[.]")[[1]][1]) %in% c(group1, group2))
    keepOldData <- oldData[, keepDataCol]
    
    data <- DelayedArray(keepOldData)
    returnDFPair <- SummarizedExperiment(assays = data, rowData = dataMetrics)
    return(returnDFPair)
}