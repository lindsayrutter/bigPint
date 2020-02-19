convertSEPair <- function(dataSE, group1, group2){
    
    dfFormat <- as.data.frame(rowData(dataSE))
    dataMetrics <- lapply(split.default(dfFormat[-1], sub("\\..*", "",names(dfFormat[-1]))), function(x) cbind(dfFormat[1], setNames(x, sub(".*\\.", "", names(x)))))
    
    oldData <- as.data.frame(assay(dataSE)) #changed from as.data.frame(assay(dataSE))
    keepDataCol <- which(sapply(colnames(oldData), function(x) strsplit(x, "[.]")[[1]][1]) %in% c(group1, group2))
    keepOldData <- oldData[, keepDataCol]
    
    data <- DelayedArray(keepOldData)
    returnDFPair <- SummarizedExperiment(assays = data)
    
    if (length(dataMetrics) > 0){    
        
        listNames <- names(dataMetrics)
        keepList <- which(sapply(listNames, function(x) strsplit(x, "[.]")[[1]][1]) %in% c(paste0(group1, "_", group2)))
        dataMetrics = dataMetrics[keepList]
    
        for (k in 1:length(dataMetrics)){
            colnames(dataMetrics[[k]])[1] = "ID"   
        }
        returnDFPair <- SummarizedExperiment(assays = data, rowData = dataMetrics)
    }
    
    return(returnDFPair)
}