#' @title Convert SummarizedExperiment to contain defined treatment pair
#' 
#' @description Reduce a SummarizedExperiment object that initially contains
#' more than two treatment groups to now only contain a user-specified subset
#' pair of treatment groups. Note that this function is only necesary for
#' users using the SummarizedExperiment input object (instead of the
#' combination of data and dataMetrics input objects.)
#' 
#' @param dataSE SUMMARIZEDEXPERIMENT | Summarized experiment format that
#' can be used in lieu of data and dataMetrics
#' @param group1 CHARACTER STRING | Name of one treatment group that will
#' remain in the dataSE object
#' @param group2 CHARACTER STRING | Name of second treatment group that will
#' remain in the dataSE object
#' @importFrom SummarizedExperiment rowData assay SummarizedExperiment
#' @importFrom DelayedArray DelayedArray
#' @importFrom stats setNames
#' @return A new dataSE object that is a subset of the input dataSE in that it
#' now only contains the user-specified pair of treatment groups.
#' @export
#' @examples
#' # Example: Read in example SummarizedExperiment object that contains three
#' # treatment groups (S1, S2, and S3). Reduce it to now only contain two
#' # treatment groups (S1 and S3).
#' 
#' data(se_soybean_cn_sub)
#' se_soybean_cn_sub_2 <- convertSEPair(se_soybean_cn_sub, "S1", "S3")
#' 
convertSEPair <- function(dataSE, group1, group2){
    
    dfFormat <- as.data.frame(rowData(dataSE))
    dataMetrics <- lapply(split.default(dfFormat[-1], sub("\\..*", "",
        names(dfFormat[-1]))), function(x) cbind(dfFormat[1],
        setNames(x, sub(".*\\.", "", names(x)))))
    
    oldData <- as.data.frame(assay(dataSE))
    keepDataCol <- which(vapply(colnames(oldData), function(x)
        strsplit(x, "[.]")[[1]][1]) %in% c(group1, group2))
    keepOldData <- oldData[, keepDataCol]
    
    data <- DelayedArray(keepOldData)
    returnDFPair <- SummarizedExperiment(assays = data)
    
    if (length(dataMetrics) > 0){    
        
        listNames <- names(dataMetrics)
        keepList <- which(lapply(listNames, function(x)
            strsplit(x, "[.]")[[1]][1]) %in% c(paste0(group1, "_", group2)))
        dataMetrics = dataMetrics[keepList]
    
        for (k in seq_len(length(dataMetrics))){
            colnames(dataMetrics[[k]])[1] = "ID"   
        }
        returnDFPair <- SummarizedExperiment(assays = data,
            rowData = dataMetrics)
    }
    
    return(returnDFPair)
}
