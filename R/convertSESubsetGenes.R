#' @title Convert SummarizedExperiment to only contain defined genes
#' 
#' @description Reduce a SummarizedExperiment object so that it only now
#' contains a subset of genes.
#' 
#' @param dataSE SUMMARIZEDEXPERIMENT | Summarized experiment format that
#' can be used in lieu of data and dataMetrics
#' @param geneList CHARACTER ARRAY | List of gene IDs to remain in dataSE
#' @return A new dataSE object that is a subset of the input dataSE in that it
#' now only contains the user-specified list of genes.
#' @export
#' @examples
#' # Example: Read in example SummarizedExperiment object that originally
#' # contains 5604 genes. Reduce it to now only contain ten genes (specifically
#' # the ones with the lowest FDR).
#' 
#' data(se_soybean_ir_sub)
#' geneList <- as.data.frame(rowData(se_soybean_ir_sub)) %>%
#'     arrange(N_P.FDR) %>% filter(row_number() <= 10)
#' geneList <- geneList[,1]
#' se_soybean_ir_sub_2 <- convertSESubsetGenes(se_soybean_ir_sub, geneList)
#' 
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