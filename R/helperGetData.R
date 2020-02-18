helperGetData <- function(x){
    tryCatch(
        expr = {
            cbind(ID = rownames(x), as.data.frame(assay(x)), 
                  stringsAsFactors = FALSE)
        },
        error = function(e){
            message('Be sure the dataSE object has rownames with gene IDs
            and assay() functions from SummarizedExperiment.')
        }
    )
    cbind(ID = rownames(x), as.data.frame(assay(x)), 
          stringsAsFactors = FALSE)
}
