helperTestHaveData <- function(){
    tryCatch(
        # expr = {
        #     cbind(ID = rownames(x), as.data.frame(assay(x)), 
        #           stringsAsFactors = FALSE)
        # },
        error = function(e){
            message('Be sure you have either data or dataSE object defined.')
        }
    )
}
