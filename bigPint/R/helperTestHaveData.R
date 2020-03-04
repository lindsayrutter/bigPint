helperTestHaveData <- function(){
    tryCatch(
        error = function(e){
            message('Be sure you have either data or dataSE object defined.')
        }
    )
}
