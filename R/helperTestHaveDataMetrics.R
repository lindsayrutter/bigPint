helperTestHaveDataMetrics <- function(){
    tryCatch(
        error = function(e){
            message('Be sure you have either dataMetrics or dataSE object defined.')
        }
    )
}
