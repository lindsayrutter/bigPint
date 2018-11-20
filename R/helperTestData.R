helperTestData <- function(data){
  
  logicArray = vapply(data[,-1], function(x) class(x) %in% c("numeric",
  "integer"), logical(length=1))
  
  if (!class(data) == "data.frame"){ 
    stop("Error: Data object must be data frame.")
  }
  else if (colnames(data)[1] != "ID"){
    stop("Error: First column of data object must be called 'ID'.")
  }
  else if (!class(data[,1]) == "character"){
    stop("Error: First column of data object must be of type 'character'.")    
  }
  else if (anyDuplicated(data[,1])>0){
    stop("Error: First column of data object must contain unique IDs.") 
  }
  else if (ncol(data) < 5){
    stop("Error: There must be at least five columns in the data object.
         At a minimum, there must be one column for IDs, two columns for
         replicates of one treatment group, and two columns for replicates
         of a second treatment group.") 
  }
  else if (all(logicArray == FALSE)){
    stop("Error: All columns but the first must be of class integer or
         numeric in the data object.")
  }
  
}