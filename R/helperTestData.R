helperTestData <- function(data){
  
  generalMessage = "For more information about formatting the dataMetrics
  objects, see https://lrutter.github.io/bigPint/articles/data.html"
  
  logicClass = vapply(data[,-1], function(x) class(x) %in% c("numeric",
  "integer"), logical(length=1))
  
  colNames = colnames(data[,-1])

  logicPerl = grep("^[a-zA-Z0-9]+\\.[0-9]+", colNames, perl=TRUE)
  
  if (all(logicPerl == c(1:length(colNames)))){
    colGroups = c()
    colReps = c()
    seqVec <- seq(1,length(colNames))
    for (i in seq_along(seqVec)){
      colGroups[i] <- strsplit(colNames[i],"[.]")[[1]][1]
      colReps[i] <- strsplit(colNames[i],"[.]")[[1]][2]
    }
    uGroups = unique(colGroups)
    nGroups = length(unique(colGroups))
    
    logicReps = vapply(uGroups, function(x) length(which(colGroups %in% x))>1,
    logical(length=1))
  }
  
  if (!class(data) == "data.frame"){ 
    stop(paste0("Data object must be of class 'data.frame'. ", generalMessage))
  }
  else if (colnames(data)[1] != "ID"){
    stop(paste0("First column of data object must be called 'ID'. ", generalMessage))
  }
  else if (!class(data[,1]) == "character"){
    stop(paste0("First column of data object must be of class 'character'." , generalMessage))    
  }
  else if (anyDuplicated(data[,1])>0){
    stop(paste0("First column of data object must contain unique IDs. ", generalMessage))
  }
  else if (ncol(data) < 5){
    stop(paste0("There must be at least five columns in the data object.
         At a minimum, there must be one column for IDs, two columns for
         replicates of one treatment group, and two columns for replicates
         of a second treatment group. ", generalMessage)) 
  }
  else if (!all(logicClass == TRUE)){
    stop(paste0("All columns but the first must be of class 'integer' or
         'numeric' in the data object. ", generalMessage))
  }
  else if (all(logicPerl != c(1:length(colNames)))){
    stop(paste0("In the data object, the names of all columns but the first
         must match the Perl expression '^[a-zA-Z0-9]+\\.[0-9]+'", generalMessage))
  }
  else if (length(uGroups) < 2){
    stop(paste0("There must be at least two treatment groups in the
    data object. ", generalMessage))
  }
  else if (!all(logicReps == TRUE)){
    stop(paste0("Each treatment group must have at least two replicates
    in the data object. ", generalMessage))
  }
}
