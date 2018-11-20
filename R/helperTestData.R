helperTestData <- function(data){
  
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
  else if (!all(logicClass == TRUE)){
    stop("Error: All columns but the first must be of class 'integer' or
         'numeric' in the data object.")
  }
  else if (all(logicPerl != c(1:length(colNames)))){
    stop("Error: In the data object, the names of all columns but the first
         must match the Perl expression '^[a-zA-Z0-9]+\\.[0-9]+'")
  }
  else if (uGroups < 2){
    stop("Error: There must be at least two treatment groups in the
    data object.")
  }
  else if (!all(logicReps == TRUE)){
    step("Error: Each treatment group must ahve at least two replicates
    in the dat object.")
  }
}