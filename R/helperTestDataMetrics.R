helperTestData <- function(data, dataMetrics){
  
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
  
  metricNames = names(dataMetrics)
  combnMetrics = (nGroups * (nGroups-1))/2
  logicDF <- lapply(dataMetrics, function(x) class(x) == "data.frame")
  logicID <- lapply(dataMetrics, function(x) colnames(x)[1] == "ID")
  logicIDChar <- lapply(dataMetrics, function(x) class(x[,1]) == "character")
  logicIDUniq <- lapply(dataMetrics, function(x) anyDuplicated(x[,1])>0)
  logicListName = grep("^[a-zA-Z0-9]+_[a-zA-Z0-9]+", metricNames, perl=TRUE)
  
  if (all(logicListName == seq(1:length(metricNames)))){
    metric1 = c()
    metric2 = c()
    metricNotSame = c()
    seqVec <- seq(1,length(metricNames))
    for (i in seq_along(seqVec)){
      metric1[i] <- strsplit(metricNames[i],"[_]")[[1]][1]
      metric2[i] <- strsplit(metricNames[i],"[_]")[[1]][2]
      metricNotSame[i] <- (metric1[i] != metric2[i])
    }
  }
  metric12 = c(metric1, metric2)
  metrict = table(metric12)
  numListName = sum(metrict==(combnMetrics-1))
  
  if (!class(dataMetrics) == "list"){ 
    stop("Data metrics object must be of class 'list'.")
  }
  else if (length(dataMetrics) != combnMetrics){
    stop(paste0("There should be ", combnMetrics, "list elements
    in the data metrics object to represent each pairwise combination of the "
    , nGroups, "treatment groups in the data object."))
  }
  else if (!all(logicDF == TRUE)){
    stop("Each list element in data metrics object must be of class 'data.frame'.")    
  }
  else if (!all(logicID == TRUE)){
    stop("The first column of each list element in the data metrics object must be
    called 'ID'.") 
  }
  else if (!all(logicIDChar == TRUE)){
    stop("The first column of each list element in the data metrics object must be
    of class 'character'.") 
  }  
  else if (!all(logicIDUniq == TRUE)){
    stop("The first column of each list element in the data metrics object must
    contain unique items.") 
  }   
  else if (!all(logicListName == TRUE)){
    stop("The name of each list element in the data metrics object must match the
    Perl expression ^[a-zA-Z0-9]+_[a-zA-Z0-9]+.") 
  }   
  else if (numListName != combnMetrics){
    stop("The name of each list element in the data metrics object must match the
    Perl expression ^[a-zA-Z0-9]+_[a-zA-Z0-9]+. Each pattern [a-zA-Z0-9]
    should be the alphanumeric name of a treatment group in the data object.")    
  }
  else if (!all(logicListName == TRUE)){
    stop("The name of at least one of the list elements in the data metrics
    object repeats the same treatment group name on both sides of the underscore
    (for example: 'A_A'). The names of each list elements in the data metrics
    object should have different treatment gruops names on both sides of the
    underscore (for example: 'A_B'.")    
  }

}
