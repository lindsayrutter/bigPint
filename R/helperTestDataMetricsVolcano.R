helperTestData <- function(data, dataMetrics, threshVar, PValue, logFC){
  
  generalMessage = "For more information about formatting the dataMetrics
  objects, see https://lrutter.github.io/bigPint/articles/dataMetrics.html.
  Note that volcano plots require that each element in the dataMetrics object
  has additional two columns, a PValue column and a logFC column."
  
  colNames = colnames(data[,-1])
  
  logicClass = vapply(data[,-1], function(x) class(x) %in% c("numeric",
  "integer"), logical(length=1))

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
  logicThresh <- lapply(dataMetrics, function(x) threshVar %in% colnames(x))
  logicPValue <- lapply(dataMetrics, function(x) PValue %in% colnames(x))
  logicFC <- lapply(dataMetrics, function(x) logFC %in% colnames(x))
    
  if (all(logicThresh == TRUE) && all(logicPValue == TRUE) && all(logicFC == TRUE)){
    logicThreshQuant = c()
    logicPValueQuant = c()
    logicFCQuant = c()
    seqVec <- seq(1,length(metricNames))
    for (i in seq_along(seqVec)){
      indexThresh <- which(colnames(dataMetrics[[i]]) %in% threshVar)
      logicThreshQuant[i] <- class(dataMetrics[[i]][[indexThresh]]) %in%
      c("numeric", "integer")
      indexPValue <- which(colnames(dataMetrics[[i]]) %in% PValue)
      logicPValueQuant[i] <- class(dataMetrics[[i]][[indexPValue]]) %in%
      c("numeric", "integer")
      indexFC <- which(colnames(dataMetrics[[i]]) %in% logFC)
      logicFCQuant[i] <- class(dataMetrics[[i]][[indexFC]]) %in%
      c("numeric", "integer")
    }
  }
  
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
    stop(paste0("Data metrics object must be of class 'list'. ", generalMessage))
  }
  else if (length(dataMetrics) != combnMetrics){
    stop(paste0("There should be ", combnMetrics, "list elements
    in the data metrics object to represent each pairwise combination of the "
    , nGroups, "treatment groups in the data object. ", generalMessage))
  }
  else if (!all(logicDF == TRUE)){
    stop(paste0("Each list element in data metrics object must be of class
    'data.frame'. ", generalMessage))
  }
  else if (!all(logicID == TRUE)){
    stop(paste0("The first column of each list element in the data metrics
    object must be called 'ID'. ", generalMessage))
  }
  else if (!all(logicIDChar == TRUE)){
    stop(paste0("The first column of each list element in the data metrics
    object must be of class 'character'. ", generalMessage))
  }  
  else if (!all(logicIDUniq == TRUE)){
    stop(paste0("The first column of each list element in the data metrics
    object must contain unique items. ", generalMessage))
  }   
  else if (!all(logicListName == TRUE)){
    stop(paste0("The name of each list element in the data metrics object
    must match the Perl expression ^[a-zA-Z0-9]+_[a-zA-Z0-9]+. ", generalMessage))
  }   
  else if (numListName != combnMetrics){
    stop(paste0("The name of each list element in the data metrics object must
    match the Perl expression ^[a-zA-Z0-9]+_[a-zA-Z0-9]+. Each pattern [a-zA-Z0-9]
    should be the alphanumeric name of a treatment group in the data object. ",
    generalMessage))
  }
  else if (!all(logicListName == TRUE)){
    stop(paste0("The name of at least one of the list elements in the data metrics
    object repeats the same treatment group name on both sides of the underscore
    (for example: 'A_A'). The names of each list element in the data metrics
    object should have different treatment groups names on both sides of the
    underscore (for example: 'A_B'). ", generalMessage))
  }
  else if (!all(logicThresh == TRUE)){
    stop(paste0("At least one column in each list element in the data metrics object
    should have the same name as the threshVar object. ", generalMessage))
  }
  else if (!all(logicFC == TRUE)){
    stop(paste0("For volcano plots, at least one column in each list element in
    the data metrics object should have the same name as the logFC object. ",
    generalMessage))
  }
  else if (!all(logicPValue == TRUE)){
    stop(paste0("For volcano plots, at least one column in each list element in
    the data metrics object should have the same name as the PValue object. ",
    generalMessage))
  }
  else if (!all(logicThreshQuant == TRUE)){
    stop(paste0("The column in each list element in the data metrics object that has
    the same name as the threshVar object should be of class 'numeric' or 
    'integer'. ", generalMessage))    
  }
  else if (!all(logicFCQuant == TRUE)){
    stop(paste0("For volcano plots, the column in each list element in the data
    metrics object that has the same name as the logFC object should be of class
    'numeric' or'integer'. ", generalMessage))    
  }
  else if (!all(logicPValueQuant == TRUE)){
    stop(paste0("For volcano plots, the column in each list element in the data
    metrics object that has the same name as the PValue object should be of class
    'numeric' or 'integer'. ", generalMessage))    
  }
}
