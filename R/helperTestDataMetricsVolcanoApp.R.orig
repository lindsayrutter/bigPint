#' @importFrom methods is
helperTestDataMetricsVolcanoApp <- function(data, dataMetrics, PValue,
logFC){

generalMessage = "For more information about formatting the dataMetrics
objects, see https://lindsayrutter.github.io/bigPint/articles/dataMetrics.html.
Note that volcano plots require that each element in the dataMetrics object
has additional two columns, a PValue column and a logFC column."

if (!methods::is(dataMetrics, "list")){ 
    stop(paste0("Data metrics object must be of class 'list'. ",
    generalMessage))
}

colNames = colnames(data[,-1])
seqVec <- seq_along(colNames)

logicClass = vapply(data[,-1], function(x) methods::is(x, "numeric") ||
methods::is(x, "integer"), logical(length=1))

logicPerl = grep("^[a-zA-Z0-9]+\\.[0-9]+", colNames, perl=TRUE)

if (all(logicPerl == seq_along(colNames))){
    colGroups <- vapply(seqVec, function(i){
        strsplit(colNames[i],"[.]")[[1]][1]
    }, character(1))
    colReps <- vapply(seqVec, function(i){
        strsplit(colNames[i],"[.]")[[1]][2]
    }, character(1))
    uGroups = unique(colGroups)
    nGroups = length(unique(colGroups))
    logicReps = vapply(uGroups, function(x)
    length(which(colGroups %in% x))>1,
    logical(length=1))
}

metricNames = names(dataMetrics)
combnMetrics = (nGroups * (nGroups-1))/2
logicDF <- lapply(dataMetrics, function(x) methods::is(x, "data.frame"))

if (!all(logicDF == TRUE)){
    stop(paste0("Each list element in data metrics object must be of class
    'data.frame'. ", generalMessage))
}

logicID <- lapply(dataMetrics, function(x) colnames(x)[1] == "ID")
logicIDChar <- lapply(dataMetrics, function(x) methods::is(x[,1],
"character"))
logicIDDup <- lapply(dataMetrics, function(x) anyDuplicated(x[,1])>0)
logicListName = grep("^[a-zA-Z0-9]+_[a-zA-Z0-9]+", metricNames, perl=TRUE)
logicPValue <- lapply(dataMetrics, function(x) PValue %in% colnames(x))
logicFC <- lapply(dataMetrics, function(x) logFC %in% colnames(x))

refID = sort(data$ID)
refIDs = lapply(dataMetrics, function(x) all(sort(x[,1]) == refID))

if (all(logicPValue == TRUE) && all(logicFC == TRUE)){
    seqVec <- seq_along(metricNames)
    
    logicPValueQuant <- vapply(seqVec, function(i){
        indexPValue <- which(colnames(dataMetrics[[i]]) %in% PValue);
        methods::is(dataMetrics[[i]][[indexPValue]],
        "numeric") || methods::is(dataMetrics[[i]][[indexPValue]], "integer")
    }, logical(1))
    
    logicFCQuant <- vapply(seqVec, function(i){
        indexFC <- which(colnames(dataMetrics[[i]]) %in% logFC);
        methods::is(dataMetrics[[i]][[indexFC]],
        "numeric") || methods::is(dataMetrics[[i]][[indexFC]], "integer")
    }, logical(1))
}

if (all(logicListName == seq_along(metricNames))){
    seqVec <- seq_along(metricNames)
    metric1 <- vapply(seqVec, function(i){
        strsplit(metricNames[i],"[_]")[[1]][1]
    }, character(1))
    metric2 <- vapply(seqVec, function(i){
        strsplit(metricNames[i],"[_]")[[1]][2]
    }, character(1))
    metricNotSame <- vapply(seqVec, function(i){
        metric1[i] != metric2[i]
    }, logical(1))
}
metric12 = c(metric1, metric2)
metrict = table(metric12)
ddMSame = sort(unique(metric12)) == sort(uGroups)
numListName = sum(metrict==(nGroups-1))

if (length(dataMetrics) != combnMetrics){
    stop(paste0("There should be ", combnMetrics, " list elements in the
    data metrics object to represent each pairwise combination of the ",
    nGroups, " treatment groups in the data object. ", generalMessage))
}
else if (!all(logicID == TRUE)){
    stop(paste0("The first column of each list element in the data metrics
    object must be called 'ID'. ", generalMessage))
}
else if (!all(logicIDChar == TRUE)){
    stop(paste0("The first column of each list element in the data metrics
    object must be of class 'character'. ", generalMessage))
}  
else if (!all(logicIDDup != TRUE)){
    stop(paste0("The first column of each list element in the data metrics
    object must contain unique items. ", generalMessage))
}   
else if (length(metricNames) != length(logicListName)){
    stop(paste0("The name of each list element in the data metrics object
    must match the Perl expression ^[a-zA-Z0-9]+_[a-zA-Z0-9]+. ",
    generalMessage))
}   
else if (numListName != nGroups){
    stop(paste0("The name of each list element in the data metrics object
    must match the Perl expression ^[a-zA-Z0-9]+_[a-zA-Z0-9]+. Each pattern
    [a-zA-Z0-9] should be the alphanumeric name of a treatment group in the
    data object. ", generalMessage))
}
else if (!all(ddMSame == TRUE)){
    stop(paste0("The names of the list elements in the data metrics object
    must include the treatment groups from the data object. ",
    generalMessage))    
}
else if (!all(logicFC == TRUE)){
    stop(paste0("For volcano plots, at least one column in each list element
    in the data metrics object should have the same name as the logFC
    object. ", generalMessage))
}
else if (!all(logicPValue == TRUE)){
    stop(paste0("For volcano plots, at least one column in each list element
    in the data metrics object should have the same name as the PValue
    object. ", generalMessage))
}
else if (!all(logicFCQuant == TRUE)){
    stop(paste0("For volcano plots, the column in each list element in the
    data metrics object that has the same name as the logFC object should be
    of class 'numeric' or'integer'. ", generalMessage))    
}
else if (!all(logicPValueQuant == TRUE)){
    stop(paste0("For volcano plots, the column in each list element in the
    data metrics object that has the same name as the PValue object should
    be of class 'numeric' or 'integer'. ", generalMessage))    
}
else if (!all(refIDs == TRUE)){
    stop(paste0("The ID column in each list element in the data metrics
    object must contain the same IDs (regardless of order) as the ID column
    in the data object. ", generalMessage))
}
else if (!all(metricNotSame)){
    stop(paste0("The name of at least one of the list elements in the data
    metrics object repeats the same treatment group name on both sides of
    the underscore (for example: 'A_A'). The names of each list element in
    the data metrics object should have different treatment groups names on
    both sides of the underscore (for example: 'A_B'). ", generalMessage))
}
}
