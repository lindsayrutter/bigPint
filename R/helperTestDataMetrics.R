helperTestDataMetrics <- function(data, dataMetrics, threshVar){

colNames = colnames(data[,-1])
seqVec <- seq_along(colNames)
    
generalMessage = "For more information about formatting the dataMetrics
object, see https://lindsayrutter.github.io/bigPint/articles/dataMetrics.html"

if (!methods::is(dataMetrics, "list")){ 
    stop(paste0("Data metrics object must be of class 'list'. ",
    generalMessage))
}

logicClass = vapply(data[,-1], function(x) methods::is(x, "numeric") ||
methods::is(x, "integer"), logical(length=1))

logicPerl = grep("^[a-zA-Z0-9]+\\.[0-9]+", colNames, perl=TRUE)

if (all(logicPerl == seqVec)){
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
logicThresh <- lapply(dataMetrics, function(x) threshVar %in% colnames(x))

refID = sort(data$ID)
refIDs = lapply(dataMetrics, function(x) all(sort(x[,1]) == refID))

if (all(logicThresh == TRUE)){
    seqVec <- seq_along(metricNames)
    logicThreshQuant <- vapply(seqVec, function(i){
    indexThresh <- which(colnames(dataMetrics[[i]]) %in% threshVar);
    methods::is(dataMetrics[[i]][[indexThresh]],
    "numeric") || methods::is(dataMetrics[[i]][[indexThresh]], "integer")
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
else if (!all(logicThresh == TRUE)){
    stop(paste0("At least one column in each list element in the data
    metrics object should have the same name as the threshVar object. ",
    generalMessage))
}
else if (!all(logicThreshQuant == TRUE)){
    stop(paste0("The column in each list element in the data metrics object
    that has the same name as the threshVar object should be of class
    'numeric' or 'integer'. ", generalMessage))    
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
