helperTestData <- function(data){
    
colNames = colnames(data[,-1])
seqVec <- seq_along(colNames)
    
generalMessage = "For more information about formatting the data
object, see https://lindsayrutter.github.io/bigPint/articles/data.html"

if (!methods::is(data, "data.frame")){ 
    stop(paste0("Data object must be of class 'data.frame'. If you are
    using a SummarizedExperiment input data, be sure to defined it as
    the 'dataSE' object. ",
    generalMessage))
}

# Convert to data.frame in case tibble
data <- as.data.frame(data)

logicClass = vapply(data[,-1], function(x) methods::is(x, "numeric") ||
methods::is(x, "integer"), logical(length=1))

logicPerl = grep("^[a-zA-Z0-9]+\\.[0-9]+", colNames, perl=TRUE)

if (!all(logicPerl == seqVec)){
    stop(paste0("In the data object, the names of all columns but the first
    must match the Perl expression '^[a-zA-Z0-9]+\\.[0-9]+'", generalMessage))
}

if (all(logicPerl == seqVec)){
    colGroups <- vapply(seqVec, function(i){
        strsplit(colNames[i],"[.]")[[1]][1]
    }, character(1))
    colReps <- vapply(seqVec, function(i){
        strsplit(colNames[i],"[.]")[[1]][2]
    }, character(1))
    
    
    uGroups = unique(colGroups)
    nGroups = length(unique(colGroups))
    
    logicReps = vapply(uGroups, function(x) length(which(colGroups %in% x))>1,
    logical(length=1))
}

if (colnames(data)[1] != "ID"){
    stop(paste0("First column of data object must be called 'ID'. ",
    generalMessage))
}
else if (!methods::is(data[,1], "character")){
    stop(paste0("First column of data object must be of class 'character'. ",
    generalMessage))    
}
else if (anyDuplicated(data[,1])>0){
    stop(paste0("First column of data object must contain unique IDs. ",
    generalMessage))
}
else if (ncol(data) < 5){
    stop(paste0("There must be at least five columns in the data object.
    At a minimum, there must be one column for IDs, two columns for
    replicates of one treatment group, and two columns for replicates of a
    second treatment group. ", generalMessage))
}
else if (!all(logicClass == TRUE)){
    stop(paste0("All columns but the first must be of class 'integer' or
    'numeric' in the data object. ", generalMessage))
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
