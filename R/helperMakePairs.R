helperMakePairs <- function(data){
colNames <- colnames(data)

seqVec <- seq_along(colNames)
colGroups <- vapply(seqVec, function(i){
    strsplit(colNames[i],"[.]")[[1]][1]
}, character(1))

myPairs <- unique(colGroups)
if ("ID" %in% colNames){
    myPairs <- myPairs[-which(myPairs=="ID")] 
}
return(list(myPairs = myPairs, colGroups = colGroups))
}