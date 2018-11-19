helperMakePairs <- function(data){
colNames <- colnames(data)
colGroups <- c()
seqVec <- seq(1,length(colNames))
for (i in seq_along(seqVec)){
    colGroups[i] <- strsplit(colNames[i],"[.]")[[1]][1]
}
myPairs <- unique(colGroups)
if ("ID" %in% colNames){
    myPairs <- myPairs[-which(myPairs=="ID")] 
}
return(list(myPairs = myPairs, colGroups = colGroups))
}