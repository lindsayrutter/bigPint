helperSMHex <- function(data, dataMetrics, outDir, pointSize, pointColor,
xbins, threshVar, threshVal, geneList){

counts <- hexID <- ID <- NULL
myPairs <- helperMakePairs(data)[["myPairs"]]
colGroups <- helperMakePairs(data)[["colGroups"]]
ifelse(!dir.exists(outDir), dir.create(outDir), FALSE)

maxVal = max(data[,-1])
minVal = min(data[,-1])
maxRange = c(minVal, maxVal)
xbins=xbins
buffer = (maxRange[2]-maxRange[1])/(xbins/2)

cols.combn <- combn(myPairs, 2, simplify = FALSE)
data_list <- lapply(cols.combn, function(x) {data %>% select(ID,
starts_with(x[1]), starts_with(x[2]))})
names_list <- lapply(cols.combn, function(x) {paste0(x[1], "_", x[2])})

my_fn <- function(data, mapping, ...){
xChar = as.character(mapping$x)[2]
yChar = as.character(mapping$y)[2]
x = data[,c(xChar)]
y = data[,c(yChar)]

h <- hexbin(x=x, y=y, xbins=xbins, shape=1, IDs=TRUE, xbnds=maxRange,
ybnds=maxRange)
hexdf <- data.frame (hcell2xy (h),  hexID = h@cell, counts = h@count)
attr(hexdf, "cID") <- h@cID

p <- ggplot(hexdf, aes(x=x, y=y, fill = counts, hexID=hexID)) +
geom_hex(stat="identity") + geom_abline(intercept = 0, color = "red",
size = 0.25) + coord_cartesian(xlim = c(maxRange[1]-1*buffer,
maxRange[2]+buffer), ylim = c(maxRange[1]-1*buffer, maxRange[2]+buffer))

if (!is.null(geneList)){
    degData <- data[which(rownames(data) %in% geneList),]
    p <- p + geom_point(data = degData, aes_string(x=xChar, y=yChar),
    inherit.aes = FALSE, color = pointColor, size = pointSize)
}
if (!is.null(dataMetrics)){
    colGroups <- helperMakePairs(data)[["colGroups"]]
    group1 = unique(colGroups)[1]
    group2 = unique(colGroups)[2]
    
    rowDEG1 <- which(dataMetrics[[paste0(group1,"_",group2)]]
    [threshVar] < threshVal)
    rowDEG2 <- which(dataMetrics[[paste0(group2,"_",group1)]]
    [threshVar] < threshVal)
    rowDEG <- c(rowDEG1, rowDEG2)
    degID1 <- dataMetrics[[paste0(group1,"_",group2)]][rowDEG,]$ID
    degID2 <- dataMetrics[[paste0(group2,"_",group1)]][rowDEG,]$ID
    degID <- c(degID1, degID2)
    dataID = cbind(ID=rownames(data), data)
    degData <- dataID[which(dataID$ID %in% degID),]
    p <- p + geom_point(data = degData, aes_string(x=xChar, y=yChar),
    inherit.aes = FALSE, color = pointColor, size = pointSize)
}
p
}

my_fn2 <- function(data){
rownames(data) = data$ID
p <- ggpairs(data %>% select(- ID), 
lower = list(continuous = my_fn), 
upper = list(continuous = wrap("cor", size = 4)))
return(p)
}

ret <- lapply(data_list, function(x) my_fn2(x))
names(ret) <- names_list
return(ret)
}  
