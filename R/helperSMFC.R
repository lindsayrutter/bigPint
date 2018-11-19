helperSMFC <- function(data, dataMetrics, outDir, pointSize, threshFC,
threshVar, threshVal){

lwr <- upr <- ID <- NULL
rownames(data) <- data$ID
minLine = 0
maxLine = max(data[,-1])
inc = (maxLine-minLine)/100
xv = seq(minLine, maxLine, inc)
uyv = xv*(threshFC+1)
lyv = xv/(threshFC+1)
lineDF = data.frame(xv=xv, uy=uyv, lyv=lyv)
myPairs <- helperMakePairs(data)[["myPairs"]]
colGroups <- helperMakePairs(data)[["colGroups"]]
ifelse(!dir.exists(outDir), dir.create(outDir), FALSE)

k=1
names_list = list()
data_list = list()
seqVec <- seq(1, length(myPairs)-1)
for (i in seq_along(seqVec)){
    for (j in (i+1):length(myPairs)){
        group1 = myPairs[i]
        group2 = myPairs[j]
        datSel <- cbind(ID=data$ID, data[,which(colGroups %in%
        c(group1, group2))])
        data_list[[k]] <- datSel
        names_list[[k]] <- paste0(group1,"_",group2)
        k = k +1
    }
}

my_fn <- function(data, mapping, ...){
xChar = as.character(mapping$x)[2]
yChar = as.character(mapping$y)[2]
x = data[,c(xChar)]
y = data[,c(yChar)]
indexPoints=c()
seqVec = seq(1, length(x))
for (i in seq_along(seqVec)){
    fract = x[i]/y[i]
    if (!is.nan(fract)){
        if(fract > (threshFC + 1) || fract < (1/(threshFC+1))){
            indexPoints = c(indexPoints, i)
        }
    }
}
plotPoints = data[indexPoints,]

colGroups <- helperMakePairs(data)[["colGroups"]]
group1 = unique(colGroups)[1]
group2 = unique(colGroups)[2]

if (!is.null(dataMetrics)){
    rowDEG1 <- which(dataMetrics[[paste0(group1,"_",group2)]][threshVar]
    < threshVal)
    rowDEG2 <- which(dataMetrics[[paste0(group2,"_",group1)]][threshVar]
    < threshVal)
    rowDEG <- c(rowDEG1, rowDEG2)
    degID1 <- dataMetrics[[paste0(group1,"_",group2)]][rowDEG,]$ID
    degID2 <- dataMetrics[[paste0(group2,"_",group1)]][rowDEG,]$ID
    degID <- c(degID1, degID2)
    degData <- data[which(rownames(data) %in% degID),]
    indexBoth = rownames(plotPoints) %in% rownames(degData)
    indexBlue = rownames(degData) %in% rownames(plotPoints)
    redPoints = plotPoints[indexBoth,]
    greyPoints = plotPoints[!indexBoth,]
    bluePoints = degData[!indexBlue,]
    
    p <- ggplot(lineDF, aes(x=xv, y=lyv)) + geom_line(aes(y = lyv), alpha=0.1) +
    geom_line(aes(y = uyv), alpha=0.1) + geom_ribbon(aes(ymin=lyv,ymax=uyv),
    fill="blue", alpha="0.3") + geom_point(data = bluePoints,
    aes_string(x = xChar, y = yChar), size = pointSize, alpha=0.5,
    color = "blue")+ geom_point(data = greyPoints,
    aes_string(x=xChar, y = yChar), size=pointSize, color = "darkgrey") +
    geom_point(data = redPoints, aes_string(x = xChar, y = yChar),
    size = pointSize, color = "red")
}
else{
    p <- ggplot(lineDF, aes(x=xv, y=lyv)) + geom_line(aes(y = lyv), alpha=0.1) +
    geom_line(aes(y = uyv), alpha=0.1) + geom_ribbon(aes(ymin=lyv,ymax=uyv),
    fill="blue", alpha="0.3") + geom_point(data = plotPoints,
    aes_string(x = xChar, y = yChar), size=pointSize)  
}
p
}

my_fn2 <- function(data){
    p <- ggpairs(data %>% select(- ID), lower = list(continuous = my_fn),
    upper = list(continuous = wrap("cor", size = 4)))
    return(p)
}

ret <- lapply(data_list, function(x) my_fn2(x))
names(ret) <- names_list
return(ret)
} 
