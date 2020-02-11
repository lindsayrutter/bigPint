helperMakeHex <- function(datSel, sampleIndex1, sampleIndex2, xbins){

minVal = min(datSel[,-1])
maxVal = max(datSel[,-1])
maxRange = c(minVal, maxVal)
buffer = (maxRange[2]-maxRange[1])/(xbins/2)

seqVec1 <- seq_along(sampleIndex1)
seqVec2 <- seq_along(sampleIndex2)
index1 <- sampleIndex1[seqVec1]
index2 <- sampleIndex2[seqVec2]
x <- unlist(datSel[,index1])
y <- unlist(datSel[,index2])

h <- hexbin(x=x, y=y, xbins=xbins, shape=1, IDs=TRUE, xbnds=maxRange,
ybnds=maxRange)

hexdf <- helperMakeHexDF(h)[["hexdf"]]
clrs <- helperMakeHexDF(h)[["clrs"]]
my_breaks <- helperMakeHexDF(h)[["my_breaks"]]

return(list(hexdf=hexdf, maxRange=maxRange, clrs=clrs, my_breaks=my_breaks,
x=x, y=y))
}