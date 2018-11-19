helperMakeHex <- function(datSel, sampleIndex1, sampleIndex2, xbins){

minVal = min(datSel[,-1])
maxVal = max(datSel[,-1])
maxRange = c(minVal, maxVal)
buffer = (maxRange[2]-maxRange[1])/(xbins/2)
x <- c()
y <- c()
seqVec <- seq(1, length(sampleIndex1))
for (i in seq_along(seqVec)){
    seqVec <- seq(1, length(sampleIndex2))
    for (j in seq_along(seqVec)){
        x <- c(x, unlist(datSel[,(sampleIndex1[i])]))
        y <- c(y, unlist(datSel[,(sampleIndex2[j])]))
    }
}

h <- hexbin(x=x, y=y, xbins=xbins, shape=1, IDs=TRUE, xbnds=maxRange,
ybnds=maxRange)

hexdf <- helperMakeHexDF(h)[["hexdf"]]
clrs <- helperMakeHexDF(h)[["clrs"]]
my_breaks <- helperMakeHexDF(h)[["my_breaks"]]

return(list(hexdf=hexdf, maxRange=maxRange, clrs=clrs, my_breaks=my_breaks,
x=x, y=y))
}