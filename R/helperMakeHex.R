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
hexdf <- data.frame (hcell2xy (h),  hexID = h@cell, counts = h@count)
attr(hexdf, "cID") <- h@cID

# By default, groups into six equal-sized bins
hexdf$countColor <- cut2(hexdf$counts, g=6, oneval=FALSE)
hexdf$countColor2 <- as.factor(unlist(lapply(as.character(hexdf$countColor),
function(x) substring(strsplit(gsub(" ", "", x, fixed = TRUE),",")[[1]][1],2))))
hexdf$countColor2 <- factor(hexdf$countColor2,
levels = as.character(sort(as.numeric(levels(hexdf$countColor2)))))

seqVec <- seq(1, length(levels(hexdf$countColor2))-1)
for (i in seq_along(seqVec)){
    levels(hexdf$countColor2)[i] <- paste0(levels(hexdf$countColor2)[i],
    "-", levels(hexdf$countColor2)[i+1])
}
levels(hexdf$countColor2)[length(levels(hexdf$countColor2))] <-
paste0(levels(hexdf$countColor2)[length(levels(hexdf$countColor2))], "+")

my_breaks = levels(hexdf$countColor2)
clrs <- brewer.pal(length(my_breaks)+3, "Blues")
clrs <- clrs[3:length(clrs)]

return(list(hexdf=hexdf, maxRange=maxRange, clrs=clrs, my_breaks=my_breaks,
x=x, y=y))
}