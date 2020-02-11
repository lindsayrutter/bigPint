helperMakeHexDF <- function (h){
hexdf <- data.frame (hcell2xy (h),  hexID = h@cell, counts = h@count)
attr(hexdf, "cID") <- h@cID

# By default, groups into six equal-sized bins
hexdf$countColor <- cut2(hexdf$counts, g=6, oneval=FALSE)
hexdf$countColor2 <- as.factor(unlist(lapply(as.character(hexdf$countColor),
function(x) substring(strsplit(gsub(" ", "", x, fixed = TRUE),",")[[1]][1],2))))
hexdf$countColor2 <- factor(hexdf$countColor2,
levels = as.character(sort(as.numeric(levels(hexdf$countColor2)))))

seqVec <- seq_along(levels(hexdf$countColor2)[-1])
newNames <- vapply(seqVec, function(x){paste0(levels(hexdf$countColor2)[x],
"-", levels(hexdf$countColor2)[x+1])}, character(1))
levels(hexdf$countColor2)[seq_along(levels(hexdf$countColor2)[-1])] <- newNames
levels(hexdf$countColor2)[length(levels(hexdf$countColor2))] <-
paste0(levels(hexdf$countColor2)[length(levels(hexdf$countColor2))], "+")

my_breaks = levels(hexdf$countColor2)
clrs <- brewer.pal(length(my_breaks)+3, "Blues")
clrs <- clrs[3:length(clrs)]

return(list(hexdf=hexdf, clrs=clrs, my_breaks=my_breaks))
}