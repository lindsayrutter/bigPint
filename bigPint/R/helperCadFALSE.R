helperCadFALSE <- function(cData, dataMetrics, metricPair, aggMethod, nC,
threshVar, threshVal, verbose, vxAxis, saveFile, boxDat, xAxisLabel,
yAxisLabel, lineAlpha, lineSize, plotName, outDir, colList) {
    dendo = cData
    rownames(dendo) = NULL
    d = suppressWarnings(dist(as.matrix(dendo)))
    hC = hclust(d, method=aggMethod)
    k = cutree(hC, k=nC)
    
    seqVec = seq(nC)
    plot_clusters = lapply(seq_along(seqVec), function(j){
        i = rev(order(table(k)))[j]
        x = as.data.frame(cData[which(k==i),])
        x$cluster = "color"
        x$cluster2 = factor(x$cluster)
        x$ID = factor(x$ID)
        xNames = x$ID
        
        if (!is.null(dataMetrics)){
            metricFDR = metricPair[which(as.character(metricPair$ID) %in% 
            xNames),]
            sigID = metricFDR[which(metricFDR[[threshVar]]<=threshVal),]$ID
            xSig = x[which(xNames %in% sigID),]
            xSigNames = rownames(xSig)
            nGenes = nrow(xSig) 
        }
        else{
            xSig = x
            xSigNames = rownames(x)
            nGenes = nrow(x)
        }
        
        if (verbose==TRUE){
            IDs = as.character(xNames)
            saveRDS(IDs, file = paste(outDir, "/", plotName, "_", nC, "_",
            j, ".rds", sep=""))
        }
        
        xSig$ID = xSigNames
        lastTwoIndices = c(ncol(xSig) -1, ncol(xSig))
        xSig = xSig[, -lastTwoIndices]
        
        pcpDat <- melt(xSig, id.vars="ID")
        colnames(pcpDat) <- c("ID", "Sample", "Count")
        pcpDat$Sample <- as.character(pcpDat$Sample)
        pcpDat$ID <- as.factor(pcpDat$ID)
        
        p <- ggplot(boxDat, aes_string(x = 'Sample', y = 'Count')) +
        geom_boxplot() + geom_line(data=pcpDat, aes_string(x = 'Sample',
        y = 'Count', group = 'ID'), colour = colList[j], alpha=lineAlpha,
        size = lineSize) + ylab(yAxisLabel) + xlab(xAxisLabel) +
        ggtitle(paste("Cluster ", j, " Genes (n=", format(nGenes,
        big.mark=",", scientific=FALSE), ")",sep="")) +
        theme(plot.title = element_text(hjust = 0.5, size=14,
        face="plain"), axis.text=element_text(size=11),
        axis.title=element_text(size=14))
        
        if (vxAxis == TRUE){
            p <- p + theme(axis.text.x = element_text(angle=90, hjust=1))
        }
        
        if (verbose==TRUE){
            fileName = paste(outDir, "/", plotName, "_", nC, "_", j, ".jpg",
            sep="")
            jpeg(fileName)
            plot(p)
            invisible(dev.off()) 
        }
        p
    })
    p = arrangeGrob(grobs=plot_clusters, ncol=2)
    
    if (saveFile == TRUE || verbose == TRUE){
        fileName = paste(outDir, "/", plotName, "_", nC, ".jpg", sep="")
        jpeg(fileName)
        grid.draw(p)
        invisible(dev.off())
    }
    return(p)
}