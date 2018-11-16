#' @title Plot static litre plots
#' 
#' @description Plot static litre plots.
#' 
#' @param data DATA FRAME | Read counts
#' @param dataMetrics LIST | Differential expression metrics
#' @param option CHARACTER STRING ["hexagon" | "allPoints"] | The background of plot; default "hexagon"
#' @param saveFile BOOLEAN [TRUE | FALSE] | Save file to outDir; default TRUE
#' @param outDir CHARACTER STRING | Output directory to save all plots; default current directory
#' @param pointSize INTEGER | Size of plotted points; default 2
#' @param pointColor CHARACTER STRING | Color of gene superimposed on litre plot; default "orange"
#' @param xbins INTEGER | Number of bins partitioning the range of the plot; default 10
#'@param threshVar CHARACTER STRING | Name of column in dataMetrics object that is used to threshold significance; default "FDR"
#' @param threshVal INTEGER | Maximum value to threshold significance from threshVar object; default 0.05
#' @param geneList CHARACTER ARRAY | List of ID values of genes to be drawn from data as parallel coordinate lines. Use this parameter if you have predetermined genes to be drawn. Otherwise, use dataMetrics, threshVar, and threshVal to create genes to be drawn; default NULL
#' 
#' @importFrom dplyr filter select %>%
#' @importFrom GGally ggpairs wrap
#' @importFrom ggplot2 ggplot aes_string aes geom_point xlim ylim geom_hex coord_cartesian xlab ylab geom_ribbon geom_boxplot geom_line geom_abline theme_gray ggtitle scale_fill_manual coord_fixed labs element_text
#' @importFrom grDevices jpeg dev.off
#' @importFrom hexbin hexbin hcell2xy
#' @importFrom htmlwidgets onRender
#' @importFrom plotly plotlyOutput ggplotly renderPlotly layout
#' @importFrom shiny verbatimTextOutput fluidPage reactive renderPrint shinyApp
#' @importFrom stats lm predict
#' @importFrom tidyr gather crossing
#' @importFrom utils str
#' @importFrom Hmisc cut2
#' @importFrom RColorBrewer brewer.pal
#' @export
#' @examples
#' # Example 1: Create one litre plot for each of the 61 genes that have FDR < 1e-10 
#' # and examine the first plot (gene "N_P_Glyma.19G168700.Wm82.a2.v1")
#' 
#' data(soybean_ir_sub)
#' soybean_ir_sub[,-1] <- log(soybean_ir_sub[,-1]+1)
#' data(soybean_ir_sub_metrics)
#' ret <- plotLitre(data = soybean_ir_sub, dataMetrics = soybean_ir_sub_metrics,
#'   threshVal = 1e-10, saveFile = FALSE)
#' length(ret)
#' names(ret)[1]
#' ret[[1]]
#' 
#' # Example 2: Create one litre plot for each of the five most significant genes 
#' # (lowest FDR values) and view the plot for gene"N_P_Glyma.19G168700.Wm82.a2.v1".
#' 
#' geneList = soybean_ir_sub_metrics[["N_P"]][1:5,]$ID
#' ret <- plotLitre(data = soybean_ir_sub, geneList = geneList, pointColor = "deeppink")
#' names(ret)
#' ret[["N_P_Glyma.19G168700.Wm82.a2.v1"]]
#' 
#' # Example 3: Create one litre plot for each of the five most significant genes 
#' # (lowest FDR values) and view the plot for gene "N_P_Glyma.19G168700.Wm82.a2.v1". 
#' # Use points instead of the default hexagons as the background.
#' 
#' ret <- plotLitre(data = soybean_ir_sub, geneList = geneList, pointColor = "deeppink",
#'   option = "allPoints")
#' names(ret)
#' ret[["N_P_Glyma.19G168700.Wm82.a2.v1"]]
#' 

plotLitre = function(data=data, dataMetrics=NULL, outDir=getwd(), pointSize=2, pointColor = "orange", xbins=10, threshVar="FDR", threshVal=0.05, geneList = NULL, saveFile = TRUE, option = "hexagon"){
  
  hexID <- counts <- countColor2 <- ID <- NULL
  
  dat <- data
  rm(data)
  
  colNames <- colnames(dat)
  myPairs <- unique(sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1]))
  myPairs <- myPairs[-which(myPairs=="ID")]
  colGroups <- sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1])
  
  ifelse(!dir.exists(outDir), dir.create(outDir), FALSE)
  
  ret = list()
  for (i in 1:(length(myPairs)-1)){
    for (j in (i+1):length(myPairs)){
      group1 = myPairs[i]
      group2 = myPairs[j]
      
      sampleIndex <- which(sapply(colnames(dat), function(x) unlist(strsplit(x,"[.]"))[1]) %in% c(group1, group2))
      datSel = dat[,c(1, sampleIndex)]
      
      sampleIndex1 <- which(sapply(colnames(datSel), function(x) unlist(strsplit(x,"[.]"))[1]) %in% c(group1))
      sampleIndex2 <- which(sapply(colnames(datSel), function(x) unlist(strsplit(x,"[.]"))[1]) %in% c(group2))
      
      minVal = min(datSel[,-1])
      maxVal = max(datSel[,-1])
      maxRange = c(minVal, maxVal)
      buffer = (maxRange[2]-maxRange[1])/(xbins/2)
      x <- c()
      y <- c()
      for (i in 1:length(sampleIndex1)){
        for (j in 1:length(sampleIndex2)){
          x <- c(x, unlist(datSel[,(sampleIndex1[i])]))
          y <- c(y, unlist(datSel[,(sampleIndex2[j])]))
        }
      }
      
      h <- hexbin(x=x, y=y, xbins=xbins, shape=1, IDs=TRUE, xbnds=maxRange, ybnds=maxRange)
      hexdf <- data.frame (hcell2xy (h),  hexID = h@cell, counts = h@count)
      attr(hexdf, "cID") <- h@cID
      
      # By default, groups into six equal-sized bins
      hexdf$countColor <- cut2(hexdf$counts, g=6, oneval=FALSE)
      hexdf$countColor2 <- as.factor(unlist(lapply(as.character(hexdf$countColor), function(x) substring(strsplit(gsub(" ", "", x, fixed = TRUE), ",")[[1]][1], 2))))
      hexdf$countColor2 <- factor(hexdf$countColor2, levels = as.character(sort(as.numeric(levels(hexdf$countColor2)))))
      
      for (i in 1:(length(levels(hexdf$countColor2))-1)){
        levels(hexdf$countColor2)[i] <- paste0(levels(hexdf$countColor2)[i],"-",levels(hexdf$countColor2)[i+1])
      }
      levels(hexdf$countColor2)[length(levels(hexdf$countColor2))] <- paste0(levels(hexdf$countColor2)[length(levels(hexdf$countColor2))], "+")
    
      my_breaks = levels(hexdf$countColor2)
      clrs <- brewer.pal(length(my_breaks)+3, "Blues")
      clrs <- clrs[3:length(clrs)]
      
      # Generate background hexagon plot for this pair of treatments
      
      if (option == "hexagon"){
        p <- ggplot(hexdf, aes(x=x, y=y, hexID=hexID, counts=counts, fill=countColor2)) + geom_hex(stat="identity") + scale_fill_manual(labels = as.character(my_breaks), values = rev(clrs), name = "Gene count") + geom_abline(intercept = 0, color = "red", size = 0.25) + labs(x = paste0("Read count ", "(", group1, ")"), y = paste0("Read count ", "(", group2, ")")) + theme(axis.text=element_text(size=15), axis.title=element_text(size=15), legend.title=element_text(size=15), legend.text=element_text(size=15)) + coord_fixed(ratio=1)   
      }
      else{
        mainPoints = data.frame(x=x, y=y)
        p <- ggplot(mainPoints, aes(x=x, y=y)) + geom_point(size = pointSize) + geom_abline(intercept = 0, color = "red", size = 0.25) + labs(x = paste0("Read count ", "(", group1, ")"), y = paste0("Read count ", "(", group2, ")")) + theme(axis.text=element_text(size=15), axis.title=element_text(size=15), legend.title=element_text(size=15), legend.text=element_text(size=15)) + coord_fixed(ratio=1)
      }
    
      if (is.null(geneList)){
        rowDEG1 <- which(dataMetrics[[paste0(myPairs[1],"_",myPairs[2])]][threshVar] < threshVal)
        rowDEG2 <- which(dataMetrics[[paste0(myPairs[2],"_",myPairs[1])]][threshVar] < threshVal)
        geneList <- dataMetrics[[paste0(myPairs[1],"_",myPairs[2])]][c(rowDEG1, rowDEG2),1]
      }
      
      for (k in 1:length(geneList)){
        currID = geneList[k]
        currGene = dat %>% filter(ID == currID)
        
        sampleComb = as.data.frame(crossing(as.numeric(currGene[sampleIndex1]), as.numeric(currGene[sampleIndex2])))
        colnames(sampleComb) = c("x", "y")
        
        pg <- p + geom_point(data = sampleComb, aes(x=x, y=y), inherit.aes = FALSE, color = pointColor, size = pointSize) + ggtitle(currGene$ID)
        
        if (saveFile == TRUE){
          jpeg(filename=paste0(outDir, "/", group1, "_", group2, "_", currGene$ID, "_litre.jpg"), height=700, width=1100)
          print(pg)
          dev.off()
        }
        ret[[paste0(group1, "_", group2, "_", currGene$ID)]] <- pg
      }
    }
  }
  invisible(ret)
}
