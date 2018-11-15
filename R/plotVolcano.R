#' @title Plot static volcano plot
#' 
#' @description Plot static volcano plot.
#' 
#' @param data DATA FRAME | Read counts
#' @param dataMetrics LIST | Differential expression metrics
#' @param option CHARACTER STRING ["hexagon" | "allPoints"] | The background of plot; default "hexagon"
#' @param threshVar CHARACTER STRING | Name of column in dataMetrics object that is used to threshold significance; default "FDR"
#' @param threshVal INTEGER | Maximum value to threshold significance from threshVar object; default 0.05
#' @param geneList CHARACTER ARRAY | List of gene IDs to be drawn onto the scatterplot matrix of all data. Use this parameter if you have predetermined subset of genes to be superimposed. Otherwise, dataMetrics, threshVar, and threshVal will be used to create genes to be superimposed onto the volcano plot; default NULL
#' @param saveFile BOOLEAN [TRUE | FALSE] | Save file to outDir; default TRUE
#' @param outDir CHARACTER STRING | Output directory to save all plots; default current directory
#' @param xbins INTEGER | Number of bins partitioning the range of the plot; default 10
#' @param pointColor CHARACTER STRING | Color of overlaid points on scatterplot matrix; default "orange"
#' @param pointSize INTEGER | Size of plotted points; default 0.5
#' @param logFC CHARACTER STRING | Name of column in dataMetrics object that contains log fold change values; default "logFC"
#' @param PValue CHARACTER STRING | Name of column in dataMetrics object that contains p-values; default "PValue"
#' @param hover BOOLEAN [TRUE | FALSE] | Allow to hover over points to identify IDs; default FALSE
#' 
#' @importFrom dplyr filter select %>%
#' @importFrom GGally ggpairs wrap
#' @importFrom ggplot2 ggplot aes_string aes geom_point xlim ylim geom_hex coord_cartesian xlab ylab geom_ribbon geom_boxplot geom_line geom_abline theme_gray ggtitle labs element_text scale_fill_gradientn
#' @importFrom grDevices jpeg dev.off
#' @importFrom hexbin hexbin hcell2xy
#' @importFrom htmlwidgets onRender
#' @importFrom plotly plotlyOutput ggplotly renderPlotly layout
#' @importFrom shiny verbatimTextOutput fluidPage reactive renderPrint shinyApp
#' @importFrom stats lm predict
#' @importFrom tidyr gather
#' @importFrom utils str
#' @importFrom plyr mapvalues
#' 
#' @export
#' @examples
#' # Example 1: Plot volcano plot with default settings for overlaid points (FDR < 0.05).
#' 
#' data(soybean_ir_sub)
#' data(soybean_ir_sub_metrics)
#' ret <- plotVolcano(soybean_ir_sub, soybean_ir_sub_metrics, pointSize = 1, saveFile = FALSE)
#' ret[[1]]
#' 
#' # Example 2: Plot volcano plot and overlay points with PValue < 1e-15.
#' 
#' ret <- plotVolcano(soybean_ir_sub, soybean_ir_sub_metrics, pointColor = "red", pointSize = 1,
#'   threshVar = "PValue", threshVal = 1e-15, saveFile = FALSE)
#' ret[[1]]
#' 
#' # Example 3: Plot volcano plot and overlay points with PValue < 1e-15. This time, plot all 
#' points (instead of hexagons) for the background.
#' 
#' ret <- plotVolcano(soybean_ir_sub, soybean_ir_sub_metrics, pointColor = "red", pointSize = 1,
#'   threshVar = "PValue", threshVal = 1e-15, option = "allPoints", saveFile = FALSE)
#' ret[[1]]
#' 
#' # Example 4: Plot volcano plot with points in background and overlay points with
#' PValue < 1e-15. This time, use a value of TRUE for the hover parameter so that you can hover 
#' over overlaid points and determine their IDs .
#' 
#' ret <- plotVolcano(soybean_ir_sub, soybean_ir_sub_metrics, pointColor = "red", pointSize = 1,
#'   threshVar = "PValue", threshVal = 1e-15, option = "allPoints", saveFile = FALSE, 
#'   hover = TRUE)
#' ret[[1]]
#' 

plotVolcano = function(data = data, dataMetrics = dataMetrics, threshVar = "FDR", threshVal = 0.05, geneList = NULL, saveFile=TRUE, outDir = getwd(), pointColor = "orange", pointSize = 0.5, xbins = 10, logFC = "logFC", PValue = "PValue", option = "hexagon", hover = FALSE){
  
  countColor2 <- counts <- hexID <- ID <- NULL
  
  colNames <- colnames(data)
  myPairs <- unique(sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1]))
  myPairs <- myPairs[-which(myPairs=="ID")]
  colGroups <- sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1])
  
  ifelse(!dir.exists(outDir), dir.create(outDir), FALSE)
  
  ret = list()
  for (k in 1:(length(myPairs)-1)){
    for (j in (k+1):length(myPairs)){
      group1 = myPairs[k]
      group2 = myPairs[j]
      datSel = cbind(ID=data$ID, data[,which(colGroups %in% c(group1, group2))])
      
      curPairDF1 = dataMetrics[[paste0(group1, "_", group2)]]
      curPairDF2 = dataMetrics[[paste0(group2, "_", group1)]]
      
      if (!is.null(curPairDF1)){
        curPairDF = curPairDF1
      }
      if (!is.null(curPairDF2)){
        curPairDF = curPairDF2
      }
      
      # Set any PValues=0 to be next lowest value
      cpd0 = which(curPairDF[[PValue]]==0)
      curPairDF[[PValue]][cpd0] = sort(unique(curPairDF[[PValue]]))[2]
      
      xMax = max(curPairDF[[logFC]])
      xMin = min(curPairDF[[logFC]])
      yMax = -log(min(curPairDF[[PValue]])) #+1e-10
      yMin = -log(max(curPairDF[[PValue]]))#+1e-10
      fcMax = ceiling(max(exp(xMax), 1/exp(xMin)))
      
      curPairSel = curPairDF[which(curPairDF[[threshVar]] < threshVal),]
      degData = filter(datSel, ID %in% curPairSel$ID)
      if (!is.null(geneList)){
        curPairSel = curPairDF[which(curPairDF$ID %in% geneList),]
        degData = filter(datSel, ID %in% geneList)
      }
      
      x = curPairDF[[logFC]]
      y = -log(curPairDF[[PValue]])#+1e-10
      x2 = curPairSel[[logFC]]
      y2 = -log(curPairSel[[PValue]])#+1e-10
      h = hexbin(x=x, y=y, xbins=xbins, shape=3, IDs=TRUE, xbnds=c(xMin, xMax), ybnds=c(yMin, yMax))
      hexdf = data.frame (hcell2xy (h),  hexID = h@cell, counts = h@count)
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
      
      datSp1 <- sapply(strsplit(levels(hexdf$countColor2), "-"), function(x) x[1])
      datSp2 <- sapply(strsplit(datSp1, "\\+"), function(x) x[1])
      
      bin <- mapvalues(hexdf$countColor2, from = levels(hexdf$countColor2), to = datSp2)
      
      hexdf$bin <- as.integer(bin)
        
      if (option == "allPoints"){
        mainPoints = data.frame(x=x, y=y)
        overlayPoints = data.frame(x=x2, y=y2)
        p <- ggplot(mainPoints, aes(x=x, y=y)) + geom_point(size = pointSize) + geom_point(data = overlayPoints, aes_string(x=x2, y=y2), colour = pointColor, size = pointSize) + theme(axis.text=element_text(size=15), axis.title=element_text(size=15), legend.title=element_text(size=15), legend.text=element_text(size=15)) + coord_cartesian(xlim = c(xMin, xMax), ylim = c(yMin, yMax)) + xlab(logFC) + ylab(paste0("-log10(", PValue, ")"))
        
        if (hover == TRUE){
          IDs=curPairSel$ID
          gP <- ggplotly(p)
          gP[["x"]][["data"]][[1]][["hoverinfo"]] <- "none"
          newText = IDs
          gP[["x"]][["data"]][[2]][["text"]] <- newText
        }
      }
      else{
        overlayPoints = data.frame(x=x2, y=y2)
        p <- ggplot(hexdf, aes(x=x, y=y)) + geom_hex(stat="identity", aes(fill=bin)) + scale_fill_gradientn(colors=rev(clrs[-1]), guide="legend", labels=levels(hexdf$countColor2), name="Count") + theme(axis.text=element_text(size=15), axis.title=element_text(size=15), legend.title=element_text(size=15), legend.text=element_text(size=15)) + coord_cartesian(xlim = c(xMin, xMax), ylim = c(yMin, yMax)) + geom_point(data=overlayPoints, aes_string(x=x2, y=y2), color = pointColor, size = pointSize, inherit.aes = FALSE) + xlab(logFC) + ylab(paste0("-log10(", PValue, ")"))
        
        if (hover == TRUE){
          IDs=curPairSel$ID
          gP <- ggplotly(p)

          for (l in 1:(length(gP[["x"]][["data"]])-1)){
            gP[["x"]][["data"]][[l]][["hoverinfo"]] <- "none"  
          }
          newText = IDs
          gP[["x"]][["data"]][[length(gP[["x"]][["data"]])]][["text"]] <- newText 
        }
      }
   
      if (saveFile == TRUE){
        jpeg(filename=paste0(outDir, "/", group1, "_", group2, "_degVolcano.jpg"), height=700, width=1100)
        print(p)
        dev.off()
      }
      if (hover ==FALSE){
        ret[[paste0(group1, "_", group2)]] <- p     
      }
      else{
        ret[[paste0(group1, "_", group2)]] <- gP        
      }

    }
  }
  invisible(ret)
}
