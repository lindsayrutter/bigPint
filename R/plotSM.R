#' @title Plot static scatterplot matrices
#' 
#' @description Plot static scatterplot matrix. Optionally, superimpose differentially expressed genes (DEGs) onto scatterplot matrix.
#' 
##' @details There are seven options:
##' \itemize{
##'  \item{"foldChange": }{Plots DEGs onto a scatterplot matrix of fold changes}
##'  \item{"orthogonal": }{Plots DEGs onto a scatterplot matrix of orthogonal distance}
##'  \item{"hexagon": }{Plot DEGs onto a scatterplot matrix of hexagon binning}
##'  \item{"allPoints": }{Plot DEGs onto a scatterplot matrix of all data points}
##' } 
#' 
#' @param data DATA FRAME | Read counts
#' @param dataMetrics LIST | Differential expression metrics; default NULL
#' @param option CHARACTER STRING ["foldChange" | "orthogonal" | "hexagon" | "allPoints"] | The type of plot; default "allPoints"
#' @param threshVar CHARACTER STRING | Name of column in dataMetrics object that is used to threshold significance; default "FDR"; used in all options
#' @param threshVal INTEGER | Maximum value to threshold significance from threshVar object; default 0.05; used in all options
#' @param geneList CHARACTER ARRAY | List of gene IDs to be drawn onto the scatterplot matrix of all data. Use this parameter if you have predetermined genes to be drawn. Otherwise, use dataMetrics, threshVar, and threshVal to create genes to be drawn onto the scatterplot matrix; default NULL; used in all options
#' @param saveFile BOOLEAN [TRUE | FALSE] | Save file to outDir; default TRUE; used in all options
#' @param outDir CHARACTER STRING | Output directory to save all plots; default current directory; used in all options
#' @param pointColor CHARACTER STRING | Color of overlaid points on scatterplot matrix; default "orange"; used for DEGs in "hexagon" and "allPoints" and used for all points in "foldChange" and "orthogonal"
#' @param pointSize INTEGER | Size of plotted points; default 0.5; used for DEGs in "hexagon" and "allPoints" and used for all points in "foldChange" and "orthogonal"
#' @param xbins INTEGER | Number of bins partitioning the range of the plot; default 10; used in option "hexagon"
#' @param threshFC INTEGER | Threshold of fold change; default 3; used in option "foldChange"
#' @param threshOrth INTEGER | Threshold of orthogonal distance; default 3; used in option "orthogonal"
#' 
#' @importFrom dplyr filter select %>%
#' @importFrom GGally ggpairs wrap
#' @importFrom ggplot2 ggplot aes_string aes geom_point xlim ylim geom_hex coord_cartesian xlab ylab geom_ribbon geom_boxplot geom_line geom_abline theme_gray ggtitle
#' @importFrom grDevices jpeg dev.off
#' @importFrom hexbin hexbin hcell2xy
#' @importFrom htmlwidgets onRender
#' @importFrom plotly plotlyOutput ggplotly renderPlotly layout
#' @importFrom shiny verbatimTextOutput fluidPage reactive renderPrint shinyApp
#' @importFrom stats lm predict
#' @importFrom tidyr gather
#' @importFrom utils str
#' 
#' @export
#' @examples
#' # Read in data and metrics (need for all examples)
#' data(soybean_cn_sub)
#' data(soybean_cn_sub_metrics)
#' data(soybean_ir_sub)
#' data(soybean_ir_sub_metrics)
#' 
#' # Create standardized version of data (need for some examples)
#' library(matrixStats)
#' library(ggplot2)
#' soybean_cn_sub_st <- as.data.frame(t(apply(as.matrix(soybean_cn_sub[,-1]), 1, scale)))
#' soybean_cn_sub_st$ID <- as.character(soybean_cn_sub$ID)
#' soybean_cn_sub_st <- soybean_cn_sub_st[,c(length(soybean_cn_sub_st),
#'   1:length(soybean_cn_sub_st)-1)]
#' colnames(soybean_cn_sub_st) <- colnames(soybean_cn_sub)
#' nID <- which(is.nan(soybean_cn_sub_st[,2]))
#' soybean_cn_sub_st[nID,2:length(soybean_cn_sub_st)] <- 0
#' 
#' \dontrun{
#' # Example 1: Plot scatterplot matrix of points. Saves three plots to outDir because saveFile 
#' equals TRUE by default.
#' 
#' plotSM(soybean_cn_sub, soybean_cn_sub_metrics)
#' }
#' 
#' # Example 2: Plot scatterplot matrix of points. Return list of plots so user can tailor them 
#' (such as add title) and does not save to outDir because saveFile equals FALSE.
#' 
#' ret <- plotSM(soybean_cn_sub, soybean_cn_sub_metrics, pointColor = "pink", saveFile = FALSE)
#' # Determine names of plots in returned list
#' names(ret)
#' ret[["S1_S2"]] + ggtitle("S1 versus S2")
#' ret[["S1_S3"]] + ggtitle("S1 versus S3")
#' ret[["S2_S3"]] + ggtitle("S2 versus S3")
#' 
#' # Example 3: Plot standardized data as scatterplot matrix of points.
#'
#' ret <- plotSM(soybean_cn_sub_st, soybean_cn_sub_metrics, pointColor = "#00C379", saveFile = 
#'   FALSE)
#' ret[[1]] + xlab("Standardized read counts") + ylab("Standardized read counts")
#' 
#' # Example 4: Plot scatterplot matrix of hexagons.
#' 
#' ret <- plotSM(soybean_cn_sub, soybean_cn_sub_metrics, option = "hexagon", xbins = 5, pointSize 
#'   = 0.1, saveFile = FALSE)
#' ret[[2]]
#' 
#' # Example 5: Plot scatterplot matrix of orthogonal distance on the logged data, first without considering the metrics dataset and then considering it.
#' 
#' soybean_ir_sub[,-1] <- log(soybean_ir_sub[,-1] + 1) 
#' ret <- plotSM(soybean_ir_sub, option = "orthogonal", threshOrth = 2.5, pointSize = 0.2, 
#'   saveFile = FALSE)
#' ret[[1]]
#' ret <- plotSM(soybean_ir_sub, soybean_ir_sub_metrics, option = "orthogonal", threshOrth = 2.5, 
#'   pointSize = 0.2, saveFile = FALSE)
#' ret[[1]]
#'
#' # Example 6: Plot scatterplot matrix of fold change.
#' 
#' ret <- plotSM(soybean_cn_sub, soybean_cn_sub_metrics, option = "foldChange", threshFC = 0.5, 
#'   pointSize = 0.2, saveFile = FALSE)
#' ret[[1]]

plotSM = function(data=data, dataMetrics=NULL, option="allPoints", saveFile = TRUE, outDir=getwd(), pointSize=0.5, pointColor = "orange", xbins=10, threshFC=3, threshOrth=3, threshVar="FDR", threshVal=0.05, geneList = NULL){
  
  #### FOLDCHANGE ####
  # Superimpose DEGs onto scatterplot matrix fold change
  if (option=="foldChange"){
    dat <- data
    rownames(dat) <- dat$ID
    rm(data)
    
    lwr <- upr <- ID <- NULL
    
    # Use for all subplots
    minLine = 0
    maxLine = max(dat[,-1])
    inc = (maxLine-minLine)/100
    xv = seq(minLine, maxLine, inc)
    uyv = xv*(threshFC+1)
    lyv = xv/(threshFC+1)
    lineDF = data.frame(xv=xv, uy=uyv, lyv=lyv)
    
    colNames <- colnames(dat)
    myPairs <- unique(sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1]))
    myPairs <- myPairs[-which(myPairs=="ID")]
    colGroups <- sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1])
    
    ifelse(!dir.exists(outDir), dir.create(outDir), FALSE)
    
    k=1
    names_list = list()
    data_list = list()
    for (i in 1:(length(myPairs)-1)){
      for (j in (i+1):length(myPairs)){
        group1 = myPairs[i]
        group2 = myPairs[j]
        datSel <- cbind(ID=dat$ID, dat[,which(colGroups %in% c(group1, group2))])
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
      for (i in 1:length(x)){
        fract = x[i]/y[i]
        if (!is.nan(fract)){
          if(fract > (threshFC + 1) || fract < (1/(threshFC+1))){
            indexPoints = c(indexPoints, i)
          }
        }
      }
      plotPoints = data[indexPoints,]
      
      colNames <- colnames(data)
      myPairs <- unique(sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1]))
      myPairs <- myPairs[-which(myPairs=="ID")]
      colGroups <- sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1])
      group1 = unique(colGroups)[1]
      group2 = unique(colGroups)[2]
      
      if (!is.null(dataMetrics)){
        rowDEG1 <- which(dataMetrics[[paste0(group1,"_",group2)]][threshVar] < threshVal)
        rowDEG2 <- which(dataMetrics[[paste0(group2,"_",group1)]][threshVar] < threshVal)
        rowDEG <- c(rowDEG1, rowDEG2)
        degID1 <- dataMetrics[[paste0(group1,"_",group2)]][rowDEG,]$ID
        degID2 <- dataMetrics[[paste0(group2,"_",group1)]][rowDEG,]$ID
        degID <- c(degID1, degID2)
        degData <- data[which(rownames(data) %in% degID),]
        indexBoth = rownames(plotPoints) %in% rownames(degData)
        indexBlue = rownames(degData) %in% rownames(plotPoints)
        redPoints = plotPoints[indexBoth,]
        greyPoints = plotPoints[!indexBoth,] # problem if indexBoth is integer(0)
        bluePoints = degData[!indexBlue,]
        
        p <- ggplot(lineDF, aes(x=xv, y=lyv)) + geom_line(aes(y = lyv), alpha=0.1) + geom_line(aes(y = uyv), alpha=0.1) + geom_ribbon(aes(ymin=lyv,ymax=uyv), fill="blue", alpha="0.3") + geom_point(data = bluePoints, aes_string(x = xChar, y = yChar), size=pointSize, alpha=0.5, color = "blue") + geom_point(data = greyPoints, aes_string(x=xChar, y = yChar), size=pointSize, color = "darkgrey") + geom_point(data = redPoints, aes_string(x = xChar, y = yChar), size=pointSize, color = "red")
      }
      else{
        p <- ggplot(lineDF, aes(x=xv, y=lyv)) + geom_line(aes(y = lyv), alpha=0.1) + geom_line(aes(y = uyv), alpha=0.1) + geom_ribbon(aes(ymin=lyv,ymax=uyv), fill="blue", alpha="0.3") + geom_point(data = plotPoints, aes_string(x = xChar, y = yChar), size=pointSize)  
      }
      p
    }
    
    my_fn2 <- function(data){
      p <- ggpairs(data %>% select(- ID), lower = list(continuous = my_fn), upper = list(continuous = wrap("cor", size = 4)))
      return(p)
    }
    
    ret <- lapply(data_list, function(x) my_fn2(x))
    names(ret) <- names_list
    
    if (saveFile == TRUE){
      for (i in 1:length(ret)){
        fileName = paste0(outDir, "/", names(ret[i]), "_degSM_", threshFC, "_FC.jpg")
        jpeg(filename=fileName, height=900, width=900)
        print(ret[[i]])
        dev.off()
      } 
    }
    invisible(ret)
  } 
  
  #### HEXAGON ####
  # Superimpose DEGs onto scatterplot matrix of hexagons
  else if (option=="hexagon"){
    
    dat <- data
    rm(data)
    
    counts <- hexID <- ID <- NULL
    
    colNames <- colnames(dat)
    myPairs <- unique(sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1]))
    myPairs <- myPairs[-which(myPairs=="ID")]
    colGroups <- sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1])
    
    ifelse(!dir.exists(outDir), dir.create(outDir), FALSE)
    
    maxVal = max(dat[,-1])
    minVal = min(dat[,-1])
    maxRange = c(minVal, maxVal)
    xbins=xbins
    #buffer = maxRange[2]/xbins
    buffer = (maxRange[2]-maxRange[1])/(xbins/2)
    
    k=1
    names_list = list()
    data_list = list()
    for (i in 1:(length(myPairs)-1)){
      for (j in (i+1):length(myPairs)){
        group1 = myPairs[i]
        group2 = myPairs[j]
        datSel <- cbind(ID=dat$ID, dat[,which(colGroups %in% c(group1, group2))])
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
      
      h <- hexbin(x=x, y=y, xbins=xbins, shape=1, IDs=TRUE, xbnds=maxRange, ybnds=maxRange)
      hexdf <- data.frame (hcell2xy (h),  hexID = h@cell, counts = h@count)
      attr(hexdf, "cID") <- h@cID
      
      p <- ggplot(hexdf, aes(x=x, y=y, fill = counts, hexID=hexID)) + geom_hex(stat="identity") + geom_abline(intercept = 0, color = "red", size = 0.25) + coord_cartesian(xlim = c(maxRange[1]-1*buffer, maxRange[2]+buffer), ylim = c(maxRange[1]-1*buffer, maxRange[2]+buffer))
      
      if (!is.null(geneList)){
        degData <- data[which(rownames(data) %in% geneList),]
        p <- p + geom_point(data = degData, aes_string(x=xChar, y=yChar), inherit.aes = FALSE, color = pointColor, size = pointSize)
      }
      if (!is.null(dataMetrics)){
        colNames <- colnames(data)
        myPairs <- unique(sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1]))
        group1 = myPairs[1]
        group2 = myPairs[2]  
        rowDEG1 <- which(dataMetrics[[paste0(group1,"_",group2)]][threshVar] < threshVal)
        rowDEG2 <- which(dataMetrics[[paste0(group2,"_",group1)]][threshVar] < threshVal)
        rowDEG <- c(rowDEG1, rowDEG2)
        degID1 <- dataMetrics[[paste0(group1,"_",group2)]][rowDEG,]$ID
        degID2 <- dataMetrics[[paste0(group2,"_",group1)]][rowDEG,]$ID
        degID <- c(degID1, degID2)
        dataID = cbind(ID=dat$ID, data)
        degData <- dataID[which(dataID$ID %in% degID),]
        p <- p + geom_point(data = degData, aes_string(x=xChar, y=yChar), inherit.aes = FALSE, color = pointColor, size = pointSize)
      }
      p
    }
    
    my_fn2 <- function(data)
    {
      rownames(data) = data$ID
      p <- ggpairs(data %>% select(- ID), 
                   lower = list(continuous = my_fn), 
                   upper = list(continuous = wrap("cor", size = 4)))
      return(p)
    }
    
    ret <- lapply(data_list, function(x) my_fn2(x))
    names(ret) <- names_list
    
    if (saveFile == TRUE){
      for (i in 1:length(ret)){
        fileName = paste0(outDir, "/", names(ret[i]), "_degSM_Hex_", threshVal, ".jpg")
        jpeg(filename=fileName, height=900, width=900)
        print(ret[[i]])
        dev.off()
      } 
    }
    
    invisible(ret)
  }
  
  #### ORTHOGONAL ####
  # Superimpose DEGs onto scatterplot matrix orthogonal distance
  else if (option=="orthogonal"){
    dat <- data
    rownames(dat) <- dat$ID
    rm(data)
    lwr <- upr <- ID <- NULL
    
    # Use for all subplots
    minLine = 0
    maxLine = max(dat[,-1])
    inc = (maxLine-minLine)/100
    xv = seq(minLine, maxLine, inc)
    uyv = xv+sqrt(2)*threshOrth
    lyv = xv-sqrt(2)*threshOrth
    lineDF = data.frame(xv=xv, uyv=uyv, lyv=lyv)
    
    colNames <- colnames(dat)
    myPairs <- unique(sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1]))
    myPairs <- myPairs[-which(myPairs=="ID")]
    colGroups <- sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1])
    
    ifelse(!dir.exists(outDir), dir.create(outDir), FALSE)
    
    k=1
    names_list = list()
    data_list = list()
    for (i in 1:(length(myPairs)-1)){
      for (j in (i+1):length(myPairs)){
        group1 = myPairs[i]
        group2 = myPairs[j]
        datSel <- cbind(ID=dat$ID, dat[,which(colGroups %in% c(group1, group2))])
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
      for (i in 1:length(x)){
        if(abs(x[i]-y[i]) > sqrt(2)*threshOrth){
          indexPoints = c(indexPoints, i)
        }
      }
      plotPoints = data[indexPoints,]
      
      colNames <- colnames(data)
      myPairs <- unique(sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1]))
      myPairs <- myPairs[-which(myPairs=="ID")] # ADDED
      colGroups <- sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1]) # ADDED
      group1 = unique(colGroups)[1]
      group2 = unique(colGroups)[2]
      
      if (!is.null(dataMetrics)){
        rowDEG1 <- which(dataMetrics[[paste0(group1,"_",group2)]][threshVar] < threshVal)
        rowDEG2 <- which(dataMetrics[[paste0(group2,"_",group1)]][threshVar] < threshVal)
        rowDEG <- c(rowDEG1, rowDEG2)
        degID1 <- dataMetrics[[paste0(group1,"_",group2)]][rowDEG,]$ID
        degID2 <- dataMetrics[[paste0(group2,"_",group1)]][rowDEG,]$ID
        degID <- c(degID1, degID2)
        degData <- data[which(rownames(data) %in% degID),]
        indexBoth = rownames(plotPoints) %in% rownames(degData)
        indexBlue = rownames(degData) %in% rownames(plotPoints)
        redPoints = plotPoints[indexBoth,]
        greyPoints = plotPoints[!indexBoth,] # problem if indexBoth is integer(0)
        bluePoints = degData[!indexBlue,]
        
        p <- ggplot(lineDF, aes(x=xv, y=lyv)) + geom_line(aes(y = lyv), alpha=0.1) + geom_line(aes(y = uyv), alpha=0.1) + geom_ribbon(aes(ymin=lyv,ymax=uyv), fill="blue", alpha="0.3") + geom_point(data = bluePoints, aes_string(x = xChar, y = yChar), size=pointSize, alpha=0.5, color = "blue") + geom_point(data = greyPoints, aes_string(x=xChar, y = yChar), size=pointSize, color = "darkgrey") + geom_point(data = redPoints, aes_string(x = xChar, y = yChar), size=pointSize, color = "red")      
      }
      else{
        p <- ggplot(lineDF, aes(x=xv, y=lyv)) + geom_line(aes(y = lyv), alpha=0.1) + geom_line(aes(y = uyv), alpha=0.1) + geom_ribbon(aes(ymin=lyv,ymax=uyv), fill="blue", alpha="0.3") + geom_point(data = plotPoints, aes_string(x = xChar, y = yChar), size=pointSize)
      }
      p
    }
    
    my_fn2 <- function(data){
      p <- ggpairs(data %>% select(- ID), lower = list(continuous = my_fn))
      return(p)
    }
    
    ret <- lapply(data_list, function(x) my_fn2(x))
    names(ret) <- names_list
    
    if (saveFile == TRUE){
      for (i in 1:length(ret)){
        fileName = paste0(outDir, "/", names(ret[i]), "_degSM_", threshOrth, "_Orth.jpg")
        jpeg(filename=fileName, height=900, width=900)
        print(ret[[i]])
        dev.off()
      } 
    }
    invisible(ret)
  }  
  
  ####ALLPOINTS ####
  # Superimpose DEGs onto scatterplot matrix of points
  else if (option=="allPoints"){
    dat <- data
    rm(data)
    sigGenes <- counts <- hexID <- ID <- NULL
    
    colNames <- colnames(dat)
    myPairs <- unique(sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1]))
    myPairs <- myPairs[-which(myPairs=="ID")]
    colGroups <- sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1])
    
    ifelse(!dir.exists(outDir), dir.create(outDir), FALSE)
    
    maxVal = max(dat[,-1])
    minVal = min(dat[,-1])
    maxRange = c(minVal, maxVal)
    
    # Utility function
    my_fn <- function(data, mapping, ...){
      
      xChar = as.character(mapping$x)[2]
      yChar = as.character(mapping$y)[2]
      x = data[,c(xChar)]
      y = data[,c(yChar)]
      
      p <- ggplot(data, aes_string(x=x, y=y)) + geom_point(size = pointSize) + geom_abline(intercept = 0, color = "red", size = 0.5) + coord_cartesian(xlim = c(maxRange[1], maxRange[2]), ylim = c(maxRange[1], maxRange[2]))
      
      if (!is.null(geneList)){
        degData <- data[which(rownames(data) %in% geneList),]
        p <- p + geom_point(data = degData, aes_string(x=xChar, y=yChar), inherit.aes = FALSE, color = pointColor, size = pointSize)      
      }
      else if (!is.null(dataMetrics)){
        colNames <- colnames(data)
        myPairs <- unique(sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1]))
        group1 = myPairs[1]
        group2 = myPairs[2]
        rowDEG1 <- which(dataMetrics[[paste0(group1,"_",group2)]][threshVar] < threshVal)
        rowDEG2 <- which(dataMetrics[[paste0(group2,"_",group1)]][threshVar] < threshVal)
        rowDEG <- c(rowDEG1, rowDEG2)
        degID1 <- as.character(dataMetrics[[paste0(group1,"_",group2)]][rowDEG,]$ID)
        degID2 <- as.character(dataMetrics[[paste0(group2,"_",group1)]][rowDEG,]$ID)
        degID <- c(degID1, degID2)
        dataID = cbind(ID=dat$ID, data)
        degData <- dataID[which(dataID$ID %in% degID),]
        
        p <- p + geom_point(data = degData, aes_string(x=xChar, y=yChar), inherit.aes = FALSE, color = pointColor, size = pointSize)  
      }
      
      p
    }
    
    k=1
    names_list = list()
    data_list = list()
    for (i in 1:(length(myPairs)-1)){
      for (j in (i+1):length(myPairs)){
        group1 = myPairs[i]
        group2 = myPairs[j]
        datSel <- cbind(ID=dat$ID, dat[,which(colGroups %in% c(group1, group2))])
        data_list[[k]] <- datSel
        names_list[[k]] <- paste0(group1,"_",group2)
        k = k +1
      }
    }
    
    my_fn2 <- function(data)
    {
      rownames(data) = data$ID
      p <- ggpairs(data %>% select(- ID), 
                   lower = list(continuous = my_fn), 
                   upper = list(continuous = wrap("cor", size = 4)))
      return(p)
    }
    
    ret <- lapply(data_list, function(x) my_fn2(x))
    names(ret) <- names_list
    
    if (saveFile == TRUE){
      for (i in 1:length(ret)){
        fileName = paste0(outDir, "/", names(ret[i]), "_degSM_allPoints_", threshVar, "_", threshVal, ".jpg")
        jpeg(filename=fileName, height=900, width=900)
        print(ret[[i]])
        dev.off()
      } 
    }
    
    invisible(ret)
  }
  
  else {
    stop("Check that you selected a valid option parameter")
  }
} 


