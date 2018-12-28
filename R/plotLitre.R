#' @title Plot static litre plots
#' 
#' @description Plot static litre plots.
#' 
#' @param data DATA FRAME | Read counts
#' @param dataMetrics LIST | Differential expression metrics; default NULL
#' @param geneList CHARACTER ARRAY | List of ID values of genes to be drawn 
#' from data as litre plots. Use this parameter if you have predetermined 
#' genes to be drawn. Otherwise, use dataMetrics, threshVar, and threshVal to 
#' create genes to be drawn; default NULL
#' @param threshVar CHARACTER STRING | Name of column in dataMetrics object 
#' that is used to threshold significance; default "FDR"
#' @param threshVal INTEGER | Maximum value to threshold significance from 
#' threshVar object; default 0.05 
#' @param option CHARACTER STRING ["hexagon" | "allPoints"] | The background 
#' of plot; default "hexagon"
#' @param pointSize INTEGER | Size of plotted points; default 2
#' @param pointColor CHARACTER STRING | Color of gene superimposed on litre 
#' plot; default "orange"
#' @param xbins INTEGER | Number of bins partitioning the range of the plot; 
#' default 10
#' @param outDir CHARACTER STRING | Output directory to save all plots; 
#' default current directory
#' @param saveFile BOOLEAN [TRUE | FALSE] | Save file to outDir; default TRUE
#' @importFrom dplyr filter select %>%
#' @importFrom GGally ggpairs wrap
#' @importFrom ggplot2 ggplot aes_string aes geom_point xlim ylim geom_hex 
#' coord_cartesian xlab ylab geom_ribbon geom_boxplot geom_line geom_abline 
#' theme_gray ggtitle scale_fill_manual coord_fixed labs element_text
#' @importFrom grDevices jpeg dev.off
#' @importFrom hexbin hexbin hcell2xy
#' @importFrom htmlwidgets onRender
#' @importFrom plotly plotlyOutput ggplotly renderPlotly layout
#' @importFrom shiny verbatimTextOutput fluidPage reactive renderPrint 
#' shinyApp
#' @importFrom stats lm predict
#' @importFrom tidyr gather crossing
#' @importFrom utils str
#' @importFrom Hmisc cut2
#' @importFrom RColorBrewer brewer.pal
#' @return List of n elements of litre plots, where n is the number of genes 
#' determined to be superimposed through the dataMetrics or geneList
#' parameter. If the saveFile parameter has a value of TRUE, then each of 
#' these litre plots is saved to the location specified in the outDir 
#' parameter as a JPG file.
#' @export
#' @examples
#' # Example 1: Create litre plots for each of the 61 genes with FDR < 1e-10. 
#' # Examine the first plot (gene "N_P_Glyma.19G168700.Wm82.a2.v1")
#' 
#' data(soybean_ir_sub)
#' soybean_ir_sub[,-1] <- log(soybean_ir_sub[,-1]+1)
#' data(soybean_ir_sub_metrics)
#' ret <- plotLitre(data = soybean_ir_sub,
#'     dataMetrics = soybean_ir_sub_metrics, threshVal = 1e-10,
#'     saveFile = FALSE)
#' length(ret)
#' names(ret)[1]
#' ret[[1]]
#' 
#' # Example 2: Create litre plots for each of the five most significant genes
#' # (low FDR values). View plot for gene "N_P_Glyma.19G168700.Wm82.a2.v1".
#' 
#' geneList = soybean_ir_sub_metrics[["N_P"]][1:5,]$ID
#' ret <- plotLitre(data = soybean_ir_sub, geneList = geneList,
#'     pointColor = "deeppink")
#' names(ret)
#' ret[["N_P_Glyma.19G168700.Wm82.a2.v1"]]
#' 
#' # Example 3: Create one litre plot for each of the five most significant 
#' # genes (low FDR values). View the plot for gene
#' # "N_P_Glyma.19G168700.Wm82.a2.v1". Use points instead of the default 
#' # hexagons as the background.
#' 
#' ret <- plotLitre(data = soybean_ir_sub, geneList = geneList,
#'     pointColor = "deeppink", option = "allPoints")
#' names(ret)
#' ret[["N_P_Glyma.19G168700.Wm82.a2.v1"]]
#' 

plotLitre = function(data=data, dataMetrics=NULL, geneList = NULL, 
    threshVar="FDR", threshVal=0.05, option = "hexagon", pointSize=2,
    pointColor = "orange", xbins=10, outDir=getwd(), saveFile = TRUE){

# Check that input parameters fit required formats
helperTestData(data)
if (is.null(geneList) && !is.null(dataMetrics)){
    helperTestDataMetrics(data, dataMetrics, threshVar)
}

hexID <- counts <- countColor2 <- ID <- NULL
myPairs <- helperMakePairs(data)[["myPairs"]]
colGroups <- helperMakePairs(data)[["colGroups"]]
cols.combn <- combn(myPairs, 2, simplify = FALSE) ### ADDED

ifelse(!dir.exists(outDir), dir.create(outDir), FALSE)

ret <- lapply(cols.combn, function(x){
    group1 = x[1]
    group2 = x[2]
    si1 <- which(colGroups %in% group1)
    si2 <- which(colGroups %in% group2)
    si <- c(si1, si2)
    datSel = data[,c(1, si)]
    
    hexdf = helperMakeHex(datSel, si1, si2, xbins)[["hexdf"]]
    maxRange = helperMakeHex(datSel, si1, si2, xbins)[["maxRange"]]
    clrs = helperMakeHex(datSel, si1, si2, xbins)[["clrs"]]
    my_breaks = helperMakeHex(datSel, si1, si2, xbins)[["my_breaks"]]
    x = helperMakeHex(datSel, si1, si2, xbins)[["x"]]
    y = helperMakeHex(datSel, si1, si2, xbins)[["y"]]
    
    if (option == "hexagon"){
        p <- ggplot(hexdf, aes(x=x, y=y, hexID=hexID, counts=counts,
        fill=countColor2)) + geom_hex(stat="identity") +
        scale_fill_manual(labels = as.character(my_breaks),
        values = rev(clrs), name = "Gene count") +
        geom_abline(intercept = 0, color = "red", size = 0.25) + 
        labs(x = paste0("Read count ", "(", group1, ")"),
        y = paste0("Read count ", "(", group2, ")")) +
        theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15)) +
        coord_fixed(ratio=1)   
    }
    else{
        mainPoints = data.frame(x=x, y=y)
        p <- ggplot(mainPoints, aes(x=x, y=y)) +
        geom_point(size = pointSize) + geom_abline(intercept = 0,
        color = "red", size = 0.25) + labs(x = paste0("Read count ", "(",
        group1, ")"), y = paste0("Read count ", "(", group2, ")")) +
        theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15)) +
        coord_fixed(ratio=1)
    }
    if (is.null(geneList)){
        rowDEG1 <- which(dataMetrics[[paste0(group1,"_",group2)]]
        [threshVar] < threshVal)
        rowDEG2 <- which(dataMetrics[[paste0(group2,"_",group1)]]
        [threshVar] < threshVal)
        geneList <- dataMetrics[[paste0(group1, "_",
        group2)]][c(rowDEG1, rowDEG2),1]
    }

    ret <- lapply(geneList, function(x) {
        currID = x
        currGene = data %>% filter(ID == currID)
        sampleComb = as.data.frame(crossing(as.numeric(currGene[si1]),
        as.numeric(currGene[si2])))
        colnames(sampleComb) = c("x", "y")
        
        ret <- p + geom_point(data = sampleComb, aes(x=x, y=y),
        inherit.aes = FALSE, color = pointColor, size = pointSize) +
        ggtitle(currGene$ID)
        
        if (saveFile == TRUE){
            jpeg(filename=paste0(outDir, "/", group1, "_", group2, "_",
            currGene$ID, "_litre.jpg"), height=700, width=1100)
            print(ret)
            dev.off()
        }
        return(list(plot = ret, name = paste0(group1, "_", group2, "_",
        currGene$ID)))
    })
})
ret <- ret[[1]]
retPlots <- lapply(ret, function(x) {x$plot})
retNames <- lapply(ret, function(x) {x$name})
names(retPlots) <- retNames
invisible(retPlots)
}