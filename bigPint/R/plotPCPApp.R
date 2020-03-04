PKGENVIR <- new.env(parent=emptyenv()) # package level envir 

#' @title Plot interactive parallel coordinate plots
#' 
#' @description Plot interactive parallel coordinate plots.
#' 
#' @param data DATA FRAME | Read counts for parallel coordinate lines
#' @param dataSE SUMMARIZEDEXPERIMENT | Summarized experiment format that
#' can be used in lieu of data; default NULL
#' @param pointColor CHARACTER STRING | Color of overlaid points on scatterplot 
#' matrix; default "orange"
#' @importFrom plotly plotlyOutput ggplotly renderPlotly layout
#' @importFrom ggplot2 ggplot aes_string geom_point xlim ylim scale_x_discrete
#' @importFrom shiny verbatimTextOutput fluidPage reactive renderPrint shinyApp 
#' bootstrapPage basicPage req
#' @importFrom htmlwidgets onRender
#' @importFrom utils str
#' @importFrom dplyr select %>% one_of
#' @return A Shiny application that shows a parallel coordinate plot and allows 
#' users to draw rectangular areas across samples and remove genes that are not 
#' inside these areas. The user can download a file that contains the gene IDs 
#' that remain.
#' @export
#' @examples
#' # The first example uses data and dataMetrics objects as
#' # input. The last example creates the same plot now using the
#' # SummarizedExperiment (i.e. dataSE) object input.
#' 
#' # Example: Create interactive parallel coordinate plot for genes that have
#' # FDR < 0.01 and logFC < -4. Standardize genes to have an average of zero
#' # and a standard deviation of one.
#' 
#' data(soybean_ir_sub)
#' data(soybean_ir_sub_metrics)
#' 
#' # Create standardized version of data
#' library(matrixStats)
#' soybean_ir_sub_st = as.data.frame(t(apply(as.matrix(soybean_ir_sub[,-1]), 1,
#'     scale)))
#' soybean_ir_sub_st$ID = as.character(soybean_ir_sub$ID)
#' soybean_ir_sub_st = soybean_ir_sub_st[,c(length(soybean_ir_sub_st),
#'     1:length(soybean_ir_sub_st)-1)]
#' colnames(soybean_ir_sub_st) = colnames(soybean_ir_sub)
#' nID = which(is.nan(soybean_ir_sub_st[,2]))
#' soybean_ir_sub_st[nID,2:length(soybean_ir_sub_st)] = 0
#' 
#' library(dplyr, warn.conflicts = FALSE)
#' plotGenes = filter(soybean_ir_sub_metrics[["N_P"]], FDR < 0.01,
#'     logFC < -4) %>% select(ID)
#' pcpDat = filter(soybean_ir_sub_st, ID %in% plotGenes[,1])
#' 
#' app <- plotPCPApp(data = pcpDat, pointColor = "purple")
#' if (interactive()) {
#'     shiny::runApp(app, display.mode = "normal")
#' }
#' 
#' # Below is the same example, only now using the
#' # SummarizedExperiment (i.e. dataSE) object as input.
#' 
#' # Example: Create interactive parallel coordinate plot for genes that have
#' # FDR < 0.01 and logFC < -4. Standardize genes to have an average of zero
#' # and a standard deviation of one.
#' 
#' \dontrun{
#' data(se_soybean_ir_sub)
#' 
#' # Create standardized version of data
#' library(matrixStats)
#' se_soybean_ir_sub_st = se_soybean_ir_sub
#' assay(se_soybean_ir_sub_st) <-as.data.frame(t(apply(as.matrix(as.data.frame(
#'     assay(se_soybean_ir_sub))), 1, scale)))
#' nID <- which(is.nan(as.data.frame(assay(se_soybean_ir_sub_st))[,1]))
#' assay(se_soybean_ir_sub_st)[nID,] <- 0
#' 
#' # To subset our SummarizedExperiment data by a list of genes, we can
#' # invoke the convertSESubsetGenes() method.
#' 
#' library(dplyr, warn.conflicts = FALSE)
#' geneList <- as.data.frame(rowData(se_soybean_ir_sub_st)) %>%
#'     filter(N_P.FDR <= 0.01) %>% filter(N_P.logFC < -4)
#' geneList <- geneList[,1]
#' pcpDat <- convertSESubsetGenes(se_soybean_ir_sub_st, geneList)
#' 
#' app <- plotPCPApp(dataSE = pcpDat, pointColor = "purple")
#' if (interactive()) {
#'     shiny::runApp(app, display.mode = "normal")
#' }
#' }
#' 

plotPCPApp = function(data = data, dataSE=NULL, pointColor = "orange"){

if (is.null(dataSE) && is.null(data)){
    helperTestHaveData()
}

if (!is.null(dataSE)){
    #Reverse engineer data
    data <- helperGetData(dataSE)
}
    
helperTestData(data)
appDir <- system.file("shiny-examples", "plotPCPApp", package = "bigPint")
if (appDir == "") {
    stop("Could not find example directory. Try re-installing `bigPint`.",
    call. = FALSE)
}
PKGENVIR$DATA <- data # put the data into envir
PKGENVIR$POINTCOLOR <- pointColor # put the pointColor into envir
return(appDir)
}