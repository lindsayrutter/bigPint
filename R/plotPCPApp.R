PKGENVIR <- new.env(parent=emptyenv()) # package level envir 

#' @title Plot interactive parallel coordinate plots
#' 
#' @description Plot interactive parallel coordinate plots.
#' 
#' @param data DATA FRAME | Read counts for parallel coordinate lines
#' @param pointColor CHARACTER STRING | Color of overlaid points on scatterplot matrix; default "orange"
#' @importFrom plotly plotlyOutput ggplotly renderPlotly layout
#' @importFrom ggplot2 ggplot aes_string geom_point xlim ylim scale_x_discrete
#' @importFrom shiny verbatimTextOutput fluidPage reactive renderPrint shinyApp bootstrapPage basicPage req
#' @importFrom htmlwidgets onRender
#' @importFrom utils str
#' @importFrom dplyr select %>% one_of
#' @export
#' @examples
#' \dontrun{
#' # Example: Create interactive parallel coordinate plot for genes that have
#' # FDR < 0.01 and logFC < -4.
#' 
#' data(soybean_ir_sub)
#' data(soybean_ir_sub_metrics)
#' 
#' # Create standardized version of data
#' library(matrixStats)
#' soybean_ir_sub_st = as.data.frame(t(apply(as.matrix(soybean_ir_sub[,-1]), 1,
#'   scale)))
#' soybean_ir_sub_st$ID = as.character(soybean_ir_sub$ID)
#' soybean_ir_sub_st = soybean_ir_sub_st[,c(length(soybean_ir_sub_st),
#'   1:length(soybean_ir_sub_st)-1)]
#' colnames(soybean_ir_sub_st) = colnames(soybean_ir_sub)
#' nID = which(is.nan(soybean_ir_sub_st[,2]))
#' soybean_ir_sub_st[nID,2:length(soybean_ir_sub_st)] = 0
#' 
#' library(dplyr, warn.conflicts = FALSE)
#' plotGenes = filter(soybean_ir_sub_metrics[["N_P"]], FDR < 0.01, logFC < -4) %>% 
#'   select(ID)
#' pcpDat = filter(soybean_ir_sub_st, ID %in% plotGenes[,1])
#' plotPCPApp(data = pcpDat, pointColor = "purple")
#' }

plotPCPApp = function(data = data, pointColor = "orange"){
  appDir <- system.file("shiny-examples", "plotPCPApp", package = "bigPint")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `bigPint`.", call. = FALSE)
  }
  PKGENVIR$DATA <- data # put the data into envir
  PKGENVIR$POINTCOLOR <- pointColor # put the pointColor into envir
  shiny::runApp(appDir, display.mode = "normal")
}