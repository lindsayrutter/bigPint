PKGENVIR <- new.env(parent=emptyenv()) # package level envir

#' @title Plot interactive litre plots
#' 
#' @description Plot interactive litre plots.
#' 
#' @param data DATA FRAME | Read counts
#' @param dataMetrics LIST | Differential expression metrics (required)
#' @param option CHARACTER STRING ["hexagon" | "allPoints"] | The background of plot; default "hexagon"; "allPoints" may be too slow depending on data
#' @param geneList CHARACTER ARRAY | List of gene IDs to be drawn onto the litre. Use this parameter if you have predetermined subset of genes to be drawn. Otherwise, all genes in the data object can be superimposed on the litre plot; default NULL
#' @param pointColor CHARACTER STRING | Color of overlaid points on scatterplot matrix; default "orange"
#' @importFrom plotly plotlyOutput ggplotly renderPlotly config
#' @importFrom ggplot2 ggplot aes_string aes xlim ylim geom_boxplot theme
#' @importFrom shiny verbatimTextOutput fluidPage reactive renderPrint shinyUI sliderInput shinyServer shinyApp HTML br reactiveValues strong em div p img observeEvent selectInput selectizeInput numericInput actionButton 
#' @importFrom htmlwidgets onRender
#' @importFrom utils str
#' @importFrom tidyr gather
#' @importFrom stats qt lm coef
#' @importFrom hexbin hexbin hcell2xy
#' @importFrom stringr str_replace str_trim
#' @importFrom dplyr %>% select
#' @importFrom shinycssloaders withSpinner
#' @importFrom shinydashboard menuItem tabItem dashboardPage dashboardHeader dashboardSidebar sidebarMenu tabItems box
#' @importFrom Hmisc cut2
#' @importFrom RColorBrewer brewer.pal
#' 
#' @export
#' @examples
#' \dontrun{
#' # Example 1: Create an interactive litre plot for the logged data using default background of hexagons
#' data(soybean_ir_sub)
#' data(soybean_ir_sub_metrics)
#' soybean_ir_sub_log <- soybean_ir_sub
#' soybean_ir_sub_log[,-1] <- log(soybean_ir_sub[,-1]+1)
#' plotLitreApp(data = soybean_ir_sub_log, dataMetrics = soybean_ir_sub_metrics)
#' 
#' # Example 2: Repeat the same process, only now plot background data as individual points. Note this may be too slow now that all points are drawn in the background.
#' plotLitreApp(data = soybean_ir_sub_log, dataMetrics = soybean_ir_sub_metrics, option = "allPoints", pointColor = "red")
#' }

plotLitreApp = function(data=data, dataMetrics=dataMetrics, geneList = NULL, pointColor = "orange", option = "hexagon"){
  appDir <- system.file("shiny-examples", "plotLitreApp", package = "bigPint")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `bigPint`.", call. = FALSE)
  }
  PKGENVIR$DATA <- data # put the data into envir
  PKGENVIR$DATAMETRICS <- dataMetrics # put the data into envir
  PKGENVIR$GENELIST <- geneList # put the data into envir
  PKGENVIR$POINTCOLOR <- pointColor # put the data into envir
  PKGENVIR$OPTION <- option # put the data into envir
  shiny::runApp(appDir, display.mode = "normal")
}
