PKGENVIR <- new.env(parent=emptyenv()) # package level envir

#' @title Plot interactive volcano plots
#' 
#' @description Plot interactive volcano plots.
#' 
#' @param data DATA FRAME | Read counts
#' @param dataMetrics LIST | Differential expression metrics. This object must 
#' contain one column named "logFC" and one column named "PValue".
#' @param option CHARACTER STRING ["hexagon" | "allPoints"] | The background of 
#' plot; default "hexagon"
#' @param pointColor CHARACTER STRING | Color of overlaid points on scatterplot 
#' matrix; default "orange"
#' @importFrom plotly plotlyOutput ggplotly renderPlotly config
#' @importFrom ggplot2 ggplot aes_string aes xlim ylim geom_boxplot theme
#' @importFrom shiny verbatimTextOutput fluidPage reactive renderPrint shinyUI 
#' sliderInput shinyServer shinyApp HTML br reactiveValues strong em div p img 
#' observeEvent selectInput selectizeInput numericInput actionButton 
#' @importFrom htmlwidgets onRender
#' @importFrom utils str
#' @importFrom tidyr gather
#' @importFrom stats qt lm coef
#' @importFrom hexbin hexbin hcell2xy
#' @importFrom stringr str_replace str_trim
#' @importFrom dplyr %>% select
#' @importFrom shinycssloaders withSpinner
#' @importFrom shinydashboard menuItem tabItem dashboardPage dashboardHeader 
#' dashboardSidebar sidebarMenu tabItems box
#' @importFrom Hmisc cut2
#' @importFrom RColorBrewer brewer.pal
#' @return A Shiny application that shows a volcano plot and allows users to 
#' overlay genes depending on two values, usually a statistical value (such as
#' P-value) and a magnitude change value (such as log fold change). The user
#' can download a file that contains the gene IDs that pass these thresholds.
#' @export
#' @examples
#' # Example 1: Create interactive volcano plot of logged data using hexagon
#' # bins for the background.
#' 
#' data(soybean_cn_sub)
#' data(soybean_cn_sub_metrics)
#' app <- plotVolcanoApp(data = soybean_cn_sub,
#'     dataMetrics = soybean_cn_sub_metrics)
#' if (interactive()) {
#'     shiny::runApp(app)
#' }
#' 
#' # Example 2: Create interactive volcano plot of logged data using points for 
#' # the background.
#' 
#' app <- plotVolcanoApp(data = soybean_cn_sub,
#'     dataMetrics = soybean_cn_sub_metrics, option = "allPoints",
#'     pointColor = "magenta")
#' if (interactive()) {
#'     shiny::runApp(app)
#' }

plotVolcanoApp = function(data=data, dataMetrics=dataMetrics, option="hexagon",
pointColor = "orange"){

helperTestData(data)
if (!is.null(dataMetrics)){
    helperTestDataMetricsVolcanoApp(data, dataMetrics, "PValue", "logFC")
}
  
appDir <- system.file("shiny-examples", "plotVolcanoApp", package = "bigPint")
if (appDir == "") {
    stop("Could not find example directory. Try re-installing `bigPint`.",
    call. = FALSE)
}
PKGENVIR$DATA <- data # put the data into envir
PKGENVIR$DATAMETRICS <- dataMetrics # put the data into envir
PKGENVIR$OPTION <- option # put the option into envir
PKGENVIR$POINTCOLOR <- pointColor # put the pointColor into envir  
return(appDir)
}