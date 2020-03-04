PKGENVIR <- new.env(parent=emptyenv()) # package level envir

#' @title Plot interactive litre plots
#' 
#' @description Plot interactive litre plots.
#' 
#' @param data DATA FRAME | Read counts
#' @param dataMetrics LIST | Differential expression metrics
#' @param dataSE SUMMARIZEDEXPERIMENT | Summarized experiment format that
#' can be used in lieu of data; default NULL
#' @param geneList CHARACTER ARRAY | List of gene IDs to be drawn onto the
#' litre. Use this parameter if you have predetermined subset of genes to be
#' drawn. Otherwise, all genes in the data object can be superimposed on the
#' litre plot; default NULL
#' @param pointColor CHARACTER STRING | Color of overlaid points on scatterplot 
#' matrix; default "orange"
#' @param option CHARACTER STRING ["hexagon" | "allPoints"] | The background
#' of plot; default "hexagon"; "allPoints" may be too slow depending on data
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
#' @importFrom stats setNames
#' @return A Shiny application that shows a litre plot background and allows
#' users to superimpose the subset of genes determined to be superimposed
#' through the dataMetrics or geneList parameter. The application allows users
#' to order how to sequentially superimpose the genes by columns in the
#' dataMetrics parameter. 
#' @export
#' @examples
#' # The first pair of examples use data and dataMetrics
#' # objects as input. The last pair of examples create the same plots now
#' # using the SummarizedExperiment (i.e. dataSE) object input.
#' 
#' # Example 1: Create an interactive litre plot for the logged data using
#' # default background of hexagons.
#' 
#' data(soybean_ir_sub)
#' data(soybean_ir_sub_metrics)
#' soybean_ir_sub_log <- soybean_ir_sub
#' soybean_ir_sub_log[,-1] <- log(soybean_ir_sub[,-1]+1)
#' app <- plotLitreApp(data = soybean_ir_sub_log,
#'     dataMetrics = soybean_ir_sub_metrics)
#' if (interactive()) {
#'     shiny::runApp(app, port = 1234, launch.browser = TRUE)
#' }
#' 
#' # Example 2: Repeat the same process, only now plot background data as 
#' # individual points. Note this may be too slow now that all points are drawn
#' # in the background.
#' 
#' app <- plotLitreApp(data = soybean_ir_sub_log,
#'     dataMetrics = soybean_ir_sub_metrics, option = "allPoints",
#'     pointColor = "red")
#' if (interactive()) {
#'     shiny::runApp(app)
#' }
#' 
#' # Below are the same pair of examples, only now using the
#' # SummarizedExperiment (i.e. dataSE) object as input.
#' 
#' # Example 1: Create an interactive litre plot for the logged data using
#' # default background of hexagons.
#' 
#' \dontrun{
#' data(se_soybean_ir_sub)
#' se_soybean_ir_sub_log <- se_soybean_ir_sub
#' assay(se_soybean_ir_sub_log) <-
#'    log(as.data.frame(assay(se_soybean_ir_sub_log))+1)
#' app <- plotLitreApp(dataSE = se_soybean_ir_sub_log)
#' if (interactive()) {
#'     shiny::runApp(app, port = 1234, launch.browser = TRUE)
#' }
#' }
#' 
#' # Example 2: Repeat the same process, only now plot background data as 
#' # individual points. Note this may be too slow now that all points are
#' # drawn in the background.
#' 
#' \dontrun{
#' app <- plotLitreApp(dataSE = se_soybean_ir_sub_log, option = "allPoints",
#'     pointColor = "red")
#' if (interactive()) {
#'     shiny::runApp(app)
#' }
#' }
#' 

plotLitreApp = function(data=data, dataMetrics=dataMetrics, dataSE=NULL,
geneList = NULL, pointColor = "orange", option = c("hexagon", "allPoints")){

option <- match.arg(option)

if (is.null(dataSE) && is.null(data)){
    helperTestHaveData()
}

if (!is.null(dataSE)){
    #Reverse engineer data
    data <- helperGetData(dataSE)
    
    if (ncol(rowData(dataSE))>0){
        #Reverse engineer dataMetrics
        reDataMetrics <- as.data.frame(rowData(dataSE))
        dataMetrics <- lapply(split.default(reDataMetrics[-1], 
        sub("\\..*", "",names(reDataMetrics[-1]))), function(x)
        cbind(reDataMetrics[1], setNames(x, sub(".*\\.", "", names(x)))))
        for (k in seq_len(length(dataMetrics))){
            colnames(dataMetrics[[k]])[1] = "ID"   
        }
    }
}

if (is.null(dataMetrics)){
    helperTestHaveDataMetrics()
}

helperTestData(data)
if (is.null(geneList) && !is.null(dataMetrics)){
    helperTestDataMetricsLitreApp(data, dataMetrics)
}

appDir <- system.file("shiny-examples", "plotLitreApp", package = "bigPint")
if (appDir == "") {
    stop("Could not find example directory. Try re-installing `bigPint`.",
    call. = FALSE)
}
PKGENVIR$DATA <- data # put the data into envir
PKGENVIR$DATAMETRICS <- dataMetrics # put the data into envir
PKGENVIR$GENELIST <- geneList # put the data into envir
PKGENVIR$POINTCOLOR <- pointColor # put the data into envir
PKGENVIR$OPTION <- option # put the data into envir
return(appDir)
}
