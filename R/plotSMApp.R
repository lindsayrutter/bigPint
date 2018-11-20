PKGENVIR <- new.env(parent=emptyenv()) # package level envir  

#' @title Plot interactive scatterplot matrices
#' 
#' @description Plot interactive scatterplot matrices.
#' 
#' @param data DATA FRAME | Read counts
#' @param xbins INTEGER | Number of bins partitioning the range of the plot; 
#' default 10 
#' @importFrom plotly plotlyOutput ggplotly renderPlotly config
#' @importFrom ggplot2 ggplot aes_string aes xlim ylim geom_boxplot theme
#' geom_hex geom_abline coord_cartesian
#' @importFrom shiny verbatimTextOutput fluidPage reactive renderPrint shinyUI 
#' sliderInput shinyServer shinyApp HTML br reactiveValues strong em div p img
#' @importFrom htmlwidgets onRender
#' @importFrom utils str
#' @importFrom tidyr gather
#' @importFrom stats qt lm coef
#' @importFrom hexbin hexbin hcell2xy
#' @importFrom stringr str_replace str_trim
#' @importFrom dplyr %>% select
#' @importFrom shinycssloaders withSpinner
#' @importFrom shinydashboard dashboardSidebar sidebarMenu menuItem
#' dashboardBody tabItems box tabItem dashboardPage dashboardHeader
#' @importForm GGally ggpairs
#' @return A Shiny application that shows a scatterplot matrix with hexagon
#' bins and allows users to click on hexagon bins to determine how many genes
#' they each contain. The user can download a file that contains the gene IDs
#' that are located in the clicked hexagon bin.
#' @export
#' @examples
#' # Example: Create interactive scatterplot matrix for first two treatment
#' # groups of data.
#' 
#' data(soybean_cn_sub)
#' soybean_cn_sub <- soybean_cn_sub[,1:7]
#' app <- plotSMApp(data=soybean_cn_sub)
#' if (interactive()) {
#'     shiny::runApp(app)
#' }

plotSMApp = function(data=data, xbins=10){
appDir <- system.file("shiny-examples", "plotSMApp", package = "bigPint")
if (appDir == "") {
    stop("Could not find example directory. Try re-installing `bigPint`.",
    call. = FALSE)
}
PKGENVIR$DATA <- data # put the data into envir
return(appDir)
}