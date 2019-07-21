library(plotly)
library(ggplot2)
library(shiny)
library(htmlwidgets)
library(utils)
library(tidyr)
library(stats)
library(hexbin)
library(stringr)
library(dplyr)
library(shinycssloaders)
library(shinydashboard)
library(shinycssloaders)
library(Hmisc)
library(RColorBrewer)

options(spinner.color.background="#F5F5F5")

data <- bigPint:::PKGENVIR$DATA ## read the data from envir
pointColor <- bigPint:::PKGENVIR$POINTCOLOR ## read the data from envir
pcpDat <- data

sidebar <- shinydashboard::dashboardSidebar(
  width = 180,
  hr(),
  shinydashboard::sidebarMenu(id="tabs",
    shinydashboard::menuItem("Application", tabName="plot1"), #from hexPlot plot1
    shinydashboard::menuItem("About", tabName = "about", selected=TRUE)
  )
)

body <- shinydashboard::dashboardBody(
  shinydashboard::tabItems(
    shinydashboard::tabItem(tabName = "plot1",
      fluidRow(
        column(width = 12,
               shinydashboard::box(width = NULL, shinycssloaders::withSpinner(plotly::plotlyOutput("plot1")), collapsible = FALSE, background = "black", title = "Parallel coordinate plot", status = "primary", solidHeader = TRUE))),
    
    shiny::fluidRow(
      shiny::column(width = 12, shinydashboard::box(width = NULL, downloadButton("downloadData", "Download selected IDs"), DT::dataTableOutput("rectdf"), collapsible = FALSE, title = "Selected genes", status = "primary", solidHeader = TRUE)))),
    
shinydashboard::tabItem(tabName = "about",
      shiny::fluidRow("This application allows users to refine/reduce a set of genes that are represented as parallel coordinate lines. The parallel coordinate lines can often represent a subset of an entire dataset (such as differentially expressed genes). The data we use for the examples below are published RNA-seq data of soybean leaves exposed to iron-rich (group N) and iron-poor (group P) hydroponic soil (Moran Lauter et al., 2016). The two treatments (N and P) each with three replicates. In this example, the subset of data plotted as parallel coordinate lines are the 120 genes that had a FDR values less than 0.01 and a log fold change values less than -4", style='padding:10px;'),
      shiny::fluidRow("As demonstrated in Figure 1, the user can begin to refine/reduce the set of genes by selecting the 'Box Select' button in the Plotly Mode Bar at the top right of the image. In this example, it seems that group P had inconsistent replicates. If, in this case, the user wishes to reduce the set of genes to only include those that have high consistency between the P group replicates, then they can use the 'Box Select' tool to draw a box as is shown in Figure 2. Then, any genes that have a count value for P.1, P.2 or P.3 that is outside of that box will be removed; the result of this is shown in Figure 3.", style='padding:10px;'),
      shiny::fluidRow("The bottom of the application will list the IDs of the remaining genes and provide a button for users to download these IDs. As shown in Figure 4, this example reduced the number of parallel coordinate lines from 120 to 9. Please note that you can download the parallel coordinate plots at any time as static .png files. You need to view this application in a web browser for this function to work. Hover over the top of the interative graphic and the Plotly Mode Bar buttons will appear. After that, simply click on the leftmost button (that has a camera icon) and this will download the static image.", style='padding:10px;'),
      shiny::fluidRow("Go ahead and test this application by switching to the 'Application' tab on the left side of the screen.", style='padding:10px;'),
      br(),
      br(),
      div(p('Figure 1'), style="text-align: center;"),
      div(img(src='Figure1.png', style="width: 75%; height: 75%"), style="text-align: center;"),
      br(),
      br(),
      div(p('Figure 2'), style="text-align: center;"),
      div(img(src='Figure2.png', style="width: 75%; height: 75%"), style="text-align: center;"),
      br(),
      br(),
      div(p('Figure 3'), style="text-align: center;"),
      div(img(src='Figure3.png', style="width: 75%; height: 75%"), style="text-align: center;"), 
      br(),
      br(),
      div(p('Figure 4'), style="text-align: center;"),
      div(img(src='Figure4.png', style="width: 75%; height: 75%"), style="text-align: center;"), 
      br(),
      shiny::fluidRow("1. Moran Lauter, A.N. and Graham, M.A. (2016) NCBI SRA bioproject accession: PRJNA318409, https://www.ncbi.nlm.nih.gov/bioproject/PRJNA318409/.", style='padding:10px;')
    )))


ui <- shinydashboard::dashboardPage(
  shinydashboard::dashboardHeader(title = "Overlaying genes", titleWidth = 180),
  sidebar,
  body
)

server <- function(input, output, session) {
  
  colNms <- colnames(pcpDat[, c(2:(ncol(pcpDat)))])
  nVar <- length(colNms)
  pairs <- unique(sapply(colNms, function(x) unlist(strsplit(x,"[.]"))[1]))
  
  emptyDat = mtcars
  emptyDat$wt[1]=nVar-0.5
  p <- ggplot(emptyDat, aes(x = wt, y = mpg)) + geom_point(alpha=0) + xlim(0,(nVar-1)) +ylim(min(pcpDat[,2:(nVar+1)]),max(pcpDat[,2:(nVar+1)])) + xlab("Sample") + ylab("Count") + scale_x_discrete(limits=colnames(pcpDat[-1]))
  gp <- ggplotly(p)
  gp[["x"]][["data"]][[1]][["text"]] <- NULL
  gp[["x"]][["data"]][[1]][["hoverinfo"]] = "none"
  
  inputRectDf <- reactive({
    req(input$rects)
    df <- input$rects
    return(df)
  })
  
  # Print the selected gene IDs
  output$rectdf = DT::renderDataTable(pcpDat %>% filter(ID %in% inputRectDf()), rownames= FALSE)
  
  session$sendCustomMessage(type = "lines", message=list(pointColor=pointColor))
  
  output$plot1 <- renderPlotly({
    gp %>% onRender("
      function(el, x, data) {

     var rects = [];
      var origPcpDat = data.pcpDat
      var pcpDat = data.pcpDat
      
      var Traces = [];
      var dLength = pcpDat.length
      var vLength = data.nVar
      var cNames = data.colNms
      for (a=0; a<dLength; a++){
      xArr = [];
      yArr = [];
      ids = [];
      for (b=0; b<vLength; b++){
      xArr.push(b+1);
      var keys = Object.values(pcpDat[a])
      ids.push(keys[0])
      yArr.push(pcpDat[a][cNames[b]]);
      }
      var pcpLine = {
      x: xArr,
      y: yArr,
      mode: 'lines',
      line: {
      color: data.pointColor,
      width: 1
      },
      text: ids,
      hoverinfo: 'text',
      opacity: 0.9,
      }
      Traces.push(pcpLine);
      }
      Plotly.addTraces(el.id, Traces);
     
      var selectTime = 1
      
      el.on('plotly_selected', function(e) {
      if (pcpDat.length>0){
      
      var dLength = pcpDat.length
      var selectedPCP = []
      var xMin = e.range.x[0]
      var xMax = e.range.x[1]
      var xMinC = Math.abs(Math.ceil(xMin))
      var xMaxF = Math.floor(xMax)
      var yMin = e.range.y[0]
      var yMax = e.range.y[1]
      var integers = []
      
      if (!((xMax<0) || (xMin>(vLength)))){
      for (a=xMinC; a<(xMaxF+1); a++){
      integers.push(a-1)
      }
      }
      var iLength = integers.length
      
      var selectedPCP = []
      var notSelectedPCP = []
      
      if (selectTime==1){
      if (iLength > 0){
      for (a=0; a<dLength; a++){
      var dat = pcpDat[a]
      var isOut = 0;
      for (b=0; b<iLength; b++){
      var yVal = dat[cNames[integers[b]]]
      if (!(yMin < yVal && yVal < yMax)){
      isOut = 1;
      }
      }
      if (isOut==1){
      selectedPCP.push(a)
      }
      }

      var updateSPCP = []
      var selectedPCPL = selectedPCP.length
      for (a=0; a<selectedPCPL; a++){
      updateSPCP[a]=selectedPCP[a]+1
      }
      
      var update = {
      line: {
      color: 'blue',
      width: 1
      }
      }
      if (selectedPCPL!=0){
      Plotly.deleteTraces(el.id, updateSPCP);
      }
      
      var newDat = []
      var selectedPCPL = selectedPCP.length
      for (a=0; a<dLength; a++){
      var equal = 0;
      for (b=0; b<selectedPCPL; b++){
      if (a==selectedPCP[b]){
      equal=1
      }
      }
      if (equal==0){
      newDat.push(pcpDat[a])
      }
      }
      pcpDat = newDat
      }
      }
      
      if (selectTime>1){
      if (iLength > 0){
      for (a=0; a<dLength; a++){
      var dat = pcpDat[a]
      var isOut = 0;
      for (b=0; b<iLength; b++){
      var yVal = dat[cNames[integers[b]]]
      if (!(yMin < yVal && yVal < yMax)){
      isOut = 1;
      }
      }
      
      if (isOut==0){
      selectedPCP.push(a)
      }
      else{
      notSelectedPCP.push(a)
      }
      }

      var updateNSPCP = []
      var notSelectedPCPL = notSelectedPCP.length
      for (a=0; a<notSelectedPCPL; a++){
      updateNSPCP[a]=notSelectedPCP[a]+1
      }

      if (notSelectedPCPL!=0){
      Plotly.deleteTraces(el.id, updateNSPCP);
      }
      
      var newDat = []
      var selectedPCPL = selectedPCP.length
      for (a=0; a<dLength; a++){
      var equal = 0;
      for (b=0; b<selectedPCPL; b++){
      if (a==selectedPCP[b]){
      equal=1
      }
      }
      if (equal==1){
      newDat.push(pcpDat[a])
      }
      }
      pcpDat = newDat
      }
      
      else{
      for (a=0; a<dLength; a++){
      notSelectedPCP.push(a)
      }
      
      var updateNSPCP = []
      var notSelectedPCPL = notSelectedPCP.length
      for (a=0; a<notSelectedPCPL; a++){
      updateNSPCP[a]=notSelectedPCP[a]+1
      }
      
      var update = {
      line: {
      color: '#FF34B3',
      width: 1
      }
      }
      if (notSelectedPCPL!=0){
      Plotly.deleteTraces(el.id, update, updateNSPCP);
      }
      pcpDat = []
      }
      }
      
      var drawRect = {
      type: 'rect',
      x0: xMin,
      y0: yMin,
      x1: xMax,
      y1: yMax,
      line: {
      color: 'gray',
      width: 1
      },
      fillcolor: 'gray',
      opacity: 0.25
      }
      var update = {
      shapes:rects
      }
      Plotly.relayout(el.id, update)

      selectTime++
      }

for (a=0; a<pcpDat.length; a++){
rects.push(pcpDat[a]['ID']);
}

Shiny.onInputChange('rects', rects);
rects = []
      })
      }", data = list(pcpDat = pcpDat, nVar = nVar, colNms = colNms, pointColor = pointColor))})
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("pcp_", pairs[1], "_", pairs[2], ".csv", sep = "")
    },
    content = function(file) {
      write.table(inputRectDf(), file, row.names = FALSE, col.names = FALSE, quote = FALSE, sep=',')
    }
  )
  
}

shiny::shinyApp(ui = ui, server = server)
