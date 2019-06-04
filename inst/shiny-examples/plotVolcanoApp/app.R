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
library(plyr)
library(shinycssloaders)
library(shinydashboard)
library(shinycssloaders)
library(Hmisc)
library(RColorBrewer)

options(spinner.color.background="#F5F5F5")
data <- bigPint:::PKGENVIR$DATA ## read the data from envir
dataMetrics <- bigPint:::PKGENVIR$DATAMETRICS ## read the dataMetrics from envir
option <- bigPint:::PKGENVIR$OPTION ## read the option from envir
pointColor <- bigPint:::PKGENVIR$POINTCOLOR ## read the pointColor from envir

dat <- data

datCol <- colnames(dat)[-which(colnames(dat) %in% "ID")]
myPairs <- unique(sapply(datCol, function(x) unlist(strsplit(x,"[.]"))[1]))
myMetrics <- colnames(dataMetrics[[1]])[-which(colnames(dataMetrics[[1]]) %in% "ID")]
values <- reactiveValues(x=0, selPair=NULL)

sidebar <- shinydashboard::dashboardSidebar(
  width = 180,
  hr(),
  shinydashboard::sidebarMenu(id="tabs",
    shinydashboard::menuItem("Application", tabName="volPlot"),
    shinydashboard::menuItem("About", tabName = "about", selected=TRUE)
  )
)

body <- shinydashboard::dashboardBody(
  shinydashboard::tabItems(
    shinydashboard::tabItem(tabName = "volPlot",
      fluidRow(
        column(width = 4, 
           shinydashboard::box(width = NULL, status = "primary", title = "Plot metrics", solidHeader = TRUE,
         shiny::selectizeInput("selPair", "Treatment pairs:", choices = myPairs, multiple = TRUE, options = list(maxItems = 2)),
         checkboxInput("absLFC", "Absolute log Fold change?", FALSE),
         uiOutput("slider"),
         sliderInput("PValue", "P-value:", min = 0, max = 1, value = 0.05, step = 0.01),
         shiny::numericInput("binSize", "Hexagon size:", value = 10, min = 1),
         shiny::numericInput("pointSize", "Point size:", value = 2, min = 1),
         shiny::actionButton("goButton", "Plot gene subset!"))),
  column(width = 8,
           shinydashboard::box(width = NULL, shinycssloaders::withSpinner(plotly::plotlyOutput("volPlot")), collapsible = FALSE, background = "black", title = "Volcano plot", status = "primary", solidHeader = TRUE))),
      
  
      fluidRow(
        column(width = 12,
               shinydashboard::box(width = NULL, shinycssloaders::withSpinner(plotly::plotlyOutput("boxPlot")), collapsible = FALSE, background = "black", title = "Parallel coordinate plot", status = "primary", solidHeader = TRUE))),
  
   shiny::fluidRow(
     shiny::column(width = 12, shinydashboard::box(width = NULL, downloadButton("downloadData", "Download selected IDs"), br(), br(), shiny::verbatimTextOutput("selectedValues1"), collapsible = TRUE, title = "Selected Gene IDs", status = "primary", solidHeader = TRUE)))),
  
shinydashboard::tabItem(tabName = "about",
      shiny::fluidRow("This application allows users to superimpose a set of genes onto a volcano plot. The volcano plot displays significance versus fold change for each observation in a dataset. It is commonly used in RNA-seq data when scientists wish to quickly identify subsets of genes that undergo meaningful changes. The data we use for the examples below are published RNA-seq data of soybean developmental stages (Brown and Hudson, 2015). They contain three treatments (S1, S2, and S3).", style='padding:10px;'),
      shiny::fluidRow("RNA-seq data tends to be large; for instance, this dataset contains 73,320 genes. Plotting each point would reduce the speed of the functionality as well as cause overplotting problems. As a result, we use hexagon bins to construct the background of the volcano plot as is shown in the right side of Figure 1 shown below.", style='padding:10px;'),
      shiny::fluidRow("This application comes with several input fields as is shown on the left side of Figure 1. The user must choose exactly two treatment groups in the 'Treatment pairs' tab. They can also select a minimum threshold for the log fold change and a maximum threshold for the p-value, which will be used to create the gene subset to be overlaid. Users can also check a box to denote whether the log fold change threshold should also be applied to the absolute value of negative values. Upon making these decisions, the user can then select the 'Plot gene subset!' button to superimpose the subset of genes that were above the designated log fold change threshold and below the designated p-value threshold.", style='padding:10px;'),
      shiny::fluidRow("In Figure 1, we see the user selected treatment pairs S2 and S3, which generated the background of blue hexagon bins. The user then overlaid any genes that had a log fold change greater than 2 and a p-value less than 0.05 as pink points. Figure 2 shows a similar situation, only now the user checked the box for the 'Absolute log fold change'. Hence, genes that had a log fold change greater than 2 or less than -2 and a p-value less than 0.05 are overlaid as pink points.", style='padding:10px;'), 
      shiny::fluidRow("The subset of genes is also superimposed as pink parallel coordinate lines on top of a side-by-side box plot (which represents the full dataset) as shown in Figure 3. Finally, the selected gene IDs are listed at the bottom of the app as demonstrated in Figure 4, and the user can download this list.", style='padding:10px;'),
      shiny::fluidRow("Please note that you can download the volcano plots at any time as static .png files. You need to view this application in a web browser for this function to work. Hover over the top of the interative graphic and the Plotly Mode Bar buttons will appear. After that, simply click on the leftmost button (that has a camera icon) and this will download the static image. Go ahead and test this application by switching to the 'Application' tab on the left side of the screen.", style='padding:10px;'),                       
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
      shiny::fluidRow("1. Brown, A.V. and Hudson, K.A. (2015) Developmental profiling of gene expression in soybean trifoliate leaves and cotyledons.", em(" BMC Plant Biology, "), strong("15"), ", 169.", style='padding:10px;')
    )))

ui <- shinydashboard::dashboardPage(
  shinydashboard::dashboardHeader(title = "Overlaying genes", titleWidth = 180),
  sidebar,
  body
)

server <- function(input, output, session) {
  
  shiny::observeEvent(input$goButton, values$x <- values$x + 1)
  shiny::observeEvent(input$selMetric, values$x <- 0)
  shiny::observeEvent(input$selOrder, values$x <- 0)
  shiny::observeEvent(input$selPair, values$x <- 0)
  shiny::observeEvent(input$binSize, values$x <- 0)
  shiny::observeEvent(input$selPair, values$selPair <- input$selPair)
  
  # Define largest fold change dynamically based on data
  fcInMax <- max(ldply(dataMetrics, rbind)[["logFC"]])
  
  # Construct dynamic input Shiny slider for fold change
  output$slider <- renderUI({
    sliderInput("logFC", "Log fold change:", min=0, max=fcInMax, value=ceiling((fcInMax)/3), step=0.1)
  })
  
  # Create data subset based on two letters user chooses
  datSel <- eventReactive(input$selPair, {
    validate(need(length(input$selPair) == 2, "Select a pair of treatments."))
    sampleIndex <- reactive(which(sapply(colnames(dat), function(x) unlist(strsplit(x,"[.]"))[1]) %in% c(input$selPair[1], input$selPair[2])))
    dat[,c(1, sampleIndex())]
  }, ignoreNULL = FALSE)
  
  # Set the main data frame to include only the pairs selected by user
  curPairDF <- eventReactive(input$selPair, {
    validate(need(length(input$selPair) == 2, "Select a pair of treatments."))
    curPairDF <- dataMetrics[[paste0(input$selPair[1], "_", input$selPair[2])]]
    
    if (is.null(curPairDF)){
      curPairDF <- dataMetrics[[paste0(input$selPair[2], "_", input$selPair[1])]]      
    }
    
    cpd0 = which(curPairDF[["PValue"]]==0)
    curPairDF[["PValue"]][cpd0] = sort(unique(curPairDF[["PValue"]]))[2]
    curPairDF
  })
  
  # Subset other data frame to include only genes with small PValue
  curPairSel <- eventReactive(values$x, {
    validate(need(values$x > 0, "Plot a gene."))
    if (input$absLFC==TRUE){
      curPairSel = curPairDF()[which(curPairDF()[["PValue"]] < input$PValue & abs(curPairDF()[["logFC"]]) > input$logFC),]
    }
    else{
      curPairSel = curPairDF()[which(curPairDF()[["PValue"]] < input$PValue & curPairDF()[["logFC"]] > input$logFC),]
    }
    curPairSel
  })
  
  # Declare shiny output volcano plot
  output$volPlot <- plotly::renderPlotly({
    
    xbins= input$binSize
    xMax = max(curPairDF()[["logFC"]])
    xMin = min(curPairDF()[["logFC"]])
    yMax = -log(min(curPairDF()[["PValue"]]))
    yMin = -log(max(curPairDF()[["PValue"]]))
    fcMax = ceiling(max(exp(xMax), 1/exp(xMin)))
    buffer = (xMax-xMin)/xbins/2
    
    x = curPairDF()[["logFC"]]
    y = -log(curPairDF()[["PValue"]])#+1e-10
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
  
    if (option=="hexagon"){
      p <- reactive(ggplot2::ggplot(hexdf, aes(x=x, y=y, hexID=hexID, counts=counts, fill=countColor2)) + geom_hex(stat="identity") + scale_fill_manual(labels = as.character(my_breaks), values = rev(clrs), name = "Count") + theme(axis.text=element_text(size=15), axis.title=element_text(size=15), legend.title=element_text(size=15), legend.text=element_text(size=15)) + coord_cartesian(xlim = c(xMin-buffer, xMax+buffer), ylim = c(yMin, yMax)) + xlab("logFC") + ylab(paste0("-log10(", "PValue", ")")))
      
      gP <- eventReactive(p(), {
        gP <- plotly::ggplotly(p(), height = 400)
        for (i in 1:(length(gP$x$data)-1)){
          info <- gP$x$data[i][[1]]$text
          info2 <- strsplit(info,"[<br/>]")
          myIndex <- which(startsWith(info2[[1]], "counts:"))
          gP$x$data[i][[1]]$text <- info2[[1]][myIndex]
        }
        gP$x$data[length(gP$x$data)][[1]]$text <- NULL
        gP
      })      
    }
    else{

      mainPoints = data.frame(x=x, y=y)
      p <- reactive(ggplot2::ggplot(mainPoints, aes(x=x, y=y)) + geom_point(size = input$pointSize) + theme(axis.text=element_text(size=15), axis.title=element_text(size=15), legend.title=element_text(size=15), legend.text=element_text(size=15)) + coord_cartesian(xlim = c(xMin, xMax), ylim = c(yMin, yMax)) + xlab("logFC") + ylab(paste0("-log10(", "PValue", ")")))
      
      gP <- eventReactive(p(), {
        gP <- plotly::ggplotly(p(), height = 400)
        gP[["x"]][["data"]][[1]][["hoverinfo"]] <- "none"
        gP
      })
      
    }

    plotlyVol <- reactive(gP())
    
    # Tailor interactivity of the plotly volcano plot object using custom JavaScript
    plotlyVol() %>% onRender("
     function(el, x, data) {
     noPoint = x.data.length;
     Shiny.addCustomMessageHandler('points', function(drawPoints) {
     ids = drawPoints.geneID
     if (x.data.length > noPoint){
     Plotly.deleteTraces(el.id, x.data.length-1);
     }
     var Traces = [];
     var trace = {
     x: drawPoints.geneX,
     y: drawPoints.geneY,
     text: ids,
     hoverinfo: 'text',
     mode: 'markers',
     marker: {
     color: drawPoints.pointColor,
     size: drawPoints.pointSize,
     },
     showlegend: false
     };
     Traces.push(trace);
     Plotly.addTraces(el.id, Traces);
     });}")
  })
  
  # If the user changes the superimposed gene
  observe({
    
    geneX <- curPairSel()[["logFC"]]
    geneY <- -log(curPairSel()[["PValue"]])
    geneID <- curPairSel()[["ID"]]
    pointSize <- input$pointSize * 4
    
    # Send x and y values of selected row into onRender() function
    session$sendCustomMessage(type = "points", message=list(geneX=geneX, geneY=geneY, geneID=geneID, pointSize = pointSize, pointColor = pointColor))
  })

  output$selectedValues1 = renderPrint({
    if ( nrow(curPairSel()) > 50) { 
      cat(paste0("Number of genes: ", nrow(curPairSel()), ". Only listing first 50 genes."))
    }
    else{
      cat(paste("Number of genes:", nrow(curPairSel())))
    }
    cat("\n")
    cat("\n")
    cat(curPairSel()[["ID"]][1:min(nrow(curPairSel()), 50)],sep="\n")
  })
  
  # Declare Shiny output boxplot
  output$boxPlot <- plotly::renderPlotly({
    nVar = reactive(ncol(datSel()))
    colNms <- reactive(colnames(datSel()[, c(2:nVar())]))
    
    boxDat <- eventReactive(datSel(), {
      boxDat <- datSel()[, c(1:nVar())] %>% gather(key, val, -c(ID))
      colnames(boxDat) <- c("ID", "Sample", "Count")
      boxDat
    })
    
    # Create reactive expression of plotly background boxplot
    BP <- reactive(ggplot2::ggplot(boxDat(), aes(x = Sample, y = Count)) + geom_boxplot() + labs(y = "Read count"))
    ggBP <- reactive(plotly::ggplotly(BP(), width=600, height = 400))
    
    # Tailor interactivity of the plotly boxplot object using custom JavaScript
    ggBP() %>% onRender("
      function(el, x, data) {
      Shiny.addCustomMessageHandler('lines', function(drawLines) {
      
      i = x.data.length
      if (i > 1){
      while (i > 1){
      Plotly.deleteTraces(el.id, (i-1));
      i--;
      }
      }
      
      var Traces = [];
      var dLength = drawLines.pcpDat.ID.length
      var dLength2 = drawLines.pcpDat['ID'].length
      var vLength = drawLines.nVar
      var cNames = drawLines.colNms
      
      for (a=0; a<dLength; a++){
      xArr = [];
      yArr = [];
      for (b=0; b<(vLength-1); b++){
      xArr.push(b+1)
      yArr.push(drawLines.pcpDat[cNames[b]][a]);
      }
      var traceLine = {
      x: xArr,
      y: yArr,
      mode: 'lines',
      text: drawLines.geneID,
      hoverinfo: 'text',

      line: {
      color: drawLines.pointColor,
      width: 1.5
      },
      opacity: 0.9,
      }
      Traces.push(traceLine);
      }
      Plotly.addTraces(el.id, Traces);
      });}")
  })
  
  observe({
    validate(need(values$x > 0, "Plot a gene."))
    degData = filter(isolate(datSel()), ID %in% curPairSel()$ID)
    nVar = ncol(degData)
    colNms <- colnames(degData[, c(2:nVar)])
    geneID <- curPairSel()[["ID"]]
    session$sendCustomMessage(type = "lines", message=list(pcpDat=degData, nVar = nVar, colNms = colNms, geneID=geneID, pointColor = pointColor))
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("volcano_", input$selPair[1], "_", input$selPair[2], ".csv", sep = "")
    },
    content = function(file) {
      write.table(curPairSel()[["ID"]], file, row.names = FALSE, col.names = FALSE, quote = FALSE, sep=',')
    }
  )
  }  

shiny::shinyApp(ui = ui, server = server)
