library(shinydashboard)
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
library(data.table)
library(GGally)
library(RColorBrewer)
library(Hmisc)
library(shinycssloaders)

options(warn=-1)

# Read data from envir
data <- bigPint:::PKGENVIR$DATA

# Create new variables based on values read in previously
datCol <- colnames(data)[-which(colnames(data) %in% "ID")]
myPairs <- unique(sapply(datCol, function(x) unlist(strsplit(x,"[.]"))[1]))
values <- reactiveValues(x=0, selPair=NULL, selMetric=NULL, selOrder=NULL)

# Initiate sidebar of Shiny dashboard
sidebar <- shinydashboard::dashboardSidebar(
  width = 180,
  shiny::hr(),
  shinydashboard::sidebarMenu(id="tabs",
    shinydashboard::menuItem("Application", tabName="scatMatPlot"), #hexPlot
    shinydashboard::menuItem("About", tabName = "about", selected=TRUE) #boxPlot
  )
)

# Initiate main body of Shiny dashboard, including Shiny input fields and application description page
body <- shinydashboard::dashboardBody(
  shinydashboard::tabItems(
    shinydashboard::tabItem(tabName = "scatMatPlot",
        shiny::fluidRow(
          shiny::column(width = 12, shinydashboard::box(width = 660, height = 660, withSpinner(plotly::plotlyOutput("scatMatPlot")), collapsible = FALSE, background = "black", title = "Binned scatterplot", status = "primary", solidHeader = TRUE))),
        
        shiny::fluidRow(
          shiny::column(width = 12, shinydashboard::box(width = NULL, withSpinner(plotly::plotlyOutput("boxPlot")), collapsible = FALSE, background = "black", title = "Boxplot", status = "primary", solidHeader = TRUE))),
        
        shiny::fluidRow(
          shiny::column(width = 12, shinydashboard::box(width = NULL, downloadButton("downloadData", "Download selected IDs"), DT::dataTableOutput("selectedValues"), collapsible = FALSE, title = "Selected genes", status = "primary", solidHeader = TRUE)))),

shinydashboard::tabItem(tabName = "about",
        shiny::fluidRow("This application allows you to examine the relationship between all variables in your dataset with an interactive scatterplot matrix. Plotting an individual point for each gene can obscure the number of genes in a given area due to overplotting. As a result, we use hexagon bins in the scatterplot matrix. If you hover over a given hexagon bin of interest, you can determine the number of genes in its area, as shown in Figure 1 below.", style='padding:10px;'),
        shiny::fluidRow("You can also click on a given hexagon bin of interest to overlay the genes it contains across all scatterplots as orange points (Figure 2). Doing so will automatically overlay these same genes as orange parallel coordinate lines across a side-by-side boxplot (which represents your full dataset) immediately below (Figure 3). Moreover, beneath that, you will see an output of the IDs of theses selected genes, which you can optionally download (Figure 4).", style='padding:10px;'),
        shiny::fluidRow("The four figures below were created in simulated data drawn from the normal distribution for didactic purposes.", style='padding:10px;'),
        shiny::fluidRow("Note that you can download the scatterplot matrix or the parallel coordinate plot at any time as static .png files. You need to view this application in a web browser for this function to work. Hover over the top of the interative graphic and the Plotly Mode Bar buttons will appear. After that, simply click on the leftmost button (that has a camera icon) and this will download the static image.", style='padding:10px;'),
        shiny::fluidRow("Please use the 'Application' tab on the left to interact with your data.", strong("When you hover over a hexagon bin, if it appears to contain thousands of genes, please avoid clicking on it."), "Plotting such a large number of genes could crash the program. In any case, hexagon bins that contain that many genes are unlikely to be of interest anyway because they usually occur close to the", em("x=y"), "line. Instead, we are typically interested in hexagon bins that deviate from the" ,em("x=y"), "line, and these bins tend to only contain a handful of genes in the first place.", style='padding:10px;'),
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
        div(img(src='Figure4.png', style="width: 75%; height: 75%"), style="text-align: center;")
    ))
)

# Combine sidebar and main body of Shiny into ui of Shiny application
ui <- shinydashboard::dashboardPage(
  shinydashboard::dashboardHeader(title = "Overlaying genes", titleWidth = 180),
  sidebar,
  body
)

# Inititate server of Shiny application
server <- function(input, output, session) {
  
  output$scatMatPlot <- renderPlotly({
    maxVal = max(abs(data[,-1]))
    maxRange = c(-1*maxVal, maxVal)
    
    # Draw hexagon geoms and red x=y line in scatterplots in lower left corner of matrix
    my_fn <- function(data, mapping, ...){
      xChar = as.character(mapping$x)[2]
      yChar = as.character(mapping$y)[2]
      x = data[,c(xChar)]
      y = data[,c(yChar)]
      h <- hexbin(x=x, y=y, xbins=15, shape=1, IDs=TRUE, xbnds=maxRange, ybnds=maxRange)
      hexdf <- data.frame (hcell2xy (h),  hexID = h@cell, counts = h@count)
      attr(hexdf, "cID") <- h@cID
      p <- ggplot(hexdf, aes(x=x, y=y, fill = counts, hexID=hexID)) + geom_hex(stat="identity") + geom_abline(intercept = 0, color = "red", size = 0.25) + coord_cartesian(xlim = c(-0.5, maxRange[2]+0.5), ylim = c(-0.5, maxRange[2]+0.5)) + theme(legend.position="none")
      p
    }
    
    # Create static scatterplot matrix
    p <- ggpairs(data[,-1], lower = list(continuous = my_fn)) + theme_gray() + theme(legend.position="none")
    # Keep two separate copies
    pS <- p
    
    # Render scatterplot matrix interactive as a plotly object
    ggPS <- ggplotly(pS, width=700, height=600)
    
    # Specify that hovering over hexagons should indicate count number
    myLength <- length(ggPS[["x"]][["data"]])
    for (i in 1:myLength){
      item =ggPS[["x"]][["data"]][[i]]$text[1]
      if (!is.null(item)){
        if (!startsWith(item, "co")){
          ggPS[["x"]][["data"]][[i]]$hoverinfo <- "none"
        }}
      hexHover = ggPS[["x"]][["data"]][[i]]$text
      if (!is.null(hexHover) && grepl("hexID", hexHover)){
        ggPS[["x"]][["data"]][[i]]$text <- strsplit(hexHover, "<")[[1]][1]
        ggPS[["x"]][["data"]][[i]]$t2 <- hexHover
        ggPS[["x"]][["data"]][[i]]$hoverinfo <- "text"
      }
    }
    
    # Specify attributes in hexagons are counts
    for(i in 2:(p$nrow)) {
      for(j in 1:(p$nrow-1)) {
        data[[paste(i,j,sep="-")]] <- attr(pS[i,j]$data, "cID")
      }
    }
    
    # Use this function to supplement the widget's built-in JavaScript rendering logic
    # with additional custom JavaScript code, just for this specific widget object.
    # Usage: onRender(x, jsCode, data = NULL)
    # x - An HTML Widget object
    # jsCode - Character vector containing JavaScript code (see Details)
    # data - An additional argument to pass to the jsCode function. This can be any R
    # object that can be serialized to JSON. If you have multiple objects to pass to the
    # function, use a named list. This is the JavaScript equivalent of the R object passed
    # into onRender as the data argument; this is an easy way to transfer e.g. data frames     
    # without having to manually do the JSON encoding. In this case, data is what the user
    # passes in as the data frame of read counts.
    # Use custom JavaScript to tailor interactivity of the plotly scatterplot matrix object
    ggPSR <- ggPS %>% onRender("
     function(el, x, data) {
     
     function range(start, stop, step){
     var a=[start], b=start;
     while(b<stop){b+=step;a.push(b)}
     return a;
     };
     
     len = Math.sqrt(document.getElementsByClassName('cartesianlayer')[0].childNodes.length);
     AxisNames = [];
     for (i = 1; i < (len+1); i++) {
     AxisNames.push(Object.keys(data[0])[i]);
     }
     noPoint = x.data.length;
     
     // If user clicks on plotly scatterplot matrix object
     el.on('plotly_click', function(e) {

     // If present, delete any old superimposed plotly geoms (orange dots)
     if (x.data.length > noPoint){
     Plotly.deleteTraces(el.id, range(noPoint, (noPoint+(len*(len-1)/2-1)), 1));
     }
     
     // Determine the IDs of genes inside selected hexagon
     xVar = (e.points[0].xaxis._id).replace(/[^0-9]/g,'')
     if (xVar.length == 0) xVar = 1
     yVar = (e.points[0].yaxis._id).replace(/[^0-9]/g,'')
     if (yVar.length == 0) yVar = 1
     myX = len + 1 - (yVar - len * (xVar - 1))
     myY = xVar
     cN = e.points[0].curveNumber
     split1 = (x.data[cN].text).split(' ')
     hexID = (x.data[cN].t2[0]).split(':')[2]
     counts = split1[1].split('<')[0]
     var selRows = [];
     data.forEach(function(row){
     if(row[myX+'-'+myY]==hexID) selRows.push(row);
     });
     selID = []
     for (a=0; a<selRows.length; a++){
     selID.push(selRows[a]['ID'])
     }

     // Save selected row IDs for parallel coordinate plot
     // 'Shiny' below refers to a JavaScript object that is provided by Shiny and is
     // available in JavaScript during the lifetime of an app. Instead of sending a message
     // from Shiny to JavaScript, we can also send messages from JavaScript to Shiny. 
     // Object with name selID, and subsequently use it to send a message back to Shiny. Here,
     // we tell it to make the message available in the R world under the name selID. That is,
     // in R we can now listen for events via input$selID
     Shiny.onInputChange('selID', selID);
     
     // Create traces for selected gene IDs as orange points that state gene names upon hovering
     var Traces = [];
     var i=0;
     var k=1;
     while ((i*len+k)<=Math.pow((len-1),2)) {
     var xArr = [];
     for (a=0; a<selRows.length; a++){
     xArr.push(selRows[a][AxisNames[i]])
     }
     while ((i+k)<len){
     var yArr = [];
     for (a=0; a<selRows.length; a++){
     yArr.push(selRows[a][AxisNames[(len-k)]])
     }
     var trace = {
     x: xArr,
     y: yArr,
     mode: 'markers',
     marker: {
     color: 'orange',
     size: 6
     },
     xaxis: 'x' + (i+1),
     yaxis: 'y' + (i*len+k),
     text: selID,
     hoverinfo: 'text',
     };
     Traces.push(trace);
     k++;
     }
     i++;
     k=1;
     }
     // Push traces to be superimposed onto the plotly scatterplot matrix object
     Plotly.addTraces(el.id, Traces);
     })}
     // Pass the R data object into the JavaScript function
     ", data = data)
    ggPSR
  })
  
  # Read selected gene IDs into Shiny
  selID <- reactive(input$selID)
  
  # Create data subset (read counts) for only the selected gene IDs
  pcpDat <- reactive(data[which(data$ID %in% selID()), ])
  
  # Print the selected gene IDs
  output$selectedValues = DT::renderDataTable(pcpDat(), rownames= FALSE)
  
  # Create static box plot of the full dataset
  colNms <- colnames(data[, -1])
  nVar <- ncol(data)
  boxDat <- data %>% gather(Sample, Count, -c(ID))
  BP <- ggplot(boxDat, aes(x = Sample, y = Count)) + geom_boxplot() + theme_gray()
  
  # Render box plot interactive as a plotly object
  ggBP <- ggplotly(BP, width=700)
  
  output$boxPlot <- renderPlotly({
    ggBP %>% onRender("
      function(el, x, data) {
      
      var Traces = [];
      
      var dLength = data.pcpDat.length
      var vLength = data.nVar
      var cNames = data.colNms
 
      for (a=0; a<dLength; a++){
      xArr = [];
      yArr = [];
      selID = [];
      for (b=0; b<vLength; b++){
      xArr.push(b+1)
      yArr.push(data.pcpDat[a][cNames[b]]);
      selID.push(data.pcpDat[a]['ID']);
      }

      var traceHiLine = {
      x: xArr,
      y: yArr,
      mode: 'lines',
      line: {
      color: 'orange',
      width: 1.5
      },
      opacity: 0.9,
      text: selID,
      hoverinfo: 'text',
      }
      Traces.push(traceHiLine);
      }
      Plotly.addTraces(el.id, Traces);
      // Pass the R objects into the JavaScript function
      }", data = list(pcpDat = pcpDat(), nVar = nVar, colNms = colNms))})
  
  # Create download button
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("SM.csv")
    },
    content = function(file) {
      write.table(pcpDat()[["ID"]], file, row.names = FALSE, col.names = FALSE, quote = FALSE, sep=',')
    }
  )
  
}

shiny::shinyApp(ui = ui, server = server)