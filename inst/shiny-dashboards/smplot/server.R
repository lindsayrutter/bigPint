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

load("soybean_cn_sub.rda")
soybean_cn_sub <- soybean_cn_sub[,1:7]
data <- soybean_cn_sub
datCol <- colnames(data)[-which(colnames(data) %in% "ID")]
myPairs <- unique(sapply(datCol, function(x) unlist(strsplit(x,"[.]"))[1]))
values <- reactiveValues(x=0, selPair=NULL, selMetric=NULL, selOrder=NULL)

server <- function(input, output, session) {
  
  output$scatMatPlot <- renderPlotly({
    
    ################################ Prepare scatterplot matrix
    ###########################################################
    
    maxVal = max(abs(data[,-1]))
    maxRange = c(-1*maxVal, maxVal)
    
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
    
    p <- ggpairs(data[,-1], lower = list(continuous = my_fn))
    pS <- p
    
    ggPS <- ggplotly(pS, width=700, height=600)
    
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
    
    for(i in 2:(p$nrow)) {
      for(j in 1:(p$nrow-1)) {
        data[[paste(i,j,sep="-")]] <- attr(pS[i,j]$data, "cID")
      }
    }
    
    ggPS2 <- ggPS %>% onRender("
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
     
     el.on('plotly_click', function(e) {
     
     if (x.data.length > noPoint){
     Plotly.deleteTraces(el.id, range(noPoint, (noPoint+(len*(len-1)/2-1)), 1));
     }
     
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
     // Save selected row IDs for PCP
     Shiny.onInputChange('selID', selID);
     
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
     Plotly.addTraces(el.id, Traces);
     })}
     ", data = data)
    
    ggPS2
    
  })
  
  selID <- reactive(input$selID)
  
  pcpDat <- reactive(data[which(data$ID %in% selID()), ])
  
  # Print the selected gene IDs
  output$selectedValues = DT::renderDataTable(pcpDat(), rownames= FALSE)
  
  colNms <- colnames(data[, -1])
  nVar <- ncol(data)
  
  boxDat <- data %>% gather(Sample, Count, -c(ID))
  BP <- ggplot(boxDat, aes(x = Sample, y = Count)) + geom_boxplot()
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
      
      }", data = list(pcpDat = pcpDat(), nVar = nVar, colNms = colNms))}) # should be p$nrow not 7

  output$downloadData <- downloadHandler(
    filename = function() {
      paste("SM.csv")
    },
    content = function(file) {
      write.table(pcpDat()[["ID"]], file, row.names = FALSE, col.names = FALSE, quote = FALSE, sep=',')
    }
  )
  
  }