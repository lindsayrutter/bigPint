#' @title Plot static parallel coordinate clusters
#' 
#' @description Perform hierarchical clustering analysis and visualize results with parallel coordinate plots. Optionally, save gene IDs within each cluster to .rds files for later use.
#' 
#' @param data DATA FRAME | Read counts
#' @param dataMetrics LIST | Differential expression metrics; default NULL
#' @param saveFile BOOLEAN [TRUE | FALSE] | Save file to outDir; default TRUE
#' @param nC INTEGER | Number of clusters; default 4
#'@param threshVar CHARACTER STRING | Name of column in dataMetrics object that is used to threshold significance; default "FDR"
#' @param threshVal INTEGER | Maximum value to threshold significance from threshVar object; default 0.05
#' @param lineSize INTEGER | Size of plotted parallel coordinate lines; default 0.1
#' @param lineAlpha INTEGER | Alpha value of plotted parallel coordinate lines, default 0.5
#' @param colList CHARACTER ARRAY | List of colors for each cluster; default is rainbow(nC)
#' @param clusterAllData BOOLEAN [TRUE | FALSE] | Create clusters based on the whole dataset and then assign significant genes to those clusters; default is TRUE. If FALSE, create clusters based on just the significant genes. With either option, the side-by-side boxplot will represent the whole dataset (from data input) and the parallel coordinate lines will represent only the significant genes (those that pass threshVal for threshVar)
#' @param aggMethod CHARACTER STRING ["ward.D" | "ward.D2" | "single" | "complete" | "average" | "mcquitty" | "median" | "centroid"] | The agglomeration method to be used in the hierarchical clustering; default "ward.D"
#' @param xAxisLabel CHARACTER STRING | Horizontal axis label; default "Sample"
#' @param yAxisLabel CHARACTER STRING | Vertical axis label; default "Count"
#' @param vxAxis BOOLEAN [TRUE | FALSE] | Flip x-axis text labels to vertical orientation; default FALSE
#' @param verbose BOOLEAN [TRUE | FALSE] | Print each cluster from each cluster size into separate files and print the associated IDs of each cluster from each cluster size into separate .rds files; default is FALSE
#' @param geneList CHARACTER ARRAY | List of ID values of genes to be drawn from data as parallel coordinate lines. Use this parameter if you have predetermined genes to be drawn. Otherwise, use dataMetrics, threshVar, and threshVal to create genes to be overlaid as parallel coordinate lines; default NULL
#' @param outDir CHARACTER STRING | Output directory to save all images; default current directory
#' @importFrom dplyr %>%
#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot element_text
#' @importFrom gridExtra arrangeGrob
#' @importFrom reshape melt
#' @importFrom grDevices rainbow
#' @importFrom graphics plot
#' @importFrom utils write.table
#' @importFrom grid grid.draw
#' @seealso
#' \code{\link[stats]{hclust}}
#' @export
#' @examples
#' # Example 1: Perform hierarchical clustering of size four using the default agglomeration 
#' # method "ward.D". Cluster only on the genes that have FDR < 1e-7 (n = 113) and overlay 
#' # these genes.
#' 
#' library(grid)
#' library(matrixStats)
#' library(ggplot2)
#' data(soybean_ir_sub)
#' soybean_ir_sub[,-1] <- log(soybean_ir_sub[-1]+1)
#' data(soybean_ir_sub_metrics)
#' colList = c("#00A600FF", rainbow(5)[c(1,4,5)])
#' ret <- plotClusters(data=soybean_ir_sub, dataMetrics = soybean_ir_sub_metrics, nC=4, 
#'   colList = colList, clusterAllData = FALSE, threshVal = 1e-7, saveFile = FALSE)
#' plot(ret[["N_P_4"]])
#' 
#' # Example 2: Perform the same analysis, only now create the four groups by clustering on all 
#' # genes in the data (n = 5,604). Then, overlay the genes that have FDR < 1e-7 (n = 113) into 
#' # their corresponding clusters.
#' 
#' ret <- plotClusters(data=soybean_ir_sub, dataMetrics = soybean_ir_sub_metrics, nC=4, 
#'   colList = colList, clusterAllData = TRUE, threshVal = 1e-7, saveFile = FALSE)
#' plot(ret[["N_P_4"]])
#' 
#' # Example 3: Perform the same analysis, only now overlay all genes in the data by keeping the
#' # dataMetrics object as its default value of NULL.
#' 
#' ret <- plotClusters(data=soybean_ir_sub, nC=4, colList = colList, clusterAllData = TRUE,
#'   saveFile = FALSE)
#' plot(ret[["N_P_4"]])
#' 
#' # Example 4: Visualization of gene clusters is usually performed on standardized data. Here, 
#' # hierarchical clustering of size four is performed using the agglomeration method "average"
#' # on standardized data. Only genes with FDR < 0.05 are used for the clustering. Only two of
#' # the three pairwise combinations of treatment groups (S1 and S2; S1 and S3) have any genes
#' # with FDR < 0.05. The output plots for these two pairs are examined. 
#' 
#' data(soybean_cn_sub)
#' data(soybean_cn_sub_metrics)
#' soybean_cn_sub_st <- as.data.frame(t(apply(as.matrix(soybean_cn_sub[,-1]), 1, scale)))
#' soybean_cn_sub_st$ID <- as.character(soybean_cn_sub$ID)
#' soybean_cn_sub_st <- soybean_cn_sub_st[,c(length(soybean_cn_sub_st),
#'   1:length(soybean_cn_sub_st)-1)]
#' colnames(soybean_cn_sub_st) <- colnames(soybean_cn_sub)
#' nID <- which(is.nan(soybean_cn_sub_st[,2]))
#' soybean_cn_sub_st[nID,2:length(soybean_cn_sub_st)] <- 0
#' ret <- plotClusters(data=soybean_cn_sub_st, dataMetrics = soybean_cn_sub_metrics, nC=4, 
#'   colList = c("#00A600FF", "#CC00FFFF", "red", "darkorange"), lineSize = 0.5, lineAlpha = 1, 
#'   clusterAllData = FALSE, aggMethod = "average", yAxisLabel = "Standardized read count", 
#'   saveFile = FALSE)
#' names(ret)
#' plot(ret[["S1_S2_4"]])
#' plot(ret[["S1_S3_4"]])
#' 
#' # Example 5: Run the same analysis, only now set the verbose parameter to value TRUE. This 
#' # will save images of each individual cluster, .rds files that contain the IDs within each 
#' # cluster, and images of the conglomerate clusters to outDir (default current working
#' # directory).
#' 
#' \dontrun{
#' plotClusters(data=soybean_cn_sub_st, dataMetrics = soybean_cn_sub_metrics, nC=4,
#'   colList = c("#00A600FF", "#CC00FFFF", "red", "darkorange"), lineSize = 0.5, lineAlpha = 1, 
#'   clusterAllData = FALSE, aggMethod = "average", yAxisLabel = "Standardized read count", 
#'   verbose = TRUE)
#' }
#' 

plotClusters <- function(data, dataMetrics = NULL, nC = 4, threshVar="FDR", threshVal=0.05, outDir=getwd(), colList = rainbow(nC), aggMethod = "ward.D", yAxisLabel = "Count", xAxisLabel = "Sample", lineSize = 0.1, lineAlpha = 0.5, clusterAllData = TRUE, verbose=FALSE, saveFile = TRUE, vxAxis = FALSE, geneList = NULL){

  key <- val <- ID <- rainbow <- NULL
  
  colNames <- colnames(data)
  myPairs <- unique(sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1]))
  myPairs <- myPairs[-which(myPairs=="ID")]
  colGroups <- sapply(colNames, function(x) unlist(strsplit(x,"[.]"))[1])
  
  ret = list()
  for (i in 1:(length(myPairs)-1)){
    for (j in (i+1):length(myPairs)){
      group1 = myPairs[i]
      group2 = myPairs[j]
      fData <- cbind(ID=data$ID, data[,which(colGroups %in% c(group1, group2))])

      boxDat <- fData %>% gather(key, val, c(-ID))
      colnames(boxDat) <- c("ID", "Sample", "Count")
      
      userOrder <- unique(boxDat$Sample)
      boxDat$Sample <- as.factor(boxDat$Sample)
      levels(boxDat$Sample) <- userOrder
      
      if (!is.null(dataMetrics) && is.null(geneList)){
        metricPair = dataMetrics[[paste0(group1,"_",group2)]]
        metricPair = metricPair[order(metricPair[threshVar]),]
        threshID <- metricPair[which(metricPair[threshVar] <= threshVal),]$ID
        cData <- fData[which(fData$ID %in% threshID),]  
      } 
      else if (!is.null(geneList)){
        cData <- fData[which(fData$ID %in% geneList),]  
      }
      else{
        cData <- fData
      }
      
      plotName <- paste0(group1,"_",group2)
      
      # geneList variable takes precedence. Create metricPair to only select geneList IDs
      if(!is.null(geneList)){
        metricPair = data.frame(ID = data$ID, FDR = 1)
        metricPair$ID = as.character(metricPair$ID)
        metricPair[which(metricPair$ID %in% geneList), ]$FDR = 0.1
        threshVar = "FDR"
        threshVal = 0.5
        dataMetrics = 1
      }
      
      # Check if there are even any genes that pass the threshold. Then, perform clustering on whole dataset and assigned significant genes to those clusters from the full dataset.
      if (nrow(data)>=nC && clusterAllData == TRUE){
        dendo = data
        rownames(dendo) = NULL
        d = suppressWarnings(dist(as.matrix(dendo)))
        hC = hclust(d, method=aggMethod)
        
        #colList = rainbow(nC)
        k = cutree(hC, k=nC)
        ###########################
        plot_clusters = lapply(1:nC, function(j){
          i = rev(order(table(k)))[j]
          x = as.data.frame(data[which(k==i),])
          x$cluster = "color"
          x$cluster2 = factor(x$cluster)
          xNames = x$ID
          if (!is.null(dataMetrics)){
            metricFDR = metricPair[which(as.character(metricPair$ID) %in% xNames),]
            sigID = metricFDR[which(metricFDR[[threshVar]]<=threshVal),]$ID
            xSig = x[which(xNames %in% sigID),]
            xSigNames = rownames(xSig)
            nGenes = nrow(xSig)            
          }
          else{
            xSig = x
            xSigNames = rownames(x)
            nGenes = nrow(x)
          }

          if (verbose==TRUE){
            IDs = as.character(sigID)
            saveRDS(IDs, file = paste(outDir, "/", plotName, "_", nC, "_", j, ".rds", sep=""))
          }
          
          xSig$ID = xSigNames
          if (nrow(xSig)>0){
            keepCol = c(1, which(sapply(colnames(xSig), function(x) unlist(strsplit(x,"[.]"))[1]) %in% c(group1, group2)))
            xSig = xSig[, keepCol]
            pcpDat <- melt(xSig, id.vars="ID")
            colnames(pcpDat) <- c("ID", "Sample", "Count")
            pcpDat$Sample <- as.character(pcpDat$Sample)
            pcpDat$ID <- as.factor(pcpDat$ID)
            
            p <- ggplot(boxDat, aes_string(x = 'Sample', y = 'Count')) + geom_boxplot() + geom_line(data=pcpDat, aes_string(x = 'Sample', y = 'Count', group = 'ID'), colour = colList[j], alpha=lineAlpha, size = lineSize) + ylab(yAxisLabel) + xlab(xAxisLabel) + ggtitle(paste("Cluster ", j, " Genes (n=", format(nGenes, big.mark=",", scientific=FALSE), ")",sep="")) + theme(plot.title = element_text(hjust = 0.5, size=14, face="plain"), axis.text=element_text(size=11), axis.title=element_text(size=14))
          }
          else{
            p <- ggplot(boxDat, aes_string(x = 'Sample', y = 'Count')) + geom_boxplot() + ylab(yAxisLabel) + xlab(xAxisLabel) + ggtitle(paste("Cluster ", j, " Genes (n=", format(nGenes, big.mark=",", scientific=FALSE), ")",sep="")) + theme(plot.title = element_text(hjust = 0.5, size=14, face="plain"), axis.text=element_text(size=11), axis.title=element_text(size=14))            
          }
          if (vxAxis == TRUE){
            p <- p + theme(axis.text.x = element_text(angle=90, hjust=1))
          }
          if (verbose==TRUE){
            fileName = paste(outDir, "/", plotName, "_", nC, "_", j, ".jpg", sep="")
            jpeg(fileName)
            plot(p)
            invisible(dev.off()) 
          }
          p
        })
        p = arrangeGrob(grobs=plot_clusters, ncol=2)
        
        if (saveFile == TRUE || verbose == TRUE){
          fileName = paste(outDir, "/", plotName, "_", nC, ".jpg", sep="")
          jpeg(fileName)
          grid.draw(p)
          invisible(dev.off())
        }
        ret[[paste0(plotName, "_", nC)]] = p
      }
      # Check if there are even any genes that pass the threshold. Then, perform clustering on just the significant genes.
      if (nrow(cData)>=nC && clusterAllData == FALSE){
        dendo = cData
        rownames(dendo) = NULL
        d = suppressWarnings(dist(as.matrix(dendo)))
        hC = hclust(d, method=aggMethod)
        
        #colList = rainbow(nC)
        k = cutree(hC, k=nC)
        ###########################
        plot_clusters = lapply(1:nC, function(j){
          i = rev(order(table(k)))[j]
          x = as.data.frame(cData[which(k==i),])
          x$cluster = "color"
          x$cluster2 = factor(x$cluster)
          x$ID = factor(x$ID)
          xNames = x$ID

          if (!is.null(dataMetrics)){
            metricFDR = metricPair[which(as.character(metricPair$ID) %in% xNames),]
            sigID = metricFDR[which(metricFDR[[threshVar]]<=threshVal),]$ID
            xSig = x[which(xNames %in% sigID),]
            xSigNames = rownames(xSig)
            nGenes = nrow(xSig)            
          }
          else{
            xSig = x
            xSigNames = rownames(x)
            nGenes = nrow(x)
          }
          
          if (verbose==TRUE){
            IDs = as.character(xNames)
            saveRDS(IDs, file = paste(outDir, "/", plotName, "_", nC, "_", j, ".rds", sep=""))
          }
          
          xSig$ID = xSigNames
          
          keepCol = c(1, which(sapply(colnames(xSig), function(x) unlist(strsplit(x,"[.]"))[1]) %in% c(group1, group2)))
          xSig = xSig[, keepCol]
          
          pcpDat <- melt(xSig, id.vars="ID")
          colnames(pcpDat) <- c("ID", "Sample", "Count")
          pcpDat$Sample <- as.character(pcpDat$Sample)
          pcpDat$ID <- as.factor(pcpDat$ID)
          
          p <- ggplot(boxDat, aes_string(x = 'Sample', y = 'Count')) + geom_boxplot() + geom_line(data=pcpDat, aes_string(x = 'Sample', y = 'Count', group = 'ID'), colour = colList[j], alpha=lineAlpha, size = lineSize) + ylab(yAxisLabel) + xlab(xAxisLabel) + ggtitle(paste("Cluster ", j, " Genes (n=", format(nGenes, big.mark=",", scientific=FALSE), ")",sep="")) + theme(plot.title = element_text(hjust = 0.5, size=14, face="plain"), axis.text=element_text(size=11), axis.title=element_text(size=14))
          
          if (vxAxis == TRUE){
            p <- p + theme(axis.text.x = element_text(angle=90, hjust=1))
          }
          
          if (verbose==TRUE){
            fileName = paste(outDir, "/", plotName, "_", nC, "_", j, ".jpg", sep="")
            jpeg(fileName)
            plot(p)
            invisible(dev.off()) 
          }
          p
        })
        p = arrangeGrob(grobs=plot_clusters, ncol=2)
        
        if (saveFile == TRUE || verbose == TRUE){
          fileName = paste(outDir, "/", plotName, "_", nC, ".jpg", sep="")
          jpeg(fileName)
          grid.draw(p)
          invisible(dev.off())
        }
        ret[[paste0(plotName, "_", nC)]] = p
      }
      
      # Indicate if no significant genes existed
      if (nrow(data)<nC && clusterAllData == TRUE){
        print("Not enough data to cluster by that many groups")
      }
      
      # Indicate if not enough significant genes existed
      if (nrow(cData)==0){
        print(paste0(group1, "_", group2, ": There were no significant genes"))
      }
      else if (nrow(cData)<nC && clusterAllData == FALSE){
        print(paste0(group1, "_", group2, ": Not enough significant genes (", nrow(cData), ",) to cluster by that many groups (", nC, ")"))
      }
      
    }
  }
  invisible(ret)
}
