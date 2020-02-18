#' @title Plot static parallel coordinate clusters
#' 
#' @description Perform hierarchical clustering analysis and visualize 
#' results with parallel coordinate plots. Optionally, save gene IDs within 
#' each cluster to .rds files for later use.
#' 
#' @param data DATA FRAME | Read counts
#' @param dataMetrics LIST | Differential expression metrics; default NULL
#' @param dataSE SUMMARIZEDEXPERIMENT | Summarized experiment format that
#' can be used in lieu of data and dataMetrics; default NULL
#' @param geneList CHARACTER ARRAY | Array of ID values of genes to be drawn 
#' from data as parallel coordinate lines. Use this parameter if you have 
#' predetermined genes to be drawn. These genes will be clustered. Otherwise, 
#' use dataMetrics, threshVar, and threshVal to create clusters to be
#' overlaid as parallel coordinate lines; default NULL. See package website
#' for examples
#' @param geneLists LIST | List of ID values of genes already clustered to be #' drawn from data as parallel coordinate lines. Each list item is an array
#' of genes ID values that are already grouped as a cluster. Unlike the 
#' singular geneList object, the plural geneLists object is not be clustered.
#' If you instead wish to cluster genes, use dataMetrics, threshVar, and 
#' threshVal or geneList to create clusters to be overlaid as parallel
#' coordinate lines; default NULL. See package website for examples
#' @param threshVar CHARACTER STRING | Name of column in dataMetrics object 
#' that is used to threshold significance; default "FDR"
#' @param threshVal INTEGER | Maximum value to threshold significance from 
#' threshVar object; default 0.05
#' @param clusterAllData BOOLEAN [TRUE | FALSE] | Create clusters based on 
#' the whole dataset and then assign significant genes to those clusters; 
#' default is TRUE. If FALSE, create clusters based on just the significant 
#' genes. With either option, the side-by-side boxplot will represent the 
#' whole dataset (from data input) and the parallel coordinate lines will 
#' represent only the significant genes (those that pass threshVal for 
#' threshVar)
#' @param nC INTEGER | Number of clusters; default 4
#' @param colList CHARACTER ARRAY | List of colors for each cluster; default 
#' is rainbow(nC)
#' @param aggMethod CHARACTER STRING ["ward.D" | "ward.D2" | "single" | 
#' "complete" | "average" | "mcquitty" | "median" | "centroid"] | The 
#' agglomeration method to be used in the hierarchical clustering; default 
#' "ward.D"
#' @param yAxisLabel CHARACTER STRING | Vertical axis label; default "Count"
#' @param xAxisLabel CHARACTER STRING | Horizontal axis label; default 
#' "Sample"
#' @param lineSize INTEGER | Size of plotted parallel coordinate lines; 
#' default 0.1
#' @param lineAlpha INTEGER | Alpha value of plotted parallel coordinate 
#' lines, default 0.5
#' @param vxAxis BOOLEAN [TRUE | FALSE] | Flip x-axis text labels to vertical
#' orientation; default FALSE
#' @param outDir CHARACTER STRING | Output directory to save all images; 
#' default tempdir()
#' @param saveFile BOOLEAN [TRUE | FALSE] | Save file to outDir; default TRUE
#' @param verbose BOOLEAN [TRUE | FALSE] | Print each cluster from each 
#' cluster size into separate files and print the associated IDs of each 
#' cluster from each cluster size into separate .rds files; default is FALSE
#' @importFrom dplyr %>%
#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot element_text
#' @importFrom gridExtra arrangeGrob
#' @importFrom reshape melt
#' @importFrom grDevices rainbow
#' @importFrom graphics plot
#' @importFrom utils write.table
#' @importFrom grid grid.draw
#' @importFrom utils combn
#' @seealso
#' \code{\link[stats]{hclust}}
#' \url{https://lindsayrutter.github.io/bigPint/articles/clusters.html}
#' @return List of n elements each containing a grid of parallel coordinate 
#' plots, where n is the number of treatment pair combinations in the data 
#' object. If the saveFile parameter has a value of TRUE, then each grid of 
#' parallel coordinate plots is saved to the location specified in the outDir
#' parameter as a JPG file. If the verbose parameter has a value of TRUE, 
#' then a JPG file for each parallel coordinate plot in each grid, RDS file 
#' containing the superimposed IDs for each parallel coordinate plot in each 
#' grid, and the JPG file of each grid of parallel coordinate plots is saved 
#' to the location specified in the outDir parameter.
#' @export
#' @examples
#' # Example 1: Perform hierarchical clustering of size four using the 
#' # default agglomeration method "ward.D". Cluster only on the genes that have
#' # FDR < 1e-7 (n = 113) and overlay these genes.
#' 
#' library(grid)
#' library(matrixStats)
#' library(ggplot2)
#' data(soybean_ir_sub)
#' soybean_ir_sub[,-1] <- log(soybean_ir_sub[-1]+1)
#' data(soybean_ir_sub_metrics)
#' colList = c("#00A600FF", rainbow(5)[c(1,4,5)])
#' ret <- plotClusters(data=soybean_ir_sub,
#'     dataMetrics = soybean_ir_sub_metrics, nC=4, colList = colList,
#'     clusterAllData = FALSE, threshVal = 1e-7, saveFile = FALSE)
#' grid.draw(ret[["N_P_4"]])
#' 
#' # Example 2: Perform the same analysis, only now create the four groups by 
#' # clustering on all genes in the data (n = 5,604). Then, overlay the genes 
#' # that have FDR < 1e-7 (n = 113) into their corresponding clusters.
#' 
#' ret <- plotClusters(data=soybean_ir_sub,
#'     dataMetrics = soybean_ir_sub_metrics, nC=4, colList = colList,
#'     clusterAllData = TRUE, threshVal = 1e-7, saveFile = FALSE)
#' grid.draw(ret[["N_P_4"]])
#' 
#' # Example 3: Perform the same analysis, only now overlay all genes in the 
#' # data by keeping the dataMetrics object as its default value of NULL.
#' 
#' ret <- plotClusters(data=soybean_ir_sub, nC=4, colList = colList,
#'     clusterAllData = TRUE, saveFile = FALSE)
#' grid.draw(ret[["N_P_4"]])
#' 
#' # Example 4: Visualization of gene clusters is usually performed on
#' # standardized data. Here, hierarchical clustering of size four is performed
#' # using the agglomeration method "average" on standardized data. Only genes 
#' # with FDR < 0.05 are used for the clustering. Only two of the three 
#' # pairwise combinations of treatment groups (S1 and S2; S1 and S3) have any 
#' # genes with FDR < 0.05. The output plots for these two pairs are examined. 
#' 
#' data(soybean_cn_sub)
#' data(soybean_cn_sub_metrics)
#' soybean_cn_sub_st <- as.data.frame(t(apply(as.matrix(soybean_cn_sub[,-1]),
#'     1, scale)))
#' soybean_cn_sub_st$ID <- as.character(soybean_cn_sub$ID)
#' soybean_cn_sub_st <- soybean_cn_sub_st[,c(length(soybean_cn_sub_st),
#'     1:length(soybean_cn_sub_st)-1)]
#' colnames(soybean_cn_sub_st) <- colnames(soybean_cn_sub)
#' nID <- which(is.nan(soybean_cn_sub_st[,2]))
#' soybean_cn_sub_st[nID,2:length(soybean_cn_sub_st)] <- 0
#' ret <- plotClusters(data=soybean_cn_sub_st,
#'     dataMetrics = soybean_cn_sub_metrics, nC=4,
#'     colList = c("#00A600FF", "#CC00FFFF", "red", "darkorange"),
#'     lineSize = 0.5, lineAlpha = 1, clusterAllData = FALSE,
#'     aggMethod = "average", yAxisLabel = "Standardized read count",
#'     saveFile = FALSE)
#' names(ret)
#' grid.draw(ret[["S1_S2_4"]])
#' grid.draw(ret[["S1_S3_4"]])
#' 
#' # Example 5: Run the same analysis, only now set the verbose parameter to 
#' # value TRUE. This will save images of each individual cluster, .rds files 
#' # that contain the IDs within each cluster, and images of the conglomerate 
#' # clusters to outDir (default tempdir()).
#' 
#' \dontrun{
#' plotClusters(data=soybean_cn_sub_st, dataMetrics = soybean_cn_sub_metrics,
#'   nC=4, colList = c("#00A600FF", "#CC00FFFF", "red", "darkorange"),
#'   lineSize = 0.5, lineAlpha = 1, clusterAllData = FALSE,
#'   aggMethod = "average", yAxisLabel = "Standardized read count",
#'   verbose = TRUE)
#' }
#' 
plotClusters <- function(data, dataMetrics = NULL, dataSE=NULL, geneList = NULL,
    geneLists = NULL, threshVar="FDR", threshVal=0.05, clusterAllData = TRUE, 
    nC = 4, colList = rainbow(nC), aggMethod = c("ward.D", "ward.D2",
    "single", "complete", "average", "mcquitty", "median", "centroid"),
    yAxisLabel = "Count", xAxisLabel = "Sample", lineSize = 0.1,
    lineAlpha = 0.5, vxAxis = FALSE, outDir=tempdir(), saveFile = TRUE,
    verbose=FALSE){

aggMethod <- match.arg(aggMethod)

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
    }
}

# Check that input parameters fit required formats
helperTestData(data)
if (is.null(geneList) && !is.null(dataMetrics)){
    helperTestDataMetrics(data, dataMetrics, threshVar)
}

key <- val <- ID <- rainbow <- NULL
myPairs <- helperMakePairs(data)[["myPairs"]]
colGroups <- helperMakePairs(data)[["colGroups"]]
cols.combn <- combn(myPairs, 2, simplify = FALSE)

if (!is.null(geneLists)){
    nC = length(geneLists)
    colList = rainbow(nC)
    seqVec = seq(nC)
}

ret <- lapply(cols.combn, function(x){
    group1 = x[1]
    group2 = x[2]
    fData <- cbind(ID=data$ID, data[,which(colGroups %in% 
    c(group1, group2))])
    
    boxDat <- fData %>% gather(key, val, c(-ID))
    colnames(boxDat) <- c("ID", "Sample", "Count")
    
    userOrder <- unique(boxDat$Sample)
    boxDat$Sample <- as.factor(boxDat$Sample)
    levels(boxDat$Sample) <- userOrder
    
    plotName <- paste0(group1,"_",group2)
    
    # If geneLists is not NULL
    # if (!is.null(geneLists)){
    # nC = length(geneLists)
    # colList = rainbow(nC)
    # seqVec = seq(nC)
    # 
    # plot_clusters = lapply(seq_along(seqVec), function(j){
    #     xSigNames = geneLists[[j]]
    #     xSig = fData %>% filter(ID %in% xSigNames)
    #     nGenes = nrow(xSig)
    #     
    #     pcpDat <- melt(xSig, id.vars="ID")
    #     colnames(pcpDat) <- c("ID", "Sample", "Count")
    #     pcpDat$Sample <- as.character(pcpDat$Sample)
    #     pcpDat$ID <- as.factor(pcpDat$ID)
    #     
    #     p <- ggplot(boxDat, aes_string(x = 'Sample', y = 'Count')) + geom_boxplot() + geom_line(data=pcpDat, aes_string(x = 'Sample', y = 'Count', group = 'ID'), colour = colList[j], alpha=lineAlpha, size = lineSize) + ylab(yAxisLabel) + xlab(xAxisLabel) + ggtitle(paste("Cluster ", j, " Genes (n=", format(nGenes, big.mark=",", scientific=FALSE), ")", sep="")) + theme(plot.title = element_text(hjust = 0.5, size=14, face="plain"), axis.text=element_text(size=11), axis.title=element_text(size=14))
    #     p
    # })
    # p = arrangeGrob(grobs=plot_clusters, ncol=2)
    # if (saveFile == TRUE || verbose == TRUE){
    #     fileName = paste(outDir, "/", plotName, "_", nC, ".jpg", sep="")
    #     jpeg(fileName)
    #     grid.draw(p)
    #     invisible(dev.off())}
    # #return(p)
    # }
    
    if (!is.null(dataMetrics) && is.null(geneList)){
        metricPair = dataMetrics[[paste0(group1,"_",group2)]]
        metricPair = metricPair[order(metricPair[threshVar]),]
        threshID <- metricPair[which(metricPair[threshVar] <=
        threshVal),]$ID
        cData <- fData[which(fData$ID %in% threshID),]  
    } 
    else if (!is.null(geneList)){
        cData <- fData[which(fData$ID %in% geneList),]  
    }
    else{
        cData <- fData
    }
    
    if (!is.null(geneLists)){
    nC = length(geneLists)
    colList = rainbow(nC)
    seqVec = seq(nC)

    plot_clusters = lapply(seq_along(seqVec), function(j){
        xSigNames = geneLists[[j]]
        xSig = fData %>% filter(ID %in% xSigNames)
        nGenes = nrow(xSig)

        pcpDat <- melt(xSig, id.vars="ID")
        colnames(pcpDat) <- c("ID", "Sample", "Count")
        pcpDat$Sample <- as.character(pcpDat$Sample)
        pcpDat$ID <- as.factor(pcpDat$ID)

        p <- ggplot(boxDat, aes_string(x = 'Sample', y = 'Count')) + geom_boxplot() + geom_line(data=pcpDat, aes_string(x = 'Sample', y = 'Count', group = 'ID'), colour = colList[j], alpha=lineAlpha, size = lineSize) + ylab(yAxisLabel) + xlab(xAxisLabel) + ggtitle(paste("Cluster ", j, " Genes (n=", format(nGenes, big.mark=",", scientific=FALSE), ")", sep="")) + theme(plot.title = element_text(hjust = 0.5, size=14, face="plain"), axis.text=element_text(size=11), axis.title=element_text(size=14))
        p
    })
    p = arrangeGrob(grobs=plot_clusters, ncol=2)
    if (saveFile == TRUE || verbose == TRUE){
        fileName = paste(outDir, "/", plotName, "_", nC, ".jpg", sep="")
        jpeg(fileName)
        grid.draw(p)
        invisible(dev.off())}
    return(p)
    }
    
    
    # geneList variable takes precedence. Create metricPair to only
    # select geneList IDs
    if(!is.null(geneList)){
        metricPair = data.frame(ID = data$ID, FDR = 1)
        metricPair$ID = as.character(metricPair$ID)
        metricPair[which(metricPair$ID %in% geneList), ]$FDR = 0.1
        threshVar = "FDR"
        threshVal = 0.5
        dataMetrics = 1
    }
    
    # Check if there are even any genes that pass the threshold. 
    # Then, perform clustering on whole dataset and assigned 
    # significant genes to those clusters from the full dataset.
    if (nrow(data)>=nC && clusterAllData == TRUE){
        p <- helperCadTRUE(data, dataMetrics, metricPair, aggMethod,
        nC, threshVar, threshVal, verbose, vxAxis, saveFile, boxDat,
        xAxisLabel, yAxisLabel, lineAlpha, lineSize, plotName,
        outDir, colList)
        return(p)
    }
    
    # Check if there are even any genes that pass the threshold. Then, 
    # perform clustering on just the significant genes.
    if (nrow(cData)>=nC && clusterAllData == FALSE){
        p <- helperCadFALSE(cData, dataMetrics, metricPair,
        aggMethod, nC, threshVar, threshVal, verbose, vxAxis,
        saveFile, boxDat, xAxisLabel, yAxisLabel, lineAlpha, lineSize,
        plotName, outDir, colList)
        return(p)
    }
    
    # Indicate if no significant genes existed
    if (nrow(data)<nC && clusterAllData == TRUE){
        print("Not enough data to cluster by that many groups")
    }
    
    # Indicate if not enough significant genes existed
    if (nrow(cData)==0){
        print(paste0(group1, "_", group2,
        ": There were no significant genes"))
    }
    else if (nrow(cData)<nC && clusterAllData == FALSE){
        print(paste0(group1, "_", group2,
        ": Not enough significant genes (" ,nrow(cData),
        ",) to cluster by that many groups (", nC, ")"))
    }
})

plotNames <- lapply(cols.combn, function(x){
    group1 = x[1]
    group2 = x[2]
    paste0(group1,"_",group2,"_",nC)
})

names(ret) <- plotNames
invisible(ret)
}
