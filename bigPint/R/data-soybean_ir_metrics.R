#' Raw soybean leaves metrics
#'
#' @name soybean_ir_metrics
#' @title Raw soybean leaves iron-metabolism metrics
#' @description This data contains metrics for raw RNA-sequencing read counts
#' from a soybean dataset that compared leaves that were exposed to iron-rich
#' (iron-postive) soil conditions versus leaves that were exposed to iron
#' -poor (iron-negative) soil conditions. The data was collected 120 minutes 
#' after iron conditions were initiated. The metrics include the log fold 
#' change and the p-values for all genes and all pairwise combinations of 
#' treatment groups.
#' @docType data
#' @format a \code{RData} instance, 1 list per treatment group combination 
#' and 1 row per gene 
#' @details \itemize{
    #' \item ID gene name
    #' \item logFC log fold change
    #' \item PValue p-value
    #' }
#' @docType data
#' @keywords datasets
#' @format A nested list of length 1. The list contains the metrics for the 
#' 56,044 genes for the one treatment group combination.
#' @seealso \code{\link{soybean_ir}} for information about the treatment 
#' groups
NULL
