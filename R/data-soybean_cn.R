#' Normalized soybean cotyledon data
#' 
#' @name soybean_cn
#' @title Normalized soybean cotyledon data
#' @description This dataset contains normalized RNA-sequencing read counts 
#' from soybean cotyledon across three time stages of development. Early 
#' stage cotyledons were collected four days after planting and were green 
#' but closed. Middle stage cotyledons were collected while green and open, 
#' soon after the plant generated its first set of unifoliate leaves. Late 
#' stage cotyledons were collected immediately after the initiation of 
#' yellowing and shrinking.
#' @docType data
#' @format a \code{RData} instance, 1 row per gene
#' @details \itemize{
    #' \item ID gene name
    #' \item S1.1 early stage replicate 1 normalized read counts
    #' \item S1.2 early stage replicate 2 normalized read counts
    #' \item S1.3 early stage replicate 3 normalized read counts
    #' \item S2.1 middle stage replicate 1 normalized read counts
    #' \item S2.2 middle stage replicate 2 normalized read counts
    #' \item S2.3 middle stage replicate 3 normalized read counts
    #' \item S3.1 late stage replicate 1 normalized read counts
    #' \item S3.2 late stage replicate 2 normalized read counts
    #' \item S3.3 late stage replicate 3 normalized read counts
    #' }
#' @docType data
#' @keywords datasets
#' @format A data frame with 73,320 rows and 10 variables
#' @references
#' Brown AV, Hudson KA (2015) Developmental profiling of gene expression in 
#' soybean trifoliate leaves and cotyledons. BMC Plant Biol 15:169
NULL
