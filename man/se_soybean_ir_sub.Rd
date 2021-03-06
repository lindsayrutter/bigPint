% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data-se_soybean_ir_sub.R
\docType{data}
\name{se_soybean_ir_sub}
\alias{se_soybean_ir_sub}
\title{Raw and subsetted soybean leaves iron-metabolism SummarizedExperiment}
\format{
(1) A DelayedMatrix with 5,604 rows and 7 variables and 1 row per
gene. Accessible by \code{assay(se_soybean_ir_sub)}. (2) A List with 1
list element per treatment group combination and 1 row per gene.
Accessible by \code{rowData(se_soybean_ir_sub)}.

A data frame with 5,604 rows and 7 variables
}
\usage{
data(se_soybean_ir_sub)
}
\description{
This dataset contains raw RNA-sequencing read counts from a 
soybean dataset that compared leaves that were exposed to iron-rich (iron
-postive) soil conditions versus leaves that were exposed to iron-poor 
(iron-negative) soil conditions. The data was collected 120 minutes after 
iron conditions were initiated. To save on size, this example dataset was 
generated by obtaining a random subset of 1 out of 10 genes from the 
original resource.
}
\details{
Raw and subsetted soybean leaves data SummarizedExperiment object


\itemize{
\item \code{assay(se_soybean_ir_sub)} Structure:
\item N.1 iron-negative condition replicate 1 raw read counts
\item N.2 iron-negative condition replicate 2 raw read counts
\item N.3 iron-negative condition replicate 3 raw read counts
\item P.1 iron-positive condition replicate 1 raw read counts
\item P.2 iron-positive condition replicate 2 raw read counts
\item P.3 iron-positive condition replicate 3 raw read counts
}

\itemize{
\item \code{rowData(se_soybean_ir_sub)} list element Structure:
\item ID gene name
\item logFC log fold change
\item logCPM log counts per million
\item LR likelihood ratio
\item PValue p-value
\item FDR FDR value
}
}
\references{
Moran Lauter AN, Graham MA. NCBI SRA bioproject accession: PRJNA318409.
}
\seealso{
\code{\link{soybean_ir_sub}} and
\code{\link{soybean_ir_sub_metrics}}
}
\keyword{datasets}
