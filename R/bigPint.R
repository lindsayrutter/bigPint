#' \code{bigPint} package
#'
#' bigPint R API
#'
#' See the README on
#' \href{https://github.com/lindsayrutter/bigPint#readme}{GitHub}
#'
#' @docType package
#' @name bigPint
NULL

## quiets concerns of R CMD check regarding variables that appear in pipelines
if(getRversion() >= "2.15.1") utils::globalVariables(c("mtcars", "wt", "mpg"))