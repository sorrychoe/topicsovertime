#' @useDynLib topicsovertime, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats predict
#' @importFrom utils globalVariables
#' @importMethodsFrom Matrix coerce
"_PACKAGE"

# The .data pronoun used in autoplot.tot() is supplied by ggplot2's data mask.
globalVariables(".data")
