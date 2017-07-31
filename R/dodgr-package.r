#' dodgr.
#'
#' @name dodgr
#' @docType package
#' @importFrom igraph distances E make_directed_graph
#' @importFrom magrittr %>%
#' @importFrom osmdata add_feature getbb opq osmdata_sf
#' @importFrom rbenchmark benchmark
#' @importFrom sp bbox point.in.polygon
#' @importFrom Rcpp evalCpp
#' @useDynLib dodgr, .registration = TRUE
NULL

#' weighting_profiles
#'
#' Collection of weighting profiles used to adjust the routing process to
#' different means of transport. Original data taken from the Routino project.
#'
#' @name weighting_profiles
#' @docType data
#' @keywords datasets
#' @format \code{data.frame} with profile names, means of transport and
#' weights.
#' @references \url{https://www.routino.org/xml/routino-profiles.xml}
NULL
