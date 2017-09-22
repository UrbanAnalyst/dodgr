#' dodgr.
#'
#' Distances on dual-weighted directed graphs using priority-queue shortest
#' paths. Weighted directed graphs have weights from A to B which may differ
#' from those from B to A. Dual-weighted directed graphs have two sets of such
#' weights. A canonical example is a street network to be used for routing in
#' which routes are calculated by weighting distances according to the type of
#' way and mode of transport, yet lengths of routes must be calculated from
#' direct distances.
#'
#' @section The Main Function:
#' \itemize{
#' \item \code{\link{dodgr_dists}}: Calculate pair-wise distances between
#' specified pairs of points in a graph.
#' }
#'
#' @section Functions to Obtain Graphs:
#' \itemize{
#' \item \code{\link{dodgr_streetnet}}: Extract a street network in Simple
#' Features (\code{sf}) form.
#' \item \code{\link{weight_streetnet}}: Convert an \code{sf}-formatted street
#' network to a \code{dodgr} graph through applying specified weights to all
#' edges.
#' }
#'
#' @section Functions to Modify Graphs:
#' \itemize{
#' \item \code{\link{dodgr_components}}: Number all graph edges according to
#' their presence in distinct connected components.
#' \item \code{\link{dodgr_convert_graph}}: Convert a graph of arbitrary form
#' into a standardised, minimal form for submission to \code{dodgr} routines.
#' \item \code{\link{dodgr_contract_graph}}: Contract a graph by removing
#' redundant edges.
#' }
#'
#' @section Miscellaneous Functions:
#' \itemize{
#' \item \code{\link{dodgr_sample}}: Randomly sample a graph, returning a single
#' connected component of a defined number of vertices.
#' \item \code{\link{dodgr_vertices}}: Extract all vertices of a graph.
#' \item \code{\link{compare_heaps}}: Compare the performance of different
#' priority queue heap structures for a given type of graph.
#' }
#'
#' @name dodgr
#' @docType package
#' @importFrom igraph distances E make_directed_graph
#' @importFrom magrittr %>%
#' @importFrom methods is
#' @importFrom osmdata add_osm_feature getbb opq osmdata_sf
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

#' hampi
#'
#' A sample street network from the township of Hampi, Karnataka, India.
#'
#' @name hampi
#' @docType data
#' @keywords datasets
#' @format A Simple Features \code{sf} \code{data.frame} containing the street
#' network of Hampi.
#'
#' @note Can be re-created with the following command, which also removes 
#' extraneous columns to reduce size:
#' @examples \dontrun{
#' hampi <- dodgr_streetnet("hampi india")
#' cols <- c ("osm_id", "highway", "oneway", "geometry")
#' hampi <- hampi [, which (names (hampi) %in% cols)]
#' }
NULL
