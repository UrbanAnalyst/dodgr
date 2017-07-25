#' test
#'
#' @param graph \code{data.frame} object representing the network graph
#' @param quiet If FALSE, display progress information
#' @return square matrix of distances between nodes
#'
#' @export
test <- function(graph, quiet = TRUE) {
    d <- rcpp_get_sp (graph, quiet)
    if (max (d) > 1e30) # float max ~ 1e38
        d [d == max (d)] <- NA

    return (d)
}
