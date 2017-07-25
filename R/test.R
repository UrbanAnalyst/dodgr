#' test
#'
#' @param graph \code{data.frame} object representing the network graph
#' @param quiet If FALSE, display progress information
#' @param pq Type of priority queue to use. Options include Fibonacci Heap
#' (default; \code{FHeap}), \code{Radix}, \code{Tri}, and \code{Heap23}.
#' @return square matrix of distances between nodes
#'
#' @export
test <- function(graph, pq = 'FHeap') {
    pq <- tolower (pq)
    if (!pq %in% c ("fheap", "head23", "radix", "tri", "triext"))
        stop (paste0 ("priority queue type not recognised; ",
                      "must be one of (fheap, head23, radix, tri, triext"))

    d <- rcpp_get_sp (graph, tolower (pq))
    if (max (d) > 1e30) # float max ~ 1e38
        d [d == max (d)] <- NA

    return (d)
}
