#' test
#'
#' @param graph \code{data.frame} object representing the network graph
#' @param heap Type of heap to use in priority queue. Options include
#' Fibonacci Heap (default; \code{FHeap}), Binary Heap (\code{BHeap}),
#' \code{Radix}, Trinomial Heap (\code{Tri}), Extended Trinomial Heap
#' (\code{TriExt}, and 2-3 Heap (\code{Heap23}).
#' @return square matrix of distances between nodes
#'
#' @export
test <- function(graph, heap = 'FHeap') {
    heaps <- c ("FHeap", "BHeap", "Radix", "Tri", "TriExt", "Heap23")
    heap <- match.arg (arg = heap, choices = heaps)

    d <- rcpp_get_sp (graph, heap)
    if (max (d) > 1e30) # float max ~ 1e38
        d [d == max (d)] <- NA

    return (d)
}
