#' test
#'
#' @param graph \code{data.frame} object representing the network graph
#' @param heap Type of heap to use in priority queue. Options include
#' Fibonacci Heap (default; \code{FHeap}), Binary Heap (\code{BHeap}),
#' \code{Radix}, Trinomial Heap (\code{TriHeap}), Extended Trinomial Heap
#' (\code{TriHeapExt}, and 2-3 Heap (\code{Heap23}).
#' @return square matrix of distances between nodes
#' @export
test <- function(graph, heap = 'FHeap') {
    heaps <- c ("FHeap", "BHeap", "Radix", "TriHeap", "TriHeapExt", "Heap23")
    heap <- match.arg (arg = heap, choices = heaps)
    if (heap == "Radix")
    {
        dfr <- min (abs (c (graph$d %% 1, graph$d %% 1 - 1)))
        if (dfr > 1e-6)
        {
            message (paste0 ("RadixHeap can only be implemented for ",
                             "integer weights;\nall weights will now be ",
                             "rounded"))
            graph$d <- round (graph$d)
            graph$d_weighted <- round (graph$d_weighted)
        }
    }
    d <- rcpp_get_sp (graph, heap)

    if (max (d) > 1e30) # float max ~ 1e38
        d [d == max (d)] <- NA

    return (d)
}
