#' compare_heaps
#'
#' Perform timing comparison between different kinds of heaps as well as with
#' equivalent \code{igraph} routine \code{distances}.
#'
#' @param graph \code{data.frame} object representing the network graph
#' @param replications Number of replications to be used in comparison
#' @return Result of \code{rbenachmar::benchmark} comparison in
#' \code{data.frame} form.
#'
#' @export
compare_heaps <- function(graph, replications = 10)
{
    # set up igraph:
    edges <- cbind (graph$from_id, graph$to_id)
    nodes <- unique (as.vector (edges)) # used below in test comparison
    edges <- as.vector (t (edges))
    igr <- igraph::make_directed_graph (edges)
    igraph::E (igr)$weight <- graph$d_weighted

    rbenchmark::benchmark (
                           d <- test (graph, heap = "BHeap"),
                           d <- test (graph, heap = "FHeap"),
                           d <- test (graph, heap = "TriHeap"),
                           d <- test (graph, heap = "TriHeapExt"),
                           d <- test (graph, heap = "Heap23"),
                           igraph::distances (igr, v = nodes, to = nodes,
                                              mode = "out"),
                           replications = 10, order = "relative")
}
