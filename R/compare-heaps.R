#' compare_heaps
#'
#' Perform timing comparison between different kinds of heaps as well as with
#' equivalent \code{igraph} routine \code{distances}. To do this, a random
#' sub-graph containing a defined number of vertices is first selected.
#' Alternatively, this random sub-graph can be pre-generated with the
#' \code{dodgr_sample} function and passed directly.
#'
#' @param graph \code{data.frame} object representing the network graph (or a
#' sub-sample selected with code{dodgr_sample})
#' @param nverts Number of vertices used to generate random sub-graph
#' @param replications Number of replications to be used in comparison
#' @return Result of \code{rbenachmar::benchmark} comparison in
#' \code{data.frame} form.
#'
#' @export
compare_heaps <- function(graph, nverts = 1000, replications = 10)
{
    graph <- dodgr_sample (graph, nverts = nverts)

    # set up igraph:
    fr_col <- find_fr_id_col (graph)
    to_col <- find_to_id_col (graph)
    edges <- cbind (graph [, fr_col], graph [, to_col])
    nodes <- unique (as.vector (edges)) # used below in test comparison
    edges <- as.vector (t (edges))
    igr <- igraph::make_directed_graph (edges)
    igraph::E (igr)$weight <- graph [, find_d_col (graph)]

    rbenchmark::benchmark (
                           d <- dodgr_dists (graph, heap = "BHeap"),
                           d <- dodgr_dists (graph, heap = "FHeap"),
                           d <- dodgr_dists (graph, heap = "TriHeap"),
                           d <- dodgr_dists (graph, heap = "TriHeapExt"),
                           d <- dodgr_dists (graph, heap = "Heap23"),
                           d <- igraph::distances (igr, v = nodes, to = nodes,
                                              mode = "out"),
                           replications = 10, order = "relative")
}
