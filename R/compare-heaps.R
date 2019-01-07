#' compare_heaps
#'
#' Perform timing comparison between different kinds of heaps as well as with
#' equivalent `igraph` routine `distances`. To do this, a random
#' sub-graph containing a defined number of vertices is first selected.
#' Alternatively, this random sub-graph can be pre-generated with the
#' `dodgr_sample` function and passed directly.
#'
#' @param graph `data.frame` object representing the network graph (or a
#' sub-sample selected with code{dodgr_sample})
#' @param nverts Number of vertices used to generate random sub-graph. If a
#' non-numeric value is given, the whole graph will be used.
#' @param replications Number of replications to be used in comparison
#' @return Result of `rbenachmar::benchmark` comparison in
#' `data.frame` form.
#'
#' @note \pkg{igraph} caches intermediate results of graph processing, so
#' the \pkg{igraph} comparisons will be faster on subsequent runs. To obtain
#' fair comparisons, run only once or re-start the current R session.
#'
#' @export
#' @examples
#' graph <- weight_streetnet (hampi)
#' compare_heaps (graph, nverts = 100, replications = 1)
compare_heaps <- function(graph, nverts = 100, replications = 2)
{
    if (is.numeric (nverts))
        graph <- dodgr_sample (graph, nverts = nverts)
    graph_contracted <- dodgr_contract_graph (graph)$graph

    # route only between points on the contracted graph:
    gr_cols <- dodgr_graph_cols (graph)
    # gr_cols are (edge_id, from, to, d, w, component, xfr, yfr, xto, yto)
    from_id <- unique (graph_contracted [[gr_cols [2] ]])
    to_id <- unique (graph_contracted [[gr_cols [3] ]])

    igr <- dodgr_to_igraph (graph)

    rbenchmark::benchmark (
                           d <- dodgr_dists (graph, from = from_id, to = to_id,
                                             heap = "BHeap"),
                           d <- dodgr_dists (graph, from = from_id, to = to_id,
                                             heap = "FHeap"),
                           d <- dodgr_dists (graph, from = from_id, to = to_id,
                                             heap = "TriHeap"),
                           d <- dodgr_dists (graph, from = from_id, to = to_id,
                                             heap = "TriHeapExt"),
                           d <- dodgr_dists (graph, from = from_id, to = to_id,
                                             heap = "Heap23"),
                           d <- dodgr_dists (graph_contracted, from = from_id,
                                             to = to_id, heap = "BHeap"),
                           d <- dodgr_dists (graph_contracted, from = from_id,
                                             to = to_id, heap = "FHeap"),
                           d <- dodgr_dists (graph_contracted, from = from_id,
                                             to = to_id, heap = "TriHeap"),
                           d <- dodgr_dists (graph_contracted, from = from_id,
                                             to = to_id, heap = "TriHeapExt"),
                           d <- dodgr_dists (graph_contracted, from = from_id,
                                             to = to_id, heap = "Heap23"),
                           d <- dodgr_dists (graph_contracted, from = from_id,
                                             to = to_id, heap = "set"),
                           d <- igraph::distances (igr, v = from_id, to = to_id,
                                              mode = "out"),
                           replications = 10, order = "relative")
}
