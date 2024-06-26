#' Compare timings of different sort heaps for a given input graph.
#'
#' Perform timing comparison between different kinds of heaps as well as with
#' equivalent routines from the \pkg{igraph} package. To do this, a random
#' sub-graph containing a defined number of vertices is first selected.
#' Alternatively, this random sub-graph can be pre-generated with the
#' `dodgr_sample` function and passed directly.
#'
#' @param graph `data.frame` object representing the network graph (or a
#' sub-sample selected with `dodgr_sample`)
#' @param nverts Number of vertices used to generate random sub-graph. If a
#' non-numeric value is given, the whole graph will be used.
#' @param replications Number of replications to be used in comparison
#' @return Result of `bench::mark` comparison.
#'
#' @family misc
#' @export
#' @examples
#' graph <- weight_streetnet (hampi)
#' \dontrun{
#' compare_heaps (graph, nverts = 1000, replications = 1)
#' }
compare_heaps <- function (graph,
                           nverts = 100,
                           replications = 2) {

    requireNamespace ("bench")
    requireNamespace ("igraph")

    if (is.numeric (nverts)) {
        graph <- dodgr_sample (graph, nverts = nverts)
    }
    graph_contracted <- dodgr_contract_graph (graph)

    # route only between points on the contracted graph:
    gr_cols <- dodgr_graph_cols (graph)
    from_id <- unique (graph_contracted [[gr_cols$from]])
    to_id <- unique (graph_contracted [[gr_cols$to]])

    igr <- dodgr_to_igraph (graph)

    bench::mark (
        BHeap = dodgr_dists (
            graph,
            from = from_id,
            to = to_id,
            heap = "BHeap"
        ),
        FHeap = dodgr_dists (
            graph,
            from = from_id,
            to = to_id,
            heap = "FHeap"
        ),
        TriHeap = dodgr_dists (
            graph,
            from = from_id,
            to = to_id,
            heap = "TriHeap"
        ),
        TriHeapExt = dodgr_dists (
            graph,
            from = from_id,
            to = to_id,
            heap = "TriHeapExt"
        ),
        Heap23 = dodgr_dists (
            graph,
            from = from_id,
            to = to_id,
            heap = "Heap23"
        ),
        BHeap_contracted = dodgr_dists (
            graph_contracted,
            from = from_id,
            to = to_id,
            heap = "BHeap"
        ),
        FHeap_contracted = dodgr_dists (
            graph_contracted,
            from = from_id,
            to = to_id,
            heap = "FHeap"
        ),
        TriHeap_contracted = dodgr_dists (
            graph_contracted,
            from = from_id,
            to = to_id,
            heap = "TriHeap"
        ),
        TriHeapExt_contracted = dodgr_dists (
            graph_contracted,
            from = from_id,
            to = to_id,
            heap = "TriHeapExt"
        ),
        Heap23_contracted = dodgr_dists (
            graph_contracted,
            from = from_id,
            to = to_id,
            heap = "Heap23"
        ),
        igraph = igraph::distances (
            igr,
            v = from_id,
            to = to_id,
            mode = "out"
        ),
        check = FALSE # contracted don't necessarily equal full dists here
    )
}
