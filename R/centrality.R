#' dodgr_centrality
#'
#' Calculate betweenness centrality for a `dodgr` network, in either vertex- or
#' edge-based form.
#'
#' @param graph `data.frame` or equivalent object representing the network
#' graph (see Details)
#' @param contract If `TRUE`, centrality is calculated on contracted graph
#' before mapping back on to the original full graph. Note that for street
#' networks, in particular those obtained from the \pkg{osmdata} package, vertex
#' placement is effectively arbitrary except at junctions; centrality for such
#' graphs should only be calculated between the latter points, and thus
#' `contract` should always be `TRUE`.
#' @param edges If `TRUE`, centrality is calculated for graph edges, returning
#' the input `graph` with an additional `centrality` column; otherwise
#' centrality is calculated for vertices, returning the equivalent of
#' `dodgr_vertices(graph)`, with an additional vertex-based `centrality` column.
#' @param heap Type of heap to use in priority queue. Options include
#' Fibonacci Heap (default; `FHeap`), Binary Heap (`BHeap`),
#' `Radix`, Trinomial Heap (`TriHeap`), Extended Trinomial Heap
#' (`TriHeapExt`, and 2-3 Heap (`Heap23`).
#' @return Modified version of graph with additonal `centrality` column added.
#'
#' @export
dodgr_centrality <- function (graph, contract = TRUE, edges = TRUE, heap = "BHeap")
{
    if ("centrality" %in% names (graph))
        warning ("graph already has a 'centrality' column; ",
                  "this will be overwritten")

    hps <- get_heap (heap, graph)
    heap <- hps$heap
    graph <- hps$graph

    gr_cols <- dodgr_graph_cols (graph)

    if (contract & !methods::is (graph, "dodgr_contracted"))
    {
        graph_full <- graph
        graph <- dodgr_contract_graph (graph)
        hashc <- get_hash (graph, hash = FALSE)
        fname_c <- file.path (tempdir (),
                              paste0 ("dodgr_edge_map_", hashc, ".Rds"))
        if (!file.exists (fname_c))
            stop ("something went wrong extracting the edge_map ... ") # nocov
        edge_map <- readRDS (fname_c)
    }

    vert_map <- make_vert_map (graph, gr_cols)

    graph2 <- convert_graph (graph, gr_cols)

    centrality <- rcpp_centrality (graph2, vert_map, "BHeap", edges)

    if (edges)
    {
        graph$centrality <- centrality
        if (contract)
            graph <- uncontract_graph (graph, edge_map, graph_full)
        res <- graph
    } else
    {
        res <- dodgr_vertices (graph)
        res$centrality <- centrality
    }

    return (res)
}
