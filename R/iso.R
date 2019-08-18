#' dodgr_isodists
#'
#' Isodistances from a specified point
#'
#' @param graph `data.frame` or equivalent object representing the network
#' graph (see Notes)
#' @param from Vector or matrix of points **from** which route distances are to
#' be calculated (see Notes)
#' @param dlim Desired limit of isodistance
#' @param heap Type of heap to use in priority queue. Options include
#' Fibonacci Heap (default; `FHeap`), Binary Heap (`BHeap`),
#' `Radix`, Trinomial Heap (`TriHeap`), Extended Trinomial Heap
#' (`TriHeapExt`, and 2-3 Heap (`Heap23`).
#' @return Isodistance points
#'
#' @export 
#' @examples
#' graph <- weight_streetnet (hampi)
#' from <- sample (graph$from_id, size = 100)
#' to <- sample (graph$to_id, size = 50)
#' d <- dodgr_dists (graph, from = from, to = to)
#' # d is a 100-by-50 matrix of distances between `from` and `to`
dodgr_isodists <- function (graph, from = NULL, dlim = NULL, heap = 'BHeap')
{
    if (is.null (dlim))
        stop ("dlim must be specified")
    if (!is.numeric (dlim) | length (dlim) > 1)
        stop ("dlim must be a single numeric value")

    graph <- tbl_to_df (graph)

    hps <- get_heap (heap, graph)
    heap <- hps$heap
    graph <- hps$graph

    gr_cols <- dodgr_graph_cols (graph)
    if (is.na (gr_cols$from) | is.na (gr_cols$to))
    {
        scols <- find_spatial_cols (graph)
        graph$from_id <- scols$xy_id$xy_fr_id
        graph$to_id <- scols$xy_id$xy_to_id
        gr_cols <- dodgr_graph_cols (graph)
    }
    vert_map <- make_vert_map (graph, gr_cols, FALSE)

    tp <- attr (graph, "turn_penalty")
    tp <- ifelse (is.null (tp), 0, tp)
    if (is (graph, "dodgr_streetnet_sc") & tp > 0)
    {
        if (!is.null (from))
        {
            from <- nodes_arg_to_pts (from, graph)
            from <- remap_verts_with_turn_penalty (graph, from, from = TRUE)
        }
    }

    from_index <- get_to_from_index (graph, vert_map, gr_cols, from)

    graph <- convert_graph (graph, gr_cols)

    d <- rcpp_get_iso (graph, vert_map, from_index$index, dlim, heap)

    if (!is.null (from_index$id))
        rownames (d) <- from_index$id
    else
        rownames (d) <- vert_map$vert
    colnames (d) <- vert_map$vert

    return (d)
}
