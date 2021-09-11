#' Proportional distances along different edge categories
#'
#' @inheritParams dodgr_dists
#' @param graph `data.frame` or equivalent object representing the network
#' graph which must have a column named "edge_type" with integer types along
#' which proportional distances are to be aggregated. The value of "edge_type" =
#' 0 is the default type for which proportional distances are not aggregated.
#' @return An expanded distance matrix with number of rows equal to number of
#' "from" points, and number of columns equal to number of "to" points
#' \emph{times} the number of distinct edge types.
#' @export
dodgr_dists_proportional <- function (graph,
                                      from = NULL,
                                      to = NULL,
                                      heap = "BHeap",
                                      quiet = TRUE) {

    if (!"edge_type" %in% names (graph))
        stop ("graph must have a column named 'edge_type'")
    if (!is.integer (graph$edge_type))
        stop ("'edge_type' column must contain integer values only")
    if (any (graph$edge_type < 0L))
        stop ("'edge_type' must contain non-negative values only")

    edge_type <- graph$edge_type

    graph <- tbl_to_df (graph)

    hps <- get_heap (heap, graph)
    heap <- hps$heap
    graph <- hps$graph

    gr_cols <- dodgr_graph_cols (graph)
    if (is.na (gr_cols$from) | is.na (gr_cols$to)) {
        scols <- find_spatial_cols (graph)
        graph$from_id <- scols$xy_id$xy_fr_id
        graph$to_id <- scols$xy_id$xy_to_id
        gr_cols <- dodgr_graph_cols (graph)
    }
    is_spatial <- is_graph_spatial (graph)
    if (!is_spatial)
        stop ("proportional distances only implemented for spatial graphs")

    vert_map <- make_vert_map (graph, gr_cols, is_spatial)

    # adjust to/from for turn penalty where that exists:
    from <- to_from_with_tp (graph, from, from = TRUE)
    to <- to_from_with_tp (graph, to, from = FALSE)

    from_index <- get_to_from_index (graph, vert_map, gr_cols, from)
    to_index <- get_to_from_index (graph, vert_map, gr_cols, to)

    graph <- convert_graph (graph, gr_cols)
    graph$edge_type <- edge_type

    if (!quiet)
        message ("Calculating shortest paths ... ", appendLF = FALSE)

    d <- rcpp_get_sp_dists_proportional (graph,
                                         vert_map,
                                         from_index$index,
                                         to_index$index,
                                         heap)

    return (d)
}
