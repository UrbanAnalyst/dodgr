#' Proportional distances along different edge categories
#'
#' @inheritParams dodgr_dists
#' @param graph `data.frame` or equivalent object representing the network
#' graph which must have a column named "edge_type" with integer types along
#' which proportional distances are to be aggregated. The value of "edge_type" =
#' 0 is the default type for which proportional distances are not aggregated.
#' @return A list of distance matrices of equal dimensions (length(from),
#' length(to)), the first of which ("distance") holds the final distances, while
#' the rest are one matrix for each unique value of "edge_type", holding the
#' distances traversed along those types of edges only.
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
    et <- sort (unique (edge_type))

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

    n <-length (to)
    d0 <- list ("distances" = d [, seq (n)])
    d <- lapply (et, function (i) {
                     index <- i * n + seq (n) - 1
                     d [, index]    })
    names (d) <- et

    return (c (d0, d))
}
