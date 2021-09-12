#' Proportional distances along different edge categories
#'
#' @inheritParams dodgr_dists
#' @param graph `data.frame` or equivalent object representing the network
#' graph which must have a column named "edge_type" which labels categories of
#' edge types along which proportional distances are to be aggregated (see
#' Note).
#' @return A list of distance matrices of equal dimensions (length(from),
#' length(to)), the first of which ("distance") holds the final distances, while
#' the rest are one matrix for each unique value of "edge_type", holding the
#' distances traversed along those types of edges only.
#'
#' @note The "edge_type" column in the graph can contain any kind of discrete or
#' categorical values, although integer values of 0 are not permissible. `NA`
#' values are ignored.
#' @export
dodgr_dists_proportional <- function (graph,
                                      from = NULL,
                                      to = NULL,
                                      heap = "BHeap",
                                      quiet = TRUE) {

    if (!"edge_type" %in% names (graph))
        stop ("graph must have a column named 'edge_type'")
    if (is.integer (graph$edge_type) & any (graph$edge_type == 0L))
        stop ("graphs with integer edge_type columns may not contain 0s")

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
    edge_type_table <- table (edge_type)
    graph$edge_type <- match (edge_type, names (edge_type_table))
    graph$edge_type [is.na (graph$edge_type)] <- 0L

    if (!quiet)
        message ("Calculating shortest paths ... ", appendLF = FALSE)

    d <- rcpp_get_sp_dists_proportional (graph,
                                         vert_map,
                                         from_index$index,
                                         to_index$index,
                                         heap)

    n <-length (to)
    d0 <- list ("distances" = d [, seq (n)])
    d <- lapply (seq_along (edge_type_table), function (i) {
                     index <- i * n + seq (n) - 1
                     d [, index]    })
    names (d) <- names (edge_type_table)

    res <- c (d0, d)
    class (res) <- c (class (res), "dodgr_dists_proportional")

    return (res)
}


#' Transform a result from 'dodgr_dists_proportional' to summary statistics
#'
#' @param d A 'dodgr_dists_proportional' object
#' @return The summary statistics (invisibly)
#' @export
summary.dodgr_dists_proportional <- function (d) {

    edge_types <- names (d) [-1]

    d0 <- d$distances
    sum_d0 <- sum (d0, na.rm = TRUE)
    d <- d [-1]

    dprop <- vapply (d, function (i) sum (i, na.rm = TRUE) / sum_d0,
                     numeric (1))

    message ("Proportional distances along each kind of edge:")
    for (i in seq_along (dprop))
        message ("  ", names (dprop) [i],
                 ": ", round (dprop [i], digits = 4))

    invisible (dprop)
}
