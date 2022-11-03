# functions to process 'from' and 'to' parameters of distance, iso, and flow
# functions, mostly via the `to_from_index_with_tp` function, which
# pre-processes the `from` and `to` arguments, returning both indices into the
# vertex map, and corresponding ID values where those exist.
#
# ------------------------------------------------------------------------


#' Get lists of 'from' and 'to' indices, potentially mapped on to compound
#' junctions for graphs with turn penalties.
#'
#' This function calls the following `fill_to_from_index` to generate the actual
#' values. For graphs with turn penalties, it also returns the modified version
#' of the graph including compound junctions.
#' @noRd
to_from_index_with_tp <- function (graph, from, to) {

    gr_cols <- dodgr_graph_cols (graph)
    is_spatial <- is_graph_spatial (graph)
    vert_map <- make_vert_map (graph, gr_cols, is_spatial)

    from_index <-
        fill_to_from_index (graph, vert_map, gr_cols, from, from = TRUE)
    to_index <- fill_to_from_index (graph, vert_map, gr_cols, to, from = FALSE)

    compound <- (get_turn_penalty (graph) > 0.0)
    graph_compound <- compound_junction_map <- NULL
    if (compound) {
        if (methods::is (graph, "dodgr_contracted")) {
            warning (
                "graphs with turn penalties should be submitted in full, ",
                "not contracted form;\nsubmitting contracted graphs may ",
                "produce unexpected behaviour."
            )
        }
        res <- create_compound_junctions (graph)
        graph_compound <- res$graph
        compound_junction_map <- res$edge_map

        # remap any 'from' and 'to' vertices to compound junction versions:
        vert_map <- make_vert_map (graph_compound, gr_cols, is_spatial)

        from_index <- remap_tf_index_for_tp (from_index, vert_map, from = TRUE)
        to_index <- remap_tf_index_for_tp (to_index, vert_map, from = FALSE)
    }

    return (list (
        from = from_index, to = to_index, vert_map = vert_map,
        compound = compound, graph_compound = graph_compound,
        compound_junction_map = compound_junction_map
    ))
}

#' fill_to_from_index
#'
#' @noRd
fill_to_from_index <- function (graph,
                                vert_map,
                                gr_cols,
                                pts,
                                from = TRUE) {

    id <- NULL
    if (is.null (pts)) {
        index <- seq_len (nrow (vert_map)) - 1L
        if (!is.null (vert_map$vert)) {
            id <- vert_map$vert
        }
    } else {
        index_id <- get_index_id_cols (graph, gr_cols, vert_map, pts)
        if (any (is.na (index_id$id))) {
            stop ("Unable to match all routing points to graph vertices")
        }
        index <- index_id$index - 1L # 0-based
        id <- index_id$id
    }

    if (!is.null (id)) {

        tflab <- ifelse (from, "from", "to")

        i <- which (id %in% graph [[gr_cols [[tflab]]]])
        if (length (i) == 0) {
            i <- which (vert_map$vert [index] %in% graph [[gr_cols [[tflab]]]])
        }

        index <- index [i]
        id <- id [i]
    }
    list (index = index, id = id)
}

#' Remap 'from_index' and 'to_index' values on to the compound junctions present
#' in 'vert_map'.
#'
#' @param index Either 'from_index' or 'to_index' calculated
#' @noRd
remap_tf_index_for_tp <- function (index, vert_map, from = TRUE) {

    vert_index <- match (index$id, vert_map$vert)
    if (from) {
        no_start <- which (!grepl ("\\_start$", index$id))
        vert_index_id <- index$id
        vert_index_id [no_start] <- paste0 (index$id [no_start], "_start")
    } else {
        no_end <- which (!grepl ("\\_end$", index$id))
        vert_index_id <- index$id
        vert_index_id [no_end] <- paste0 (index$id [no_end], "_end")
    }
    vert_index_comp <- match (vert_index_id, vert_map$vert)
    na_index <- which (!is.na (vert_index_comp))
    vert_index [na_index] <- vert_index_comp [na_index]

    index$index <- vert_index - 1L # zero-based
    index$id [na_index] <- vert_index_id [na_index]

    return (index)
}
