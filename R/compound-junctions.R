#' Create compound junctions representing turn penalties and restrictions
#'
#' @param graph A weighted streetnet with a non-zero "turn_penalty" attribute,
#' and other parameters stored as graph attributes during streetnet
#' construction.
#' @return List of two items: The expanded graph included compound junctions,
#' and an edge map mapping the new compound edges back on to original graph
#' edges.
#' @noRd
create_compound_junctions <- function (graph) {

    left_side <- attr (graph, "left_side")
    wt_profile <- attr (graph, "wt_profile")
    wt_profile_file <- attr (graph, "wt_profile_file")

    res <- join_junctions_to_graph (
        graph, wt_profile, wt_profile_file,
        left_side
    )
    if (are_turns_restricted (wt_profile, wt_profile_file)) {
        res <- remove_turn_restrictions (graph, res)
    }

    return (res)
}

get_turn_penalty <- function (graph) {
    tp <- attr (graph, "turn_penalty")
    if (!is.numeric (tp)) {
        tp <- 0.0
    }
    return (tp)
}

join_junctions_to_graph <- function (graph, wt_profile, wt_profile_file,
                                     left_side = FALSE) {

    turn_penalty <- attr (graph, "turn_penalty")
    if (is.null (turn_penalty) & !is.null (wt_profile_file)) {
        turn_penalty <- get_turn_penalties (wt_profile, wt_profile_file)$turn
    }

    resbind <- edge_map <- NULL

    if (turn_penalty > 0) {

        res <- rcpp_route_times (graph, left_side, turn_penalty)
        edge_map <- data.frame (
            "edge" = res$graph$edge_,
            "e_in" = res$graph$old_edge_in,
            "e_out" = res$graph$old_edge_out,
            stringsAsFactors = FALSE
        )
        res$graph$old_edge_in <- res$graph$old_edge_out <- NULL

        index <- which (graph$.vx0 %in% res$junction_vertices)
        graph$.vx0 [index] <- paste0 (graph$.vx0 [index], "_start")
        index <- which (graph$.vx1 %in% res$junction_vertices)
        graph$.vx1 [index] <- paste0 (graph$.vx1 [index], "_end")

        # pad out extra columns of res to match any extra in original graph
        resbind <- data.frame (array (NA, dim = c (
            nrow (res$graph),
            ncol (graph)
        )))
        names (resbind) <- names (graph)
        resbind [, which (names (graph) %in% names (res$graph))] <- res$graph
        if ("object_" %in% names (graph)) {
            # SC graphs:
            index <- match (resbind$object_, graph$object_)
            index_col <- which (!names (graph) %in% names (res$graph))
            resbind [, index_col] <- graph [index, index_col]
        }
        graph <- rbind (graph, resbind)
    }
    list (graph = graph, edge_map = edge_map)
}

#' Extract data on turn restictions
#'
#' This function is called during initial graph construction in
#' `weight_streetnet`, and the data appended as attributes to the graph.
#'
#' @return A list of 2 `data.frame` objects, on containing strict turn
#' restrictions (turn = "no"), and the only containing exclusive turn
#' permissions ("only" = "left", for example).
#' @noRd
extract_turn_restrictions <- function (x) {

    rels <- x$relation_properties # x from restrictions query above!!
    restriction_rels <- rels [rels$key == "restriction", ]
    index <- which (x$relation_members$relation_ %in%
        restriction_rels$relation_)
    restriction_ways <- x$relation_members [index, ]

    rr_no <- restriction_rels [grep ("^no\\_", restriction_rels$value), ]
    rr_only <- restriction_rels [grep ("^only\\_", restriction_rels$value), ]
    rw_no <- restriction_ways [restriction_ways$relation_ %in%
        rr_no$relation_, ]
    rw_only <- restriction_ways [restriction_ways$relation_ %in%
        rr_only$relation_, ]

    r_to_df <- function (r) {

        r0 <- data.frame (matrix (nrow = 0L, ncol = 4L))
        names (r0) <- c ("relation", "node", "from", "to")
        if (is.null (r)) {
            return (r0)
        } else if (nrow (r) == 0L) {
            return (r0)
        }

        r <- lapply (
            split (r, f = factor (r$relation_)),
            function (i) {
                # Relations should have 1 node type, and 2 way types with roles
                # "from" and "to":
                node <- which (i$type == "node") [1L]
                index <- match (c ("from", "to"), i$role) [1:2]
                c (
                    i$relation_ [1],
                    i$member [node],
                    i$member [index]
                )
            }
        )
        r <- data.frame (do.call (rbind, r))
        names (r) <- c ("relation", "node", "from", "to")
        return (stats::na.omit (r))
    }

    return (list (
        rw_no = r_to_df (rw_no),
        rw_only = r_to_df (rw_only)
    ))
}

#' Remove turn restrictions and replace with compound junctions.
#'
#' @param x The original `sc` object which ,when generated from
#' `dodgr_streetnet_sc`, includes turn restriction data
#' @param graph The processed by not yet turn-contracted graph
#' @param res The result of `join_junctions_to_graph`, with turn-contracted
#' `graph` and `edge_map` components.
#' @noRd
remove_turn_restrictions <- function (graph, res) {

    # These are the attributes inserted in initial streeetnet construction via
    # the previous `extract_turn_restrictions` function.
    rw_no <- attr (graph, "turn_restrictions_no")
    rw_only <- attr (graph, "turn_restrictions_only")

    # rw_no <- dplyr::filter (rw_no, node == node_id)

    # indices of edges which have to be removed:
    index_in <- which (
        graph$object_ %in% rw_no$from & graph$.vx1 %in% rw_no$node
    )
    index_out <- which (
        graph$object_ %in% rw_no$to & graph$.vx0 %in% rw_no$node
    )
    graph_in <- graph [index_in, ]
    graph_out <- graph [index_out, ]
    edges_to_remove <- unique (graph$edge_ [c (index_in, index_out)])

    # Then construct new compound edges:
    index_out_other <- which (
        graph$.vx0 %in% rw_no$node & !graph$edge_ %in% graph_out$edge_
    )
    graph_out_compound <- graph [index_out_other, ]
    graph_in <- graph_in [which (graph_in$.vx1 %in% graph_out_compound$.vx0), ]

    index <- match (graph_out_compound$.vx0, graph_in$.vx1)
    graph_out_compound$.vx0 <- graph_in$.vx0 [index]
    graph_out_compound$object_ <-
        paste0 (graph_in$object_ [index], "_", graph_out_compound$object_)
    graph_out_compound$d <- graph_in$d [index] + graph_out_compound$d
    graph_out_compound$d_weighted <-
        graph_in$d_weighted [index] + graph_out_compound$d_weighted
    graph_out_compound$time <- graph_in$time [index] + graph_out_compound$time
    graph_out_compound$time_weighted <-
        graph_in$time_weighted [index] + graph_out_compound$time_weighted

    # Extend edge map:
    res$edge_map <- rbind (
        res$edge_map,
        data.frame (
            edge = paste0 ("j_", rcpp_gen_hash (length (index), 10)),
            e_in = graph_in$edge_ [index],
            e_out = graph_out_compound$edge_
        )
    )

    # Remove original edges from graph, and add new compound ones:
    res$graph <- rbind (
        res$graph [which (!res$graph$edge_ %in% edges_to_remove), ],
        graph_out_compound
    )

    return (res)
}


#' Aggregate edge values across compound junctions back on to original edges.
#'
#' This function is only called in the flow aggregation routines, which differ
#' in their handling of compound junctions. For flow routines, graphs are simply
#' created with compound junctions added, and subsequent flows re-mapped using
#' this function.
#'
#' @noRd
uncompound_junctions <- function (graph, new_cols,
                                  compound_junction_map = NULL) {

    tp <- get_turn_penalty (graph)

    if (is (graph, "dodgr_streetnet_sc") && tp > 0 &&
        !is.null (compound_junction_map)) {

        jm <- compound_junction_map

        gr_cols <- dodgr_graph_cols (graph)
        index_junction <- match (jm$edge, graph [[gr_cols$edge_id]])
        index_edge_in <- match (jm$e_in, graph [[gr_cols$edge_id]])
        index_edge_out <- match (jm$e_out, graph [[gr_cols$edge_id]])

        for (n in new_cols) {

            graph [[n]] [index_edge_in] <-
                graph [[n]] [index_edge_in] +
                graph [[n]] [index_junction]

            graph [[n]] [index_edge_out] <-
                graph [[n]] [index_edge_out] +
                graph [[n]] [index_junction]
        }

        # next line removes all the compound turn angle edges:
        graph <- graph [which (!graph [[gr_cols$edge_id]] %in% jm$edge), ]
        graph$.vx0 <- gsub ("_start$", "", graph$.vx0)
        graph$.vx1 <- gsub ("_end$", "", graph$.vx1)
    }

    return (graph)
}
