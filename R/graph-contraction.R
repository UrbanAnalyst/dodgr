#' Contract graph to junction vertices only.
#'
#' Removes redundant (straight-line) vertices from graph, leaving only junction
#' vertices.
#'
#' @param graph A flat table of graph edges. Must contain columns labelled
#' `from` and `to`, or `start` and `stop`. May also contain
#' similarly labelled columns of spatial coordinates (for example
#' `from_x`) or `stop_lon`).
#' @param verts Optional list of vertices to be retained as routing points.
#' These must match the `from` and `to` columns of `graph`.
#' @param nocache If `FALSE` (default), load cached version of contracted graph
#' if previously calculated and cached. If `TRUE`, then re-contract graph even
#' if previously calculated version has been stored in cache.
#'
#' @return A contracted version of the original `graph`, containing the same
#' number of columns, but with each row representing an edge between two
#' junction vertices (or between the submitted `verts`, which may or may not be
#' junctions).
#' @family modification
#' @export
#' @examples
#' graph <- weight_streetnet (hampi)
#' nrow (graph) # 5,973
#' graph <- dodgr_contract_graph (graph)
#' nrow (graph) # 662
dodgr_contract_graph <- function (graph, verts = NULL, nocache = FALSE) {

    if (nrow (graph) == 0) {
        stop ("graph is empty")
    } # nocov

    # px is the R6 processx object for initial caching
    px <- NULL
    if ("px" %in% names (attributes (graph))) {
        px <- attr (graph, "px")
        while (px$is_alive ()) {
            px$wait ()
        }
    }

    v <- dodgr_vertices (graph)

    if (!is.null (verts)) {
        if (!(length (verts) == 1 || is.vector (verts))) {
            stop ("verts must be a single value or a vector of vertex IDs")
        }
        if (!is.character (verts)) {
            verts <- paste0 (verts)
        }
        verts <- verts [which (verts %in% v$id)]
    }

    hash <- get_hash (graph, contracted = FALSE, force = TRUE)
    hashc <- get_hash (graph, verts = verts, contracted = TRUE, force = TRUE)
    fname_c <- fs::path (
        fs::path_temp (),
        paste0 ("dodgr_graphc_", hashc, ".Rds")
    )

    if (fs::file_exists (fname_c) && !nocache) {
        graph_contracted <- list (graph = readRDS (fname_c))
    } else {
        fname <- fs::path (
            fs::path_temp (),
            paste0 ("dodgr_graph_", hash, ".Rds")
        )
        if (!fs::file_exists (fname)) {
            saveRDS (graph, fname)
        }

        graph_contracted <- dodgr_contract_graph_internal (graph, v, verts)

        gr_cols <- dodgr_graph_cols (graph_contracted$graph)
        hashe <- digest::digest (graph_contracted$graph [[gr_cols$edge_id]])
        attr (graph_contracted$graph, "hash") <- hash
        attr (graph_contracted$graph, "hashc") <- hashc
        attr (graph_contracted$graph, "hashe") <- hashe

        saveRDS (graph_contracted$graph, fname_c)

        fname_e <- fs::path (
            fs::path_temp (),
            paste0 ("dodgr_edge_map_", hashc, ".Rds")
        )
        saveRDS (graph_contracted$edge_map, fname_e)

        fname_j <- fs::path (
            fs::path_temp (),
            paste0 ("dodgr_junctions_", hashc, ".Rds")
        )
        saveRDS (graph_contracted$junctions, fname_j)
    }

    # copy the processx R6 object associated with caching the original graph:
    if (!is.null (px)) {
        attr (graph_contracted$graph, "px") <- px
    }

    return (graph_contracted$graph)
}

# get junction vertices of graphs which have been re-routed for turn angles.
# These all have either "_start" or "_end" appended to vertex names
# v is result of `dodgr_vertices` functions.
get_junction_vertices <- function (v) {
    gsub ("_start|_end", "", v$id [grep ("_start|_end", v$id)])
}

dodgr_contract_graph_internal <- function (graph, v, verts = NULL) {
    classes <- class (graph)
    graph <- tbl_to_df (graph)

    junctions <- get_junction_vertices (v)

    gr_cols <- dodgr_graph_cols (graph)
    graph2 <- convert_graph (graph, gr_cols)
    graph_contracted <- rcpp_contract_graph (graph2, verts)

    graph_contracted <-
        rm_edges_with_heterogenous_data (graph, graph_contracted, gr_cols)

    # graph_contracted$graph has only 5 cols of (edge_id, from, to, d, w). These
    # have to be matched onto original graph.  This is done by using edge_map to
    # get matching indices into both contracted and original graph:
    indx_contr <- match (
        graph_contracted$edge_map$edge_new,
        graph_contracted$graph$edge_id
    )
    indx_orig <- match (
        graph_contracted$edge_map$edge_old,
        graph [, gr_cols$edge_id]
    )
    # Then reduce the latter only to the corresponding first non-repeated values
    # of the former.
    indx_orig <- indx_orig [which (!duplicated (indx_contr))]

    indx_contr <- unique (indx_contr)
    graph_refill <- graph [indx_orig, ]
    graph_refill [, gr_cols$edge_id] <-
        graph_contracted$graph$edge_id [indx_contr]
    graph_refill [, gr_cols$from] <- graph_contracted$graph$from [indx_contr]
    graph_refill [, gr_cols$to] <- graph_contracted$graph$to [indx_contr]
    graph_refill [, gr_cols$d] <- graph_contracted$graph$d [indx_contr]
    graph_refill [, gr_cols$d_weighted] <-
        graph_contracted$graph$d_weighted [indx_contr]
    if (!is.na (gr_cols$time) && !is.na (gr_cols$time_weighted)) {
        graph_refill [, gr_cols$time] <-
            graph_contracted$graph$time [indx_contr]
        graph_refill [, gr_cols$time_weighted] <-
            graph_contracted$graph$timew [indx_contr]
    }

    # Then re-insert spatial coordinates
    if (is_graph_spatial (graph)) {
        spcols <- find_spatial_cols (graph)$fr_col
        indx <- match (graph_contracted$graph$from [indx_contr], graph2$from)
        graph_refill [, spcols [1]] <- graph [indx, spcols [1]] # nolint
        graph_refill [, spcols [2]] <- graph [indx, spcols [2]] # nolint

        spcols <- find_spatial_cols (graph)$to_col
        indx <- match (graph_contracted$graph$to [indx_contr], graph2$to)
        graph_refill [, spcols [1]] <- graph [indx, spcols [1]] # nolint
        graph_refill [, spcols [2]] <- graph [indx, spcols [2]] # nolint
        # This code matches way_id values to those in original graph, but that's
        # kind of arbitrary because with duplicated ways the ID matching can
        # never be systematically controlled
        # if ("way_id" %in% names (graph))
        # {
        #    indx <- match (graph_contracted$graph$edge_id,
        #                   graph_contracted$edge_map$edge_new)
        #    indx <- graph_contracted$edge_map$edge_old [indx]
        #    indx <- match (indx, graph [, gr_cols$edge_id ])
        #    indx2 <- which (!is.na (indx))
        #    indx <- indx [indx2]
        #    graph_refill$way_id [indx2] <- graph$way_id [indx]
        # }
    }

    # and finally replicate the uncontracted edges of graph in graph_contracted
    indx_uncontr <- which (!graph [, gr_cols$edge_id] %in%
        graph_contracted$edge_map$edge_old)
    graph_refill <- rbind (graph_refill, graph [indx_uncontr, ])

    if (any (grepl ("comp", names (graph), ignore.case = TRUE))) {
        ci <- which (grepl ("comp", names (graph_refill), ignore.case = TRUE))
        cnm <- names (graph_refill) [ci]
        graph_refill [[cnm]] <- NULL
        graph_refill <- dodgr_components (graph_refill)
        ci <- which (grepl ("comp", names (graph_refill), ignore.case = TRUE))
        names (graph_refill) [ci] <- cnm
    }

    class (graph_refill) <- c (classes, "dodgr_contracted")

    return (list (
        graph = graph_refill,
        edge_map = graph_contracted$edge_map,
        junctions = junctions
    ))
}

#' Graph contraction must ignore any compound edges along which any additional
#' data columns change (see #194).
#'
#' @return A modified version of `graph_contracted`, removing any formerly
#' contracted edges which should not be, and expanding them back out to their
#' original edges. The "edge_map" component of `graph_contracted` is also
#' modified to remove the corresponding items.
#' @noRd
rm_edges_with_heterogenous_data <- function (graph, graph_contracted, gr_cols) { # nolint

    gr_cols_index <- unlist (gr_cols)
    gr_cols_index <- gr_cols_index [-which (names (gr_cols_index) == "edge_id")]
    rm_these <- c ("geom_num", "highway", "way_id")
    rm_these <- rm_these [which (rm_these %in% names (graph))]
    gr_cols_index <- sort (c (gr_cols_index, match (rm_these, names (graph))))
    graph_extra <- graph [, -gr_cols_index, drop = FALSE]
    data_index <- which (!names (graph_extra) %in%
        c ("edge_id", "edge_new", "edge_", "object_"))
    if (length (data_index) == 0L) {
        return (graph_contracted)
    }

    graph_extra <- graph_extra [
        which (graph_extra$edge_id %in% graph_contracted$edge_map$edge_old),
    ]
    index <- match (graph_extra$edge_id, graph_contracted$edge_map$edge_old)
    graph_extra$edge_new <- graph_contracted$edge_map$edge_new [index]

    graph_extra <- split (graph_extra, f = factor (graph_extra$edge_new))
    graph_extra <- lapply (graph_extra, function (i) {
        unique (i [, data_index, drop = FALSE])
    })
    lens <- vapply (graph_extra, nrow, integer (1L))
    index_heterog <- names (lens) [which (lens > 1L)]
    if (length (index_heterog) == 0L) {
        return (graph_contracted)
    }

    edge_map_heterog <- graph_contracted$edge_map [
        graph_contracted$edge_map$edge_new %in% index_heterog,
    ]
    index_out <- match (index_heterog, graph_contracted$graph$edge_id)
    index_in <- match (edge_map_heterog$edge_old, graph [, gr_cols$edge_id])
    graph2 <- convert_graph (graph, gr_cols)
    names (graph2) <- names (graph_contracted$graph)
    graph_contracted$graph <- rbind (
        graph_contracted$graph [-index_out, , drop = FALSE],
        graph2 [index_in, , drop = FALSE]
    )
    edge_map_index <-
        match (edge_map_heterog$edge_old, graph_contracted$edge_map$edge_old)
    graph_contracted$edge_map <- graph_contracted$edge_map [-edge_map_index, ]

    return (graph_contracted)
}

#' Re-expand a contracted graph.
#'
#' Revert a contracted graph created with \link{dodgr_contract_graph} back to
#' the full, uncontracted version. This function is mostly used for the side
#' effect of mapping any new columns inserted on to the contracted graph back
#' on to the original graph, as demonstrated in the example.
#'
#' @param graph A contracted graph created from \link{dodgr_contract_graph}.
#'
#' @return A single `data.frame` representing the equivalent original,
#' uncontracted graph.
#' @family modification
#' @export
#' @examples
#' graph0 <- weight_streetnet (hampi)
#' nrow (graph0) # 6,813
#' graph1 <- dodgr_contract_graph (graph0)
#' nrow (graph1) # 760
#' graph2 <- dodgr_uncontract_graph (graph1)
#' nrow (graph2) # 6,813
#'
#' # Insert new data on to the contracted graph and uncontract it:
#' graph1$new_col <- runif (nrow (graph1))
#' graph3 <- dodgr_uncontract_graph (graph1)
#' # graph3 is then the uncontracted graph which includes "new_col" as well
#' dim (graph0)
#' dim (graph3)
dodgr_uncontract_graph <- function (graph) {

    px <- NULL
    if ("px" %in% names (attributes (graph))) {
        px <- attr (graph, "px") # processx R6 object
        while (px$is_alive ()) {
            px$wait ()
        }
    }

    edge_map <- get_edge_map (graph)

    gr_cols <- dodgr_graph_cols (graph)
    hashe_ref <- attr (graph, "hashe")
    hashe <- digest::digest (graph [[gr_cols$edge_id]])

    hash <- attr (graph, "hash")
    fname <- fs::path (fs::path_temp (), paste0 ("dodgr_graph_", hash, ".Rds"))
    if (!fs::file_exists (fname)) {
        stop (paste0 (
            "Graph must have been contracted in ", # nocov
            "current R session; and have retained ", # nocov
            "the same row structure"
        ))
    } # nocov

    # used below if rows have been removed
    graph_edges <- graph [[gr_cols$edge_id]]

    graph_full <- readRDS (fname)
    attr (graph_full, "px") <- px

    graph <- uncontract_graph (graph, edge_map, graph_full)

    # Finally, remove any edges which might have been removed from contracted
    # graph:
    if (!identical (hashe, hashe_ref)) {
        index <- which (edge_map$edge_new %in% graph_edges)
        index_in <-
            which (graph [[gr_cols$edge_id]] %in% edge_map$edge_old [index])
        graph <- graph [index_in, ]
        attr (graph, "hash") <-
            get_hash (graph, contracted = FALSE, force = TRUE)
    }

    return (graph)
}

# map contracted graph with flows (or whatever else) back onto full graph
uncontract_graph <- function (graph, edge_map, graph_full) {

    gr_cols <- dodgr_graph_cols (graph_full)
    indx_to_full <- match (edge_map$edge_old, graph_full [[gr_cols$edge_id]])
    indx_to_contr <- match (edge_map$edge_new, graph [[gr_cols$edge_id]])
    # edge_map only has the contracted edges; flows from the original
    # non-contracted edges also need to be inserted
    edges <- graph [[gr_cols$edge_id]] [which (!graph [[gr_cols$edge_id]] %in%
        edge_map$edge_new)]
    indx_to_full <- c (
        indx_to_full,
        match (edges, graph_full [[gr_cols$edge_id]])
    )
    indx_to_contr <- c (
        indx_to_contr,
        match (edges, graph [[gr_cols$edge_id]])
    )

    index <- which (!names (graph) %in% names (graph_full))
    if (length (index) > 0) {
        nms <- names (graph) [index]
        graph_full [nms] <- NA
        for (n in nms) {
            graph_full [[n]] [indx_to_full] <- graph [[n]] [indx_to_contr]
        }
    }

    return (graph_full)
}
