#' dodgr_contract_graph
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
#'
#' @return A list of two items: `graph` containing contracted version of
#' the original `graph`, converted to a standardised format, and
#' `edge_map`, a two-column matrix mapping all newly contracted edges onto
#' corresponding edges in original (uncontracted) graph.
#' @export
#' @examples
#' graph <- weight_streetnet (hampi)
#' nrow (graph) # 5,845
#' graph <- dodgr_contract_graph (graph)
#' nrow (graph) # 686
dodgr_contract_graph <- function (graph, verts = NULL)
{
    if (nrow (graph) == 0)
        stop ("graph is empty") # nocov

    # px is the R6 processx object for initial caching
    px <- NULL
    if ("px" %in% names (attributes (graph)))
    {
        px <- attr (graph, "px")
        while (px$is_alive ())
            px$wait ()
    }

    v <- dodgr_vertices (graph)

    if (!is.null (verts))
    {
        if (!(length (verts) == 1 | is.vector (verts)))
            stop ("verts must be a single value or a vector of vertex IDs")
        if (!is.character (verts))
            verts <- paste0 (verts)
        verts <- verts [which (verts %in% v$id)]
    }

    hash <- get_hash (graph, hash = TRUE)
    hashc <- get_hash (graph, verts = verts, hash = FALSE)
    fname_c <- file.path (tempdir (), paste0 ("dodgr_graphc_", hashc, ".Rds"))

    if (file.exists (fname_c))
    {
        graph_contracted <- list (graph = readRDS (fname_c))
    } else
    {
        fname <- file.path (tempdir (), paste0 ("dodgr_graph_", hash, ".Rds"))
        if (!file.exists (fname))
            saveRDS (graph, fname)

        graph_contracted <- dodgr_contract_graph_internal (graph, v, verts)

        gr_cols <- dodgr_graph_cols (graph_contracted$graph)
        hashe <- digest::digest (graph_contracted$graph [[gr_cols$edge_id]])
        attr (graph_contracted$graph, "hash") <- hash
        attr (graph_contracted$graph, "hashc") <- hashc
        attr (graph_contracted$graph, "hashe") <- hashe

        saveRDS (graph_contracted$graph, fname_c)

        fname_e <- file.path (tempdir (), paste0 ("dodgr_edge_map_", hashc, ".Rds"))
        saveRDS (graph_contracted$edge_map, fname_e)

        fname_j <- file.path (tempdir (), paste0 ("dodgr_junctions_", hashc, ".Rds"))
        saveRDS (graph_contracted$junctions, fname_j)
    }

    # copy the processx R6 object associated with caching the original graph:
    if (!is.null (px))
        attr (graph_contracted$graph, "px") <- px

    return (graph_contracted$graph)
}

# get junction vertices of graphs which have been re-routed for turn angles.
# These all have either "_start" or "_end" appended to vertex names
# v is result of `dodgr_vertices` functions.
get_junction_vertices <- function (v)
{
    gsub ("_start|_end", "", v$id [grep ("_start|_end", v$id)])
}

dodgr_contract_graph_internal <- function (graph, v, verts = NULL)
{
    classes <- class (graph)
    graph <- tbl_to_df (graph)

    junctions <- get_junction_vertices (v)

    gr_cols <- dodgr_graph_cols (graph)
    graph2 <- convert_graph (graph, gr_cols)
    graph_contracted <- rcpp_contract_graph (graph2, verts)

    # graph_contracted$graph has only 5 cols of (edge_id, from, to, d, w). These
    # have to be matched onto original graph.  This is done by using edge_map to
    # get matching indices into both contracted and original graph:
    indx_contr <- match (graph_contracted$edge_map$edge_new,
                         graph_contracted$graph$edge_id)
    indx_orig <- match (graph_contracted$edge_map$edge_old,
                        graph [, gr_cols$edge_id])
    # Then reduce the latter only to the corresponding first non-repeated values
    # of the former.
    indx_orig <- indx_orig [which (!duplicated (indx_contr))]

    indx_contr <- unique (indx_contr)
    graph_refill <- graph [indx_orig, ]
    graph_refill [, gr_cols$edge_id] <- graph_contracted$graph$edge_id [indx_contr]
    graph_refill [, gr_cols$from] <- graph_contracted$graph$from [indx_contr]
    graph_refill [, gr_cols$to] <- graph_contracted$graph$to [indx_contr]
    graph_refill [, gr_cols$d] <- graph_contracted$graph$d [indx_contr]
    graph_refill [, gr_cols$w] <- graph_contracted$graph$w [indx_contr]
    if (!is.na (gr_cols$time) & !is.na (gr_cols$time_weighted))
    {
        graph_refill [, gr_cols$time] <- graph_contracted$graph$time [indx_contr]
        graph_refill [, gr_cols$time_weighted] <-
            graph_contracted$graph$timew [indx_contr]
    }

    # Then re-insert spatial coordinates
    if (is_graph_spatial (graph))
    {
        spcols <- find_spatial_cols (graph)$fr_col
        indx <- match (graph_contracted$graph$from [indx_contr], graph2$from)
        graph_refill [, spcols [1] ] <- graph [indx, spcols [1] ]
        graph_refill [, spcols [2] ] <- graph [indx, spcols [2] ]

        spcols <- find_spatial_cols (graph)$to_col
        indx <- match (graph_contracted$graph$to [indx_contr], graph2$to)
        graph_refill [, spcols [1] ] <- graph [indx, spcols [1] ]
        graph_refill [, spcols [2] ] <- graph [indx, spcols [2] ]
        # This code matches way_id values to those in original graph, but that's
        # kind of arbitrary because with duplicated ways the ID matching can
        # never be systematically controlled
        #if ("way_id" %in% names (graph))
        #{
        #    indx <- match (graph_contracted$graph$edge_id,
        #                   graph_contracted$edge_map$edge_new)
        #    indx <- graph_contracted$edge_map$edge_old [indx]
        #    indx <- match (indx, graph [, gr_cols$edge_id ])
        #    indx2 <- which (!is.na (indx))
        #    indx <- indx [indx2]
        #    graph_refill$way_id [indx2] <- graph$way_id [indx]
        #}
    }

    # and finally replicate the uncontracted edges of graph in graph_contracted 
    indx_uncontr <- which (!graph [, gr_cols$edge_id] %in%
                           graph_contracted$edge_map$edge_old)
    graph_refill <- rbind (graph_refill, graph [indx_uncontr, ])

    if (any (grepl ("comp", names (graph), ignore.case = TRUE)))
    {
        ci <- which (grepl ("comp", names (graph_refill), ignore.case = TRUE))
        cnm <- names (graph_refill) [ci]
        graph_refill [[cnm]] <- NULL
        graph_refill <- dodgr_components (graph_refill)
        ci <- which (grepl ("comp", names (graph_refill), ignore.case = TRUE))
        names (graph_refill) [ci] <- cnm
    }

    class (graph_refill) <- c (classes, "dodgr_contracted")

    return (list (graph = graph_refill,
                  edge_map = graph_contracted$edge_map,
                  junctions = junctions))
}

#' dodgr_uncontract_graph
#'
#' Revert a contracted graph created with \link{dodgr_contract_graph} back to
#' the full, uncontracted version. This function is mostly used for the side
#' effect of mapping any new columnns inserted on to the contracted graph back
#' on to the original graph, as demonstrated in the example.
#'
#' @param graph A list of two items returned from \link{dodgr_contract_graph},
#' the first ("graph") containing the contracted graph, and the second
#' ("edge_map") mapping edges in the contracted graph back to those in the
#' original graph.
#'
#' @return A single `data.frame` representing the original, uncontracted graph.
#' @export
#' @examples
#' graph0 <- weight_streetnet (hampi)
#' nrow (graph0) # 5,845
#' graph1 <- dodgr_contract_graph (graph0)
#' nrow (graph1) # 686
#' graph2 <- dodgr_uncontract_graph (graph1)
#' nrow (graph2) # 5,845
#' 
#' # Insert new data on to the contracted graph and uncontract it:
#' graph1$new_col <- runif (nrow (graph1))
#' graph3 <- dodgr_uncontract_graph (graph1)
#' # graph3 is then the uncontracted graph which includes "new_col" as well
#' dim (graph0); dim (graph3)
dodgr_uncontract_graph <- function (graph)
{
    px <- attr (graph, "px") # processx R6 object
    while (px$is_alive ())
        px$wait ()

    edge_map <- get_edge_map (graph)

    gr_cols <- dodgr_graph_cols (graph)
    hashe_ref <- attr (graph, "hashe")
    hashe <- digest::digest (graph [[gr_cols$edge_id]])
    if (!identical (hashe, hashe_ref))
        stop ("Unable to uncontract this graph because the rows ",
              "have been changed")

    hash <- attr (graph, "hash")
    fname <- file.path (tempdir (), paste0 ("dodgr_graph_", hash, ".Rds"))
    if (!file.exists (fname))
        stop (paste0 ("Graph must have been contracted in current R session; ", # nocov
                      "and have retained the same row structure"))              # nocov

    graph_full <- readRDS (fname)
    attr (graph_full, "px") <- px

    graph <- uncontract_graph (graph, edge_map, graph_full)

    tp <- attr (graph, "turn_penalty")
    tp <- ifelse (is.null (tp), 0, tp)
    if (is (graph, "dodgr_streetnet_sc") & tp > 0)
    {
        # extra code to uncontract the compound turn-angle junctions, including
        # merging extra rows such as flow from compound junctions back into
        # "normal" (non-compound) edges
        hash <- get_hash (graph, hash = TRUE)
        fname <- file.path (tempdir (), paste0 ("dodgr_edge_contractions_",
                                                hash, ".Rds"))
        if (!file.exists (fname))
            stop (paste0 ("Graph must have been contracted in current R ",      # nocov
                          "session; and have retained the same row structure")) # nocov
        ec <- readRDS (fname)

        index_junction <- match (ec$edge, graph [[gr_cols$edge_id]])
        index_edge_in <- match (ec$e_in, graph [[gr_cols$edge_id]])
        index_edge_out <- match (ec$e_out, graph [[gr_cols$edge_id]])
        new_cols <- names (graph) [which (!names (graph) %in% 
                                          names (graph_full))]
        for (n in new_cols)
        {
            graph [[n]] [index_edge_in] <- graph [[n]] [index_edge_in] +
                graph [[n]] [index_junction]
            graph [[n]] [index_edge_out] <- graph [[n]] [index_edge_out] +
                graph [[n]] [index_junction]
        }
        # next line removes all the compound turn angle edges:
        graph <- graph [which (!graph [[gr_cols$edge_id]] %in% ec$edge), ]
        graph$.vx0 <- gsub ("_start$", "", graph$.vx0)
        graph$.vx1 <- gsub ("_end$", "", graph$.vx1)
    }

    return (graph)
}

# map contracted graph with flows (or whatever else) back onto full graph
uncontract_graph <- function (graph, edge_map, graph_full)
{
    gr_cols <- dodgr_graph_cols (graph_full)
    indx_to_full <- match (edge_map$edge_old, graph_full [[gr_cols$edge_id]])
    indx_to_contr <- match (edge_map$edge_new, graph [[gr_cols$edge_id]])
    # edge_map only has the contracted edges; flows from the original
    # non-contracted edges also need to be inserted
    edges <- graph [[gr_cols$edge_id]] [which (!graph [[gr_cols$edge_id]] %in%
                                               edge_map$edge_new)]
    indx_to_full <- c (indx_to_full, match (edges, graph_full [[gr_cols$edge_id]]))
    indx_to_contr <- c (indx_to_contr, match (edges, graph [[gr_cols$edge_id]]))

    index <- which (!names (graph) %in% names (graph_full))
    if (length (index) > 0)
    {
        nms <- names (graph) [index]
        graph_full [nms] <- NA
        for (n in nms)
        {
            graph_full [[n]] [indx_to_full] <- graph [[n]] [indx_to_contr]
        }
    }

    return (graph_full)
}

