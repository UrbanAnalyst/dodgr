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
#' These must match the `from_id` and `to_id` columns of `graph`.
#'
#' @return A list of two items: `graph` containing contracted version of
#' the original `graph`, converted to a standardised format, and
#' `edge_map`, a two-column matrix mapping all newly contracted edges onto
#' corresponding edges in original (uncontracted) graph.
#' @export
#' @examples
#' graph <- weight_streetnet (hampi)
#' nrow (graph) # 5,729
#' graph <- dodgr_contract_graph (graph)
#' nrow (graph$graph) # 764
dodgr_contract_graph <- function (graph, verts = NULL)
{
    classes <- class (graph)
    graph <- tbl_to_df (graph)
    if (nrow (graph) == 0)
        stop ("graph is empty")

    v <- dodgr_vertices (graph)
    junctions <- get_junction_vertices (v)

    if (!is.null (verts))
    {
        if (!(length (verts) == 1 | is.vector (verts)))
            stop ("verts must be a single value or a vector of vertex IDs")
        if (!is.character (verts))
            verts <- paste0 (verts)
        verts <- verts [which (verts %in% v$id)]
    }

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

    dig <- digest::digest (graph_contracted$edge_map)
    fname <- file.path (tempdir (), paste0 ("graph_", dig, ".Rds"))
    saveRDS (graph, file = fname)

    class (graph_refill) <- c (classes, "dodgr_contracted")

    return (list (graph = graph_refill,
                  edge_map = graph_contracted$edge_map,
                  junctions = junctions))
}

# get junction vertices of graphs which have been re-routed for turn angles.
# These all have either "_start" or "_end" appended to vertex names
# v is result of `dodgr_vertices` functions.
get_junction_vertices <- function (v)
{
    gsub ("_start|_end", "", v$id [grep ("_start|_end", v$id)])
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
#' nrow (graph0) # 5,729
#' graph1 <- dodgr_contract_graph (graph0)
#' nrow (graph1$graph) # 764
#' graph2 <- dodgr_uncontract_graph (graph1)
#' nrow (graph2) # 5,729
#' identical (graph0, graph2) # TRUE
#' 
#' # Insert new data on to the contracted graph and uncontract it:
#' graph1$graph$new_col <- runif (nrow (graph1$graph))
#' graph3 <- dodgr_uncontract_graph (graph1)
#' # graph3 is then the uncontracted graph which includes "new_col" as well
dodgr_uncontract_graph <- function (graph)
{
    dig <- digest::digest (graph$edge_map)
    fname <- file.path (tempdir (), paste0 ("graph_", dig, ".Rds"))
    if (!file.exists (fname))
        stop ("Graph must have been contracted in current R session")
    graph_full <- readRDS (fname)

    uncontract_graph (graph$graph, graph$edge_map, graph_full)
}
