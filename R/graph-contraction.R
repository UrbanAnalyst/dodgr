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
#' nrow (graph) # 5,742
#' graph <- dodgr_contract_graph (graph)
#' nrow (graph$graph) # 2,878
dodgr_contract_graph <- function (graph, verts = NULL)
{
    if (!is.null (verts))
    {
        if (!(length (verts) == 1 | is.vector (verts)))
            stop ("verts must be a single value or a vector of vertex IDs")
        if (!is.character (verts))
            verts <- paste0 (verts)
        verts <- verts [which (verts %in% dodgr_vertices (graph)$id)]
    }

    gr_cols <- dodgr_graph_cols (graph)
    graph2 <- convert_graph (graph, gr_cols)
    graph_contracted <- rcpp_contract_graph (graph2, verts)

    # graph_contracted$graph has only 5 cols of (edge_id, from, to, d, w). These
    # have to be matched onto original graph.  This is done by using edge_map to
    # get matching indices into both contracted and original graph:
    indx_contr <- match (graph_contracted$edge_map$edge_new,
                         graph_contracted$graph$edge_id)
    indx_orig <- match (graph_contracted$edge_map$edge_old, graph$edge_id)
    # Then reduce the latter only to the corresponding first non-repeated values of
    # the former
    indx_orig <- indx_orig [which (!duplicated (indx_contr))]
    indx_contr <- unique (indx_contr)
    graph_refill <- graph [indx_orig, ]
    graph_refill [, gr_cols [1] ] <- graph_contracted$graph$edge_id [indx_contr]
    graph_refill [, gr_cols [2] ] <- graph_contracted$graph$from [indx_contr]
    graph_refill [, gr_cols [3] ] <- graph_contracted$graph$to [indx_contr]
    graph_refill [, gr_cols [4] ] <- graph_contracted$graph$d [indx_contr]
    graph_refill [, gr_cols [5] ] <- graph_contracted$graph$w [indx_contr]

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
        #    indx <- match (indx, graph [, gr_cols [1] ])
        #    indx2 <- which (!is.na (indx))
        #    indx <- indx [indx2]
        #    graph_refill$way_id [indx2] <- graph$way_id [indx]
        #}
    }

    # and finally replicate the uncontracted edges of graph in graph_contracted 
    indx_uncontr <- which (!graph$edge_id %in%
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

    return (list (graph = graph_refill, edge_map = graph_contracted$edge_map))
}
