#' dodgr_contract_graph
#'
#' Removes redundant (straight-line) vertices from graph, leaving only junction
#' vertices.
#'
#' @param graph A flat table of graph edges. Must contain columns labelled
#' \code{from} and \code{to}, or \code{start} and \code{stop}. May also contain
#' similarly labelled columns of spatial coordinates (for example
#' \code{from_x}) or \code{stop_lon}).
#' @param verts Optional list of vertices to be retained as routing points.
#' These must match the \code{from_id} and \code{to_id} columns of \code{graph}.
#'
#' @return A list of two items: \code{graph} containing contracted version of
#' the original \code{graph}, converted to a standardised format, and
#' \code{edge_map}, a two-column matrix mapping all newly contracted edges onto
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
    # have to be matched onto original graph, and the remaining columns padded
    # out.
    indx <- match (graph_contracted$graph$from, graph2$from)
    graph_refill <- graph [indx, ]
    graph_refill [, gr_cols [1] ] <- graph_contracted$graph$edge_id
    graph_refill [, gr_cols [2] ] <- graph_contracted$graph$from
    graph_refill [, gr_cols [3] ] <- graph_contracted$graph$to
    graph_refill [, gr_cols [4] ] <- graph_contracted$graph$d
    graph_refill [, gr_cols [5] ] <- graph_contracted$graph$w

    # Then re-insert spatial coordinates (only to verts; from already there)
    if (is_graph_spatial (graph))
    {
        indx <- match (graph_contracted$graph$to, graph2$to)
        spcols <- find_spatial_cols (graph)$to_col
        graph_refill [, spcols [1] ] <- graph [indx, spcols [1] ]
        graph_refill [, spcols [2] ] <- graph [indx, spcols [2] ]
        # Then just make sure way_id values match those in original graph
        if ("way_id" %in% names (graph))
        {
            indx <- match (graph_contracted$graph$edge_id,
                           graph_contracted$edge_map$edge_new)
            indx <- graph_contracted$edge_map$edge_old [indx]
            indx <- match (indx, graph [, gr_cols [1] ])
            indx2 <- which (!is.na (indx))
            indx <- indx [indx2]
            graph_refill$way_id [indx2] <- graph$way_id [indx]
        }
    }

    if (any (grepl ("comp", names (graph), ignore.case = TRUE)))
    {
        ci <- which (grepl ("comp", names (graph), ignore.case = TRUE))
        cnm <- names (graph) [ci]
        graph_refill [[cnm]] <- NULL
        graph_refill <- dodgr_components (graph_refill)
        names (graph_refill) [ci] <- cnm
    }

    return (list (graph = graph_refill, edge_map = graph_contracted$edge_map))
}
