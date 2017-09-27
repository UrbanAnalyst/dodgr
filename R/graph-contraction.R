#' dodgr_contract_graph
#'
#' Removes redundant (straight-line) vertices from graph, leaving only junction
#' vertices.
#'
#' @param graph A flat table of graph edges. Must contain columns labelled
#' \code{from} and \code{to}, or \code{start} and \code{stop}. May also contain
#' similarly labelled columns of spatial coordinates (for example
#' \code{from_x}) or \code{stop_lon}).
#' @param verts Option \code{data.frame} of vertices obtained from
#' \code{dodgr_vertices} (submitting this will simply speed up conversion to
#' compact graph).
#' @param quiet If \code{FALSE}, display progress on screen
#'
#' @return A list with the original graph (\code{$origina}), the contracted
#' graph (\code{$contracted}), and a map (\code{map}) between the IDs of edges
#' in the contracted graph and those in the original, full graph.
#' @export
#' @examples
#' graph <- weight_streetnet (hampi)
#' nrow (graph) # 5,742
#' graph <- dodgr_contract_graph (graph)
#' nrow (graph$original) # size of contracted graph = 2,878
#' nrow (graph$map) # 3,692 = number of old edges replaced by new ones
dodgr_contract_graph <- function (graph, verts = NULL, quiet = TRUE)
{
    graph_converted <- dodgr_convert_graph (graph)$graph
    graph_contracted <- rcpp_contract_graph (graph_converted, quiet = quiet)

    # re-insert spatial coordinates:
    if (is_graph_spatial (graph))
    {
        spcols <- find_spatial_cols (graph)
        fr_col <- find_fr_id_col (graph)
        indx <- match (graph_contracted$contracted$from, graph [, fr_col])
        fr_xy <- graph [indx, spcols$fr_col]
        to_col <- find_to_id_col (graph)
        indx <- match (graph_contracted$contracted$to, graph [, to_col])
        to_xy <- graph [indx, spcols$to_col]
        graph_contracted$contracted <- cbind (graph_contracted$contracted,
                                         fr_xy, to_xy)
    }

    graph_contracted$original <- graph
    if ("component" %in% names (graph))
        graph_contracted$contracted <-
            dodgr_components (graph_contracted$contracted)

    return (graph_contracted)
}

#' dodgr_reinsert_verts
#'
#' Re-inserts vertices in a contracted graph. Vertices used for routing may often
#' be redundant - that is, be intermediate vertices between just two other
#' points - and will be removed in \code{\link{dodgr_contract_graph}}. This
#' function re-inserts vertices for routing prior to submitting graph to
#' \code{\link{dodgr_dists}}.
#'
#' @param gc The result of call to \code{\link{dodgr_contract_graph}}, a list of
#' the three items of (i) the original graph; (ii) the contracted graph; and
#' (iii) an edge map related contracted to original edges.
#' @param verts List of vertices to be used for routing, or equivalent list of
#' geographical coordinates.
#'
#' @return A modified version of \code{gc}, with the contracted graph
#' (\code{$contracted}) including all vertices contained in \code{verts}.
#' @export
#' @examples
#' graph <- weight_streetnet (hampi)
#' gc <- dodgr_contract_graph (graph)
#' nrow (gc$contracted) # size of contracted graph = 3,033
#' verts <- dodgr_vertices (graph)
#' verts <- verts [sample (nrow (verts), 100), ]
#' gc <- dodgr_reinsert_verts (gc, verts = verts)
#' nrow (gc$contracted) # will be > 2,878
dodgr_reinsert_verts <- function (gc, verts = NULL)
{
    ids <- grepl ("id", names (verts), ignore.case = TRUE)
    if (any (ids))
        verts <- verts [[which (ids)]]
    else if (is.matrix (verts) && ncol (verts) > 1)
        verts <- match_pts_to_graph (dodgr_vertices (gc$original), xy = verts)
    else if (!is.vector (verts))
        stop ("format of verts not recognised.")

    verts_contr <- dodgr_vertices (gc$contracted)
    verts <- verts [!verts %in% verts_contr$id]

    if (length (verts) > 0)
    {
        orig <- dodgr_convert_graph (gc$orig)$graph
        contr <- dodgr_convert_graph (gc$contracted)$graph

        edges <- rcpp_insert_vertices (orig, contr, gc$map, verts)
        indx_insert <- match (edges$insert, orig$edge_id)
        contr <- rbind (contr, orig [indx_insert, ])
        indx_erase <- seq (nrow (contr)) [!seq (nrow (contr)) %in% 
                                          match (edges$erase, contr$edge_id)]
        gc$contracted <- contr [indx_erase, ]
    }

    # re-insert spatial columns in contracted
    if (is_graph_spatial (gc$original))
    {
        spcols <- find_spatial_cols (gc$original)
        fr_col <- find_fr_id_col (gc$original)
        indx <- match (gc$contracted$from, gc$original [, fr_col])
        fr_xy <- gc$original [indx, spcols$fr_col]
        to_col <- find_to_id_col (gc$original)
        indx <- match (gc$contracted$to, gc$original [, to_col])
        to_xy <- gc$original [indx, spcols$to_col]
        gc$contracted <- cbind (gc$contracted, fr_xy, to_xy)
    }

    if ("component" %in% names (gc$contracted))
    {
        gc$contracted$component <- NULL
        gc$contracted <- dodgr_components (gc$contracted)
    }

    return (gc)
}
