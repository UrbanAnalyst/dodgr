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
#' @return A contracted version of the original \code{graph}, converted to a
#' standardised format.
#' @export
#' @examples
#' graph <- weight_streetnet (hampi)
#' nrow (graph) # 5,742
#' graph <- dodgr_contract_graph (graph)
#' nrow (graph) # 2,878
dodgr_contract_graph <- function (graph, verts = NULL)
{
    if (!is.null (verts))
    {
        if (!(length (verts) == 1 | is.vector (verts)))
            stop ("verts must be a single value or a vector of vertex IDs")
        if (!is.character (verts))
            verts <- paste0 (verts)
    }

    graph_converted <- dodgr_convert_graph (graph)$graph
    verts <- verts [which (verts %in% dodgr_vertices (graph_converted)$id)]
    graph_contracted <- rcpp_contract_graph (graph_converted, verts)

    # re-insert spatial coordinates:
    if (is_graph_spatial (graph))
    {
        spcols <- find_spatial_cols (graph)
        fr_col <- find_fr_id_col (graph)
        indx <- match (graph_contracted$from, graph [, fr_col])
        fr_xy <- graph [indx, spcols$fr_col]
        to_col <- find_to_id_col (graph)
        indx <- match (graph_contracted$to, graph [, to_col])
        to_xy <- graph [indx, spcols$to_col]
        graph_contracted <- cbind (graph_contracted, fr_xy, to_xy)
    }

    if ("component" %in% names (graph))
        graph_contracted <- dodgr_components (graph_contracted)

    return (graph_contracted)
}
