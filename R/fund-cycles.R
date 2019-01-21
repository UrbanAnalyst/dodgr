#' dodgr_fundamental_cycles
#'
#' Calculate fundamental cycles in a graph
#'
#' @param graph `data.frame` or equivalent object representing the network
#' graph (see Details)
#' @param vertices `data.frame` returned from \link{dodgr_vertices}`(graph)`. Will
#' be calculated if not provided, but it's quicker to pass this if it has
#' already been calculated.
#' @return List of cycle paths, in terms of vertex IDs in `graph`.
#'
#' @examples 
#' net <- weight_streetnet (hampi)
#' graph <- dodgr_contract_graph (net)$graph
#' verts <- dodgr_vertices (graph)
#' cyc <- dodgr_fundamental_cycles (graph, verts)
#' @export 
dodgr_fundamental_cycles <- function (graph, vertices = NULL)
{
    if (missing (graph))
        stop ("graph must be provided")
    if (!inherits (graph, "data.frame"))
        stop ("graph must be a data.frame object")

    if (is.null (vertices))
        vertices <- dodgr_vertices (graph)
    graph <- convert_graph (graph, dodgr_graph_cols (graph))
    rcpp_fundamental_cycles (graph, vertices)
}
