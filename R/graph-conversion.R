#' dodgr_to_igraph
#'
#' Convert a `dodgr` graph to an \pkg{igraph}.
#'
#' @param graph A `dodgr` graph
#'
#' @return The `igraph` equivalent of the input
#' @export
#' @examples
#' graph <- weight_streetnet (hampi)
#' graphi <- dodgr_to_igraph (graph)
dodgr_to_igraph <- function (graph)
{
    v <- dodgr_vertices (graph)
    graph <- graph [, c (3, 6, 1:2, 4:5, 7:13)]
    igraph::graph_from_data_frame (graph, directed = TRUE, vertices = v)
}


#' dodgr_to_tidygraph
#'
#' Convert a `dodgr` graph to an \pkg{tidygraph}.
#'
#' @param graph A `dodgr` graph
#'
#' @return The `tidygraph` equivalent of the input
#' @export
#' @examples
#' graph <- weight_streetnet (hampi)
#' grapht <- dodgr_to_tidygraph (graph)
dodgr_to_tidygraph <- function (graph)
{
    if (!requireNamespace ("tidygraph"))
        stop ("dodgr_to_tidygraph requires the tidygraph package to be installed.")
    dodgr_to_igraph (graph) %>%
        tidygraph::as_tbl_graph ()
}
