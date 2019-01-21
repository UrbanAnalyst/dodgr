#' dodgr_fundamental_cycles
#'
#' Calculate fundamental cycles in a graph
#'
#' @param graph `data.frame` or equivalent object representing the contracted
#' network graph (see Details).
#' @param vertices `data.frame` returned from \link{dodgr_vertices}`(graph)`. Will
#' be calculated if not provided, but it's quicker to pass this if it has
#' already been calculated.
#' @return List of cycle paths, in terms of vertex IDs in `graph` and, for
#' spatial graphs, the corresponding coordinates.
#'
#' @note Calculation of fundamental cycles is VERY computationally demanding,
#' and this function should only be executed on CONTRACTED graphs (that is,
#' graphs returned from \link{dodgr_contract_graph}). Results for full graphs
#' can be obtained with the function \link{dodgr_full_cycles}.
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
    if (!"flow" %in% names (graph)) # makes no difference
        graph$flow <- 1
    graph <- merge_directed_flows (graph) # uses fast C++ routines
    graph$flow <- NULL
    graphc <- convert_graph (graph, dodgr_graph_cols (graph))
    res <- rcpp_fundamental_cycles (graphc, vertices)

    if (is_graph_spatial (graph))
    {
        res <- lapply (res, function (i) {
                           index <- match (i, vertices$id)
                           data.frame (id = i,
                                       x = vertices$x [index],
                                       y = vertices$y [index],
                                       stringsAsFactors = FALSE)
                                })
    }
    return (res)
}

#' dodgr_full_cycles
#'
#' Calculate fundamental cycles on a FULL (that is, non-contracted) graph.
#' @inheritParams dodgr_fundamental_cycles
#' @note This function converts the `graph` to its contracted form, calculates
#' the fundamental cycles on that version, and then expands these cycles back
#' onto the original graph. This is far more computationally efficient than
#' calculating fundamental cycles on a full (non-contracted) graph.
#' @export
dodgr_full_cycles <- function (graph)
{
    graphc <- dodgr_contract_graph (graph)
    v <- dodgr_vertices (graphc$graph)

    x <- dodgr_fundamental_cycles (graphc$graph, vertices = v)

    from_to <- paste0 (graphc$graph$from_id, "-", graphc$graph$to_id)
    ids <- lapply (x, function (i) {
           idpairs <- paste0 (i$id [-length (i$id)], "-",
                              to = i$id [-1])
           edges_c <- graphc$graph$edge_id [match (idpairs, from_to)]
           edges_new <- lapply (as.list (edges_c), function (j) {
                    if (j %in% graphc$edge_map$edge_new)
                    {
                        index <- which (graphc$edge_map$edge_new %in% j)
                        j <- as.numeric (graphc$edge_map$edge_old [index])
                        if ((graph$from_id [j [1] ] == graph$to_id [utils::tail (j, 1)]) ||
                            (graph$from_id [j [1] ] == graph$to_id [j [2] ]))
                            j <- rev (j)
                    }
                    return (as.numeric (j))
                    }) # end edges_new lapply
           unlist (edges_new)
            }) # end ids lapply
    # ids at that point is a sequences of indices into graph. This is then
    # converted to a sequence of vertex IDs, through just adding the last vertex
    # of the sequence on to close the polygon
    gr_cols <- dodgr_graph_cols (graph) # (from, to) = [, 2:3]
    res <- lapply (ids, function (i)
                   graph [c (i, i [1]), gr_cols [2]] )

    if (is_graph_spatial (graph))
    {
        vertices <- dodgr_vertices (graph)
        res <- lapply (res, function (i) {
                           index <- match (i, vertices$id)
                           data.frame (id = i,
                                       x = vertices$x [index],
                                       y = vertices$y [index],
                                       stringsAsFactors = FALSE)
                                })
    }
    return (res)
}
