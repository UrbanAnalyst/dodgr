#' dodgr_times
#'
#' Calculate matrix of pair-wise travel times between points.
#'
#' @inherit dodgr_dists
#' @param graph A `dodgr` network returned from the \link{weight_streetnet}
#' function using a network obtained with the \pkg{osmdata} `osmdata_sc`
#' function, possibly contracted with \link{dodgr_contract_graph}.
#'
#' @export 
dodgr_times <- function (graph, from = NULL, to = NULL, heap = 'BHeap')
{
    gr_cols <- dodgr_graph_cols (graph)
    if (is.na (gr_cols$time))
        stop ("graph has no time column")

    graph <- tbl_to_df (graph)

    graph [[gr_cols$d]] <- graph [[gr_cols$time]]

    dodgr_dists (graph, from, to, heap = heap)
}
