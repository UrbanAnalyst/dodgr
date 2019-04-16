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
    graph <- tbl_to_df (graph)

    #graph$d <- graph$d_weighted <- graph$time
    #graph$d <- graph$time
    graph$time <- NULL

    dodgr_dists (graph, from, to, heap = heap)
}
