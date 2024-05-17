#' Calculate matrix of pair-wise travel times between points.
#'
#' @inherit dodgr_dists
#' @param graph `data.frame` or equivalent object representing the network
#' graph (see Notes). For `dodgr` street networks, this may be a network derived
#' from either \pkg{sf} or \pkg{silicate} ("sc") data, generated with
#' \link{weight_streetnet}. Note, however, that networks derived from \pkg{sf}
#' data will generally not produce reliable estimates of times. Accurate
#' estimates can only be guaranteed by using networks derived from
#' \pkg{silicate} ("sc") data.
#' @param shortest If `TRUE`, calculate times along the \emph{shortest} rather
#' than fastest paths.
#'
#' @family distances
#' @export
dodgr_times <- function (graph,
                         from = NULL,
                         to = NULL,
                         shortest = FALSE,
                         heap = "BHeap") {

    graph <- tbl_to_df (graph)

    gr_cols <- dodgr_graph_cols (graph)
    if (is.na (gr_cols$time)) {
        stop ("graph has no time column")
    }

    graph [[gr_cols$d]] <- graph [[gr_cols$time]]

    if (!shortest) {
        if (is.na (gr_cols$time_weighted)) {
            stop (
                "Graph does not contain a weighted time column from ",
                "which to calculate fastest paths."
            )
        }
        graph [[gr_cols$d_weighted]] <- graph [[gr_cols$time_weighted]]
    }

    dodgr_dists (graph, from, to, heap = heap)
}
