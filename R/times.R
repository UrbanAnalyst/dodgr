#' dodgr_times
#'
#' Calculate matrix of pair-wise travel times between points.
#'
#' @inherit dodgr_dists
#' @param graph A `dodgr` network returned from the \link{weight_streetnet}
#' function using a network obtained with the \pkg{osmdata} `osmdata_sc`
#' function.
#'
#' @export 
dodgr_times <- function (graph, from = NULL, to = NULL, heap = 'BHeap')
{
    if (attr (graph, "turn_penalty")) # either TRUE or FALSE
    {
        # have to match any from or to points onto junctions that have been
        # re-routed for turn angles. These central junction points are all
        # marked with "_start" (for .vx0) or "_end" (for .vx1).
        v <- unique (c (graph$.vx0, graph$.vx1))
        v <- v [grep ("start|end", v)]
        v <- gsub ("_start", "", v)
        v <- gsub ("_end", "", v)
        v <- unique (v) # list of junction vertices

        from [from %in% v] <- paste0 (from [from %in% v], "_start")
        to [to %in% v] <- paste0 (to [to %in% v], "_end")
    }

    graph$d <- graph$d_weighted <- graph$time
    graph$time <- NULL

    dodgr_dists (graph, from, to, heap = heap)
}
