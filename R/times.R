#' dodgr_times
#'
#' Calculate matrix of pair-wise travel times between points.
#'
#' @inherit dodgr_dists
#' @param graph A `dodgr` network returned from the \link{weight_streetnet}
#' function using a network obtained with the \pkg{osmdata} `osmdata_sc`
#' function.
#' @param left_side A `TRUE` value indicates traffic that travels on the left
#' side of the street; `FALSE` for the right side.
#' @param turn_penalty Time penalty in seconds for turning across oncoming
#' traffic.
#'
#' @export 
dodgr_times <- function (graph, from = NULL, to = NULL, heap = 'BHeap',
                         left_side = FALSE, turn_penalty = 10) 
{
    if (is.null (from) | is.null (to))
        stop ("both from and to must be specified.")

    if (turn_penalty > 0)
    {
        attr (graph, "left_side") <- left_side
        attr (graph, "turn_penalty") <- turn_penalty
        hash <- digest::digest (graph)
        prefix <- "routetimes"
        # The cached object is not the graph itself, rather just the new bit routing
        # across junctions according to turn angles:
        if (is_graph_cached (hash, prefix))
        {
            res <- retrieve_cached_graph (hash, prefix)
        } else
        {
            res <- rcpp_route_times (graph, left_side, turn_penalty)
            cache_graph (res, hash, prefix)
        }

        # The junction vertices can still be used as routing points, but need to
        # be disconnected from the replacement turning-angle junctions. This is
        # done by seperately renaming the incoming and outgoing versions:
        index <- which (graph$.vx0 %in% res$junction_vertices)
        v_start <- graph$.vx0 [index]
        graph$.vx0 [index] <- paste0 (graph$.vx0 [index], "_start")
        index <- which (graph$.vx1 %in% res$junction_vertices)
        v_end <- graph$.vx1 [index]
        graph$.vx1 [index] <- paste0 (graph$.vx1 [index], "_end")

        graph <- rbind (graph, res$graph)

        from [from %in% v_start] <- paste0 (from [from %in% v_start], "_start")
        to [to %in% v_end] <- paste0 (to [to %in% v_end], "_end")
    }

    graph$d <- graph$d_weighted <- graph$time
    graph$time <- NULL

    dodgr_dists (graph, from, to, heap = heap)
}

is_graph_cached <- function (hash, prefix)
{
    fname <- file.path (tempdir (), paste0 (prefix, "_", hash, ".Rds"))
    file.exists (fname)
}

# hash is from a reference graph; the cached graph here is the modified version
# with turning angles
cache_graph <- function (graph, hash, prefix)
{
    fname <- file.path (tempdir (), paste0 (prefix, "_", hash, ".Rds"))
    saveRDS (graph, file = fname)
}

retrieve_cached_graph <- function (hash, prefix)
{
    fname <- file.path (tempdir (), paste0 (prefix, "_", hash, ".Rds"))
    readRDS (fname)
}
