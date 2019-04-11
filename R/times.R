#' dodgr_times
#'
#' Calculate matrix of pair-wise travel times between points.
#'
#' @inherit dodgr_dists
#' @param graph A `dodgr` network returned from the \link{weight_streetnet}
#' function using a network obtained with the \pkg{osmdata} `osmdata_sc`
#' function.
#' @param ignore_oneway Only if `travel_time = TRUE`, a `TRUE` value allows
#' travel \emph{against} the permitted direction of oneway streets.
#' @param left_side Only if `travel_time = TRUE`, a `TRUE` value indicates
#' traffic that travels on the left side of the street; `FALSE` for the right
#' side.
#' @param turn_penalty Time penalty in seconds for turning across oncoming
#' traffic.
#'
#' @export 
dodgr_times <- function (graph, from = NULL, to = NULL, heap = 'BHeap',
                         ignore_oneway = FALSE, left_side = FALSE,
                         turn_penalty = 45) 
{
    hash <- digest::digest (graph)
    prefix <- "routetimes"
    if (is_graph_cached (hash, prefix))
    {
        graph <- retrieve_cached_graph (hash, prefix)
    } else
    {
        graph <- rm_duplicated_edges (graph)
        res <- rcpp_route_times (graph, ignore_oneway, left_side, turn_penalty)
        graph <- rbind (graph, res$graph)
        graph$d <- graph$d_weighted <- graph$time
        graph$time <- NULL
        cache_graph (graph, hash, prefix)
    }
    dodgr_dists (graph, from, to, heap)
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
