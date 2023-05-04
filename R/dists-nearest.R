#' Calculate vector of shortest distances from a series of 'from' points to
#' nearest one of series of 'to' points.
#'
#' @param graph `data.frame` or equivalent object representing the network
#' graph (see Notes)
#' @param from Vector or matrix of points **from** which route distances are to
#' be calculated (see Notes)
#' @param to Vector or matrix of points **to** which shortest route distances
#' are to be calculated to nearest 'to' point only.
#' @param shortest If `FALSE`, calculate distances along the \emph{fastest}
#' rather than shortest routes (see Notes).
#' @param heap Type of heap to use in priority queue. Options include
#' Fibonacci Heap (default; `FHeap`), Binary Heap (`BHeap`),
#' `Trinomial Heap (`TriHeap`), Extended Trinomial Heap
#' (`TriHeapExt`, and 2-3 Heap (`Heap23`).
#' @param quiet If `FALSE`, display progress messages on screen.
#' @return Vector of distances, one element for each 'from' point giving the
#' distance to the nearest 'to' point.
#'
#' @note `graph` must minimally contain three columns of `from`,
#' `to`, `dist`. If an additional column named `weight` or
#' `wt` is present, shortest paths are calculated according to values
#' specified in that column; otherwise according to `dist` values. Either
#' way, final distances between `from` and `to` points are calculated
#' by default according to values of `dist`. That is, paths between any pair of
#' points will be calculated according to the minimal total sum of `weight`
#' values (if present), while reported distances will be total sums of `dist`
#' values.
#'
#' For street networks produced with \link{weight_streetnet}, distances may also
#' be calculated along the \emph{fastest} routes with the `shortest = FALSE`
#' option. Graphs must in this case have columns of `time` and `time_weighted`.
#' Note that the fastest routes will only be approximate when derived from
#' \pkg{sf}-format data generated with the \pkg{osmdata} function
#' `osmdata_sf()`, and will be much more accurate when derived from `sc`-format
#' data generated with `osmdata_sc()`. See \link{weight_streetnet} for details.
#'
#' The `from` and `to` columns of `graph` may be either single
#' columns of numeric or character values specifying the numbers or names of
#' graph vertices, or combinations to two columns specifying geographical
#' (longitude and latitude) coordinates. In the latter case, almost any sensible
#' combination of names will be accepted (for example, `fromx, fromy`,
#' `from_x, from_y`, or `fr_lat, fr_lon`.)
#'
#' `from` and `to` values can be either two-column matrices or
#' equivalent of longitude and latitude coordinates, or else single columns
#' precisely matching node numbers or names given in `graph$from` or
#' `graph$to`. If `to` is `NULL`, pairwise distances are calculated from all
#' `from` points to all other nodes in `graph`. If both `from` and `to` are
#' `NULL`, pairwise distances are calculated between all nodes in `graph`.
#'
#' Calculations are always calculated in parallel, using multiple threads.
#'
#' @family distances
#' @export
#' @examples
#' # A simple graph
#' graph <- data.frame (
#'     from = c ("A", "B", "B", "B", "C", "C", "D", "D"),
#'     to = c ("B", "A", "C", "D", "B", "D", "C", "A"),
#'     d = c (1, 2, 1, 3, 2, 1, 2, 1)
#' )
#' dodgr_dists (graph)
#'
#' # A larger example from the included [hampi()] data.
#' graph <- weight_streetnet (hampi)
#' from <- sample (graph$from_id, size = 100)
#' to <- sample (graph$to_id, size = 50)
#' d <- dodgr_dists (graph, from = from, to = to)
#' # d is a 100-by-50 matrix of distances between `from` and `to`
#'
#' \dontrun{
#' # a more complex street network example, thanks to @chrijo; see
#' # https://github.com/ATFutures/dodgr/issues/47
#'
#' xy <- rbind (
#'     c (7.005994, 51.45774), # limbeckerplatz 1 essen germany
#'     c (7.012874, 51.45041)
#' ) # hauptbahnhof essen germany
#' xy <- data.frame (lon = xy [, 1], lat = xy [, 2])
#' essen <- dodgr_streetnet (pts = xy, expand = 0.2, quiet = FALSE)
#' graph <- weight_streetnet (essen, wt_profile = "foot")
#' d <- dodgr_dists (graph, from = xy, to = xy)
#' # First reason why this does not work is because the graph has multiple,
#' # disconnected components.
#' table (graph$component)
#' # reduce to largest connected component, which is always number 1
#' graph <- graph [which (graph$component == 1), ]
#' d <- dodgr_dists (graph, from = xy, to = xy)
#' # should work, but even then note that
#' table (essen$level)
#' # There are parts of the network on different building levels (because of
#' # shopping malls and the like). These may or may not be connected, so it may
#' # be necessary to filter out particular levels
#' index <- which (!(essen$level == "-1" | essen$level == "1")) # for example
#' library (sf) # needed for following sub-select operation
#' essen <- essen [index, ]
#' graph <- weight_streetnet (essen, wt_profile = "foot")
#' graph <- graph [which (graph$component == 1), ]
#' d <- dodgr_dists (graph, from = xy, to = xy)
#' }
dodgr_dists_nearest <- function (graph,
                                 from = NULL,
                                 to = NULL,
                                 shortest = TRUE,
                                 heap = "BHeap",
                                 quiet = TRUE) {

    graph <- tbl_to_df (graph)

    hps <- get_heap (heap, graph)
    heap <- hps$heap
    graph <- hps$graph

    graph <- preprocess_spatial_cols (graph)
    gr_cols <- dodgr_graph_cols (graph)
    to_from_indices <- to_from_index_with_tp (graph, from, to)
    if (to_from_indices$compound) {
        graph <- to_from_indices$graph_compound
    }

    if (!shortest) {
        if (is.na (gr_cols$time_weighted)) {
            stop (
                "Graph does not contain a weighted time column from ",
                "which to calculate fastest paths."
            )
        }
        graph [[gr_cols$d_weighted]] <- graph [[gr_cols$time_weighted]]
    }

    graph <- convert_graph (graph, gr_cols)

    if (!quiet) {
        message ("Calculating shortest paths ... ", appendLF = FALSE)
    }

    if (heap == "TriHeapExt") {
        stop (
            "Extended TriHeaps can not be calculated in parallel.",
            call. = FALSE
        )
    }

    d <- rcpp_get_sp_dists_nearest (
        graph,
        to_from_indices$vert_map,
        to_from_indices$from$index,
        to_from_indices$to$index,
        heap
    )
    index <- seq_along (to_from_indices$from$index)
    nearest_index <- as.integer (d [index + length (index)])
    d <- d [index]
    nearest_index <- match (nearest_index, to_from_indices$to$index)
    nearest_ids <- to_from_indices$to$id [nearest_index]
    nearest_ids <- gsub ("\\_(start|end)$", "", nearest_ids)

    return (data.frame (
        from = to_from_indices$from$id,
        to = nearest_ids,
        d = d,
        stringsAsFactors = FALSE
    ))
}
