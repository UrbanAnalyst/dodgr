null_to_na <- function (x) {

    if (length (x) == 0) {
        x <- NA
    }
    return (x)
}

#' Match the heap arg and convert graph is necessary
#'
#' @param heap Name of heap as passed to `dodgr_dists`
#' @param graph `data.frame` of graph edges
#' @return List of matched heap arg and potentially converted graph
#' @noRd
get_heap <- function (heap,
                      graph) {

    heaps <- c ("FHeap", "BHeap", "TriHeap", "TriHeapExt", "Heap23", "set")
    heap <- match.arg (arg = heap, choices = heaps)

    list (heap = heap, graph = graph)
}

#' Get appropriate measure for geodist distances.
#'
#' Default measure is "cheap", but that becomes inaccurate beyond around 100km.
#' This function works out the approximate maximal graph distances, and
#' determines an appropriate measure based on that. Note that "geodesic"
#' distances are not used, as calculation times for those are enormously longer
#' than either cheap or Haversine.
#'
#' @return "cheap" if maximal distances are < 100km, otherwise "haversine".
#' @noRd
get_geodist_measure <- function (graph) {

    dmax <- max_spatial_dist (graph) / 1000
    ifelse (dmax < 100, "cheap", "haversine")
}
