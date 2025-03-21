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
#' Measures for graphs are stored in `options("dodgr_dist_measure")`, as a list
#' with each measure named after the graph hash.
#'
#' @return "cheap" if maximal distances are < 100km, otherwise "haversine".
#' @noRd
get_geodist_measure <- function (graph) {

    hash <- attr (graph, "hash")
    measure_list <- getOption ("dodgr_dist_measure", "")

    has_measure <- !is.null (hash)
    has_single_measure <- FALSE
    if ("all" %in% names (measure_list)) {
        has_single_measure <- TRUE
    } else if (has_measure) {
        has_measure <- any (nzchar (measure_list)) && hash %in% names (measure_list)
    }

    if (has_single_measure) {
        measure <- measure_list [["all"]]
    } else if (has_measure) {
        measure <- measure_list [[hash]]
    } else {
        dmax <- max_spatial_dist (graph) / 1000
        measure <- ifelse (dmax < 100, "cheap", "haversine")

        # This is also called at the start of SC construction, before graph has
        # any hash.
        if (!is.null (hash)) {
            if (!any (nzchar (measure_list))) {
                measure_list <- NULL
            }
            measure_list <- c (measure_list, measure)
            names (measure_list) [length (measure_list)] <- eval (hash)
            options ("dodgr_dist_measure" = measure_list)
        }
    }

    return (measure)
}

#' Force \link{weight_streetnet} to use geodesic distances.
#'
#' Distances by default are Mapbox "cheap" distances if maximal network
#' distances are < 100km, otherwise Haversine distances. Calling this function
#' forces all calls to \link{weight_streetnet} from that point on to use
#' geodesic distances. These are more computationally expensive to calculate,
#' and weighting networks will likely take more time.
#'
#' @param unset Calling this function with `unset = TRUE` reverts distance
#' calculations to those described above, rather than geodesic.
#' @return Nothing; the function is called for its side-effect only of setting
#' distance calculations to geodesic.
#'
#' @family extraction
#' @export
dodgr_streetnet_geodesic <- function (unset = FALSE) {

    if (unset) {
        options ("dodgr_dist_measure" = NULL)
        msg <- "revert to default measures"
    } else {
        options ("dodgr_dist_measure" = c (all = "geodesic"))
        msg <- "use the geodesic measure"
    }

    objs <- ls (envir = .GlobalEnv)
    objs_are_graphs <- vapply (objs, function (o) {
        inherits (get (o), "dodgr_streetnet")
    }, logical (1L))
    if (any (objs_are_graphs)) {
        message (
            "Only graphs created from this point on with ",
            "'weight_streetnet()' will ",
            msg
        )
    }
}
