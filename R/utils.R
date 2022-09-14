null_to_na <- function (x) {

    if (length (x) == 0) {
        x <- NA
    }
    return (x)
}

#' get_heap
#'
#' Match the heap arg and convert graph is necessary
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
