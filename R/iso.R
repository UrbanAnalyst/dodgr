#' dodgr_isodists
#'
#' Calculate isodistance contours from specified points. Function is fully
#' vectorized to calculate accept vectors of central points and vectors
#' defining multiple isodistances.
#'
#' @param graph `data.frame` or equivalent object representing the network
#' graph (see Notes)
#' @param from Vector or matrix of points **from** which isodistances are to
#' be calculated.
#' @param dlim Vector of desired limits of isodistances in metres.
#' @param contract If `TRUE`, calculate isodists only to vertices in the
#' contract graph, in other words, only to junction vertices.
#' @param heap Type of heap to use in priority queue. Options include
#' Fibonacci Heap (default; `FHeap`), Binary Heap (`BHeap`),
#' Trinomial Heap (`TriHeap`), Extended Trinomial Heap
#' (`TriHeapExt`, and 2-3 Heap (`Heap23`).
#' @return A single `data.frame` of isodistances as points sorted anticlockwise
#' around each origin (`from`) point, with columns denoting the `from` points
#' and `dlim` value(s). The isodistance contours are given as `id` values and
#' associated coordinates of the series of points from each `from` point at the
#' specified isodistances.
#'
#' @note Isodists are calculated by default using parallel computation with the
#' maximal number of available cores or threads. This number can be reduced by
#' specifying a value via
#' `RcppParallel::setThreadOptions (numThreads = <desired_number>)`.
#'
#' @family distances
#' @export
#' @examples
#' graph <- weight_streetnet (hampi)
#' from <- sample (graph$from_id, size = 100)
#' dlim <- c (1, 2, 5, 10, 20) * 100
#' d <- dodgr_isodists (graph, from = from, dlim)
dodgr_isodists <- function (graph,
                            from = NULL,
                            dlim = NULL,
                            contract = TRUE,
                            heap = "BHeap") {

    if (is.null (dlim)) {
        stop ("dlim must be specified")
    }
    if (!is.numeric (dlim)) {
        stop ("dlim must be numeric")
    }

    dat <- iso_pre (graph, from, heap, contract = contract)

    # expand dlim to an extra value to capture max boundary
    if (length (dlim) == 1) {
        dlim_exp <- c (dlim, dlim + dlim / 4) # nocov
    } else {
        dlim_exp <- c (dlim, dlim [length (dlim)] + dlim [length (dlim)] -
            dlim [length (dlim) - 1])
    }

    d <- rcpp_get_iso (
        dat$graph, dat$vert_map, dat$from_index$index,
        sort (dlim_exp), dat$heap
    )
    d [d > max (dlim)] <- NA

    vert_names <- gsub ("\\_(start|end)$", "", dat$vert_map$vert)
    from_id <- gsub ("\\_start$", "", dat$from_index$id)

    if (!is.null (dat$from_index$id)) {
        rownames (d) <- from_id
    } else {
        rownames (d) <- vert_names
    } # nocov
    colnames (d) <- vert_names

    # verts ON isohulls are flagged with negative values at isodistance;
    # terminal verts are also negative at their specified distance; all other
    # points inside any hulls are positive. The latter are removed here:
    d [d >= 0] <- NA
    d <- -d # convert negative-flagged iso-values to positive

    return (dmat_to_pts (d, from_id, dat$v, dlim))
}

iso_pre <- function (graph, from = NULL, heap = "BHeap", contract = TRUE) {

    v <- dodgr_vertices (graph)
    graph <- tbl_to_df (graph)

    hps <- get_heap (heap, graph)
    heap <- hps$heap
    graph <- hps$graph

    if (!is.null (from)) {
        from <- nodes_arg_to_pts (from, graph)
    }

    graph <- preprocess_spatial_cols (graph)
    gr_cols <- dodgr_graph_cols (graph)
    to_from_indices <- to_from_index_with_tp (graph, from, to = NULL)
    if (to_from_indices$compound) {
        graph <- to_from_indices$graph_compound
    }

    if (contract && !methods::is (graph, "dodgr_contracted")) {
        graph <- dodgr_contract_graph (graph, verts = to_from_indices$from$id)
    }

    graph <- convert_graph (graph, gr_cols)

    list (
        v = v,
        graph = graph,
        vert_map = to_from_indices$vert_map,
        from_index = to_from_indices$from,
        heap = heap
    )
}

#' dodgr_isochrones
#'
#' Calculate isochrone contours from specified points. Function is fully
#' vectorized to calculate accept vectors of central points and vectors
#' defining multiple isochrone thresholds.
#'
#' @inherit dodgr_isodists
#'
#' @param from Vector or matrix of points **from** which isochrones are to
#' be calculated.
#' @param tlim Vector of desired limits of isochrones in seconds
#' @return A single `data.frame` of isochrones as points sorted anticlockwise
#' around each origin (`from`) point, with columns denoting the `from` points
#' and `tlim` value(s). The isochrones are given as `id` values and associated
#' coordinates of the series of points from each `from` point at the specified
#' isochrone times.
#'
#' Isochrones are calculated by default using parallel computation with the
#' maximal number of available cores or threads. This number can be reduced by
#' specifying a value via `RcppParallel::setThreadOptions (numThreads =
#' <desired_number>)`.
#'
#' @family distances
#' @export
#' @examples
#' \dontrun{
#' # Use osmdata package to extract 'SC'-format data:
#' library (osmdata)
#' dat <- opq ("hampi india") %>%
#'     add_osm_feature (key = "highway") %>%
#'     osmdata_sc ()
#' graph <- weight_streetnet (dat)
#' from <- sample (graph$.vx0, size = 100)
#' tlim <- c (5, 10, 20, 30, 60) * 60 # times in seconds
#' x <- dodgr_isochrones (graph, from = from, tlim)
#' }
dodgr_isochrones <- function (graph,
                              from = NULL,
                              tlim = NULL,
                              heap = "BHeap") {

    if (!methods::is (graph, "dodgr_streetnet_sc")) {
        stop (
            "isochrones can only be calculated from SC-class networks ",
            "generated by osmdata_sc."
        )
    }
    graph$d_weighted <- graph$time_weighted
    graph$d <- graph$time

    res <- dodgr_isodists (graph, from = from, dlim = tlim, heap = heap)
    names (res) [names (res) == "dlim"] <- "tlim"
    return (res)
}

#' dodgr_isoverts
#'
#' Calculate isodistance or isochrone contours from specified points, and return
#' lists of all network vertices contained within the contours. Function is
#' fully vectorized to calculate accept vectors of central points and vectors
#' defining multiple isochrone thresholds. Provide one or more `dlim` values for
#' isodistances, or one or more `tlim` values for isochrones.
#'
#' @inheritParams dodgr_isodists
#'
#' @param from Vector or matrix of points **from** which isodistances or
#' isochrones are to be calculated.
#' @param tlim Vector of desired limits of isochrones in seconds
#' @return A single `data.frame` of vertex IDs, with columns denoting the `from`
#' points and `tlim` value(s). The isochrones are given as `id` values and
#' associated coordinates of the series of points from each `from` point at the
#' specified isochrone times.
#'
#' Isoverts are calculated by default using parallel computation with the
#' maximal number of available cores or threads. This number can be reduced by
#' specifying a value via `RcppParallel::setThreadOptions (numThreads =
#' <desired_number>)`.
#'
#' @family distances
#' @export
#' @examples
#' \dontrun{
#' # Use osmdata package to extract 'SC'-format data:
#' library (osmdata)
#' dat <- opq ("hampi india") %>%
#'     add_osm_feature (key = "highway") %>%
#'     osmdata_sc ()
#' graph <- weight_streetnet (dat)
#' from <- sample (graph$.vx0, size = 100)
#' tlim <- c (5, 10, 20, 30, 60) * 60 # times in seconds
#' x <- dodgr_isoverts (graph, from = from, tlim)
#' }
dodgr_isoverts <- function (graph,
                            from = NULL,
                            dlim = NULL,
                            tlim = NULL,
                            heap = "BHeap") {

    if (!methods::is (graph, "dodgr_streetnet_sc")) {
        stop (
            "isoverts can only be calculated from SC-class networks ",
            "generated by osmdata_sc."
        )
    }
    if (!is.null (dlim) && !is.null (tlim)) {
        stop ("Only one of dlim or tlim can be provided")
    }

    has_tlim <- FALSE
    if (!is.null (tlim)) {
        graph$d_weighted <- graph$time_weighted
        graph$d <- graph$time
        dlim <- tlim
        has_tlim <- TRUE
    }

    dat <- iso_pre (graph, from, heap)

    # expand dlim to an extra value to capture max boundary
    if (length (dlim) == 1) {
        dlim_exp <- c (dlim, dlim + dlim / 4)
    } else {
        dlim_exp <- c (dlim, dlim [length (dlim)] + dlim [length (dlim)] -
            dlim [length (dlim) - 1])
    }

    d <- rcpp_get_iso (
        dat$graph,
        dat$vert_map,
        dat$from_index$index,
        sort (dlim_exp),
        dat$heap
    )
    d [d > max (dlim)] <- NA
    index <- which (!is.na (d))
    d [index] [d [index] < 0] <- -d [index] [d [index] < 0]

    if (!is.null (dat$from_index$id)) {
        rownames (d) <- gsub ("\\_start$", "", dat$from_index$id)
    } else {
        rownames (d) <- gsub ("\\_start$", "", dat$vert_map$vert)
    }
    colnames (d) <- gsub ("\\_(start|end)$", "", dat$vert_map$vert)

    # convert d-values to the *next highest* specified dlim value. breaks need
    # to start < 0 to include 0 in lowest class
    breaks <- c (-0.001, dlim)
    na_index <- which (!is.na (d))
    f <- cut (d [na_index], breaks = breaks)
    index <- match (f, attr (f, "levels"))
    d [na_index] <- breaks [-1] [index]

    res <- dmat_to_pts (d, dat$from_index$id, dat$v, dlim)
    if (has_tlim) {
        names (res) [names (res) == "dlim"] <- "tlim"
    }
    return (res)
}

# convert distance matrix with values equal to various isodistances into list of
# lists of points ordered around the central points
dmat_to_pts <- function (d, from, v, dlim) {

    pt_names <- colnames (d)
    pts <- list ()
    for (i in seq_len (nrow (d))) {
        o <- v [match (from [i], v$id), ]
        pts [[i]] <- lapply (dlim, function (j) {
            res <- pt_names [which (d [i, ] == j)]
            # Then any additional terminal vertices
            index <- which (!d [i, ] %in% dlim &
                d [i, ] < j)
            res <- c (res, pt_names [index])
            res <- v [match (res, v$id), ]
            if (nrow (res) > 0) {
                res$from <- o$id
                res <- order_points (res, o)
                res$dlim <- j
                res <- res [, c (
                    "from",
                    "dlim",
                    "id",
                    "x",
                    "y"
                )]
            }
            return (res)
        })
        names (pts [[i]]) <- paste (dlim)
    }
    names (pts) <- rownames (d)

    # flatten lists
    pts <- do.call (rbind, lapply (pts, function (i) do.call (rbind, i)))
    rownames (pts) <- NULL

    return (pts)
}

# order points around circle
order_points <- function (pts, origin) {

    dx <- pts$x - origin$x
    dy <- pts$y - origin$y
    theta <- rep (NA, nrow (pts))

    index <- which (dx > 0 & dy >= 0)
    theta [index] <- atan (dy [index] / dx [index])
    index <- which (dx > 0 & dy < 0)
    theta [index] <- 2 * pi + atan (dy [index] / dx [index])
    index <- which (dx < 0)
    theta [index] <- pi + atan (dy [index] / dx [index])
    index <- which (dx == 0 & dy >= 0)
    theta [index] <- 0
    index <- which (dx == 0 & dy < 0)
    theta [index] <- 3 * pi / 2

    pts <- pts [order (theta), c ("from", "id", "x", "y")]
    rbind (pts, pts [1, ])
}
