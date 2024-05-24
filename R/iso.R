#' Calculate isodistance contours from specified points.
#'
#' Function is fully vectorized to calculate accept vectors of central points
#' and vectors defining multiple isodistances.
#'
#' @param graph `data.frame` or equivalent object representing the network
#' graph. For `dodgr` street networks, this may be a network derived from either
#' \pkg{sf} or \pkg{silicate} ("sc") data, generated with
#' \link{weight_streetnet}.
#' @param from Vector or matrix of points **from** which isodistances are to
#' be calculated.
#' @param dlim Vector of desired limits of isodistances in metres.
#' @param concavity A value between 0 and 1, with 0 giving (generally smoother
#' but less detailed) convex iso-contours and 1 giving highly concave (and
#' generally more detailed) contours.
#' @param length_threshold The minimal length of a segment of the iso-contour
#' to be made more convex according to the 'concavity` parameter.. Low values
#' will produce highly detailed hulls which may cause problems; if in doubt, or
#' if odd results appear, increase this value.
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
                            concavity = 0,
                            length_threshold = 0,
                            contract = TRUE,
                            heap = "BHeap") {

    if (is.null (dlim)) {
        stop ("dlim must be specified")
    }
    if (!is.numeric (dlim)) {
        stop ("dlim must be numeric")
    }

    requireNamespace ("geodist")

    concavity <- check_concavity (concavity)
    # Then adjust to inverse value:
    concavity <- 1 / max (concavity, 1e-6)

    dat <- iso_pre (graph, from, heap, contract = contract)

    d <- m_iso_calculate (dat, dlim)
    from_id <- gsub ("\\_start$", "", dat$from_index$id)

    return (dmat_to_hulls (d, from_id, dat$v, dlim, concavity, length_threshold))
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

iso_calculate <- function (dat, dlim) {

    d <- rcpp_get_iso (
        dat$graph, dat$vert_map, dat$from_index$index,
        dlim, dat$heap
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

    return (d)
}

m_iso_calculate <- memoise::memoise (iso_calculate)

#' Calculate isochrone contours from specified points.
#'
#' Function is fully vectorized to calculate accept vectors of central points
#' and vectors defining multiple isochrone thresholds.
#'
#' @inherit dodgr_isodists
#'
#' @param graph `data.frame` or equivalent object representing the network
#' graph. For `dodgr` street networks, this must be a network derived from
#' \pkg{silicate} ("sc") data, generated with \link{weight_streetnet}. This
#' function does not work with networks derived from \pkg{sf} data.
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
                              concavity = 0,
                              length_threshold = 0,
                              heap = "BHeap") {

    if (!methods::is (graph, "dodgr_streetnet_sc")) {
        stop (
            "isochrones can only be calculated from SC-class networks ",
            "generated by osmdata_sc."
        )
    }
    graph$d_weighted <- graph$time_weighted
    graph$d <- graph$time

    res <- dodgr_isodists (
        graph,
        from = from, dlim = tlim,
        concavity = concavity, length_threshold = length_threshold, heap = heap
    )
    names (res) [names (res) == "dlim"] <- "tlim"
    return (res)
}

#' Calculate isodistance or isochrone contours from specified points.
#'
#' Returns lists of all network vertices contained within the contours. Function
#' is fully vectorized to calculate accept vectors of central points and vectors
#' defining multiple isochrone thresholds. Provide one or more `dlim` values for
#' isodistances, or one or more `tlim` values for isochrones.
#'
#' @inheritParams dodgr_isodists
#'
#' @param graph `data.frame` or equivalent object representing the network
#' graph. For `dodgr` street networks, this must be a network derived from
#' \pkg{silicate} ("sc") data, generated with \link{weight_streetnet}. This
#' function does not work with networks derived from \pkg{sf} data.
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

    d <- m_iso_calculate (dat, dlim)
    from_id <- gsub ("\\_start$", "", dat$from_index$id)

    res <- dmat_to_pts (d, from_id, dat$v, dlim)
    if (has_tlim) {
        names (res) [which (names (res) == "dlim")] <- "tlim"
    }
    return (res)
}

# convert distance matrix with values equal to various isodistances into list of
# lists of points defining the iso-contour hulls ordered around the central points
dmat_to_hulls <- function (d, from, v, dlim, concavity, length_threshold) {

    pt_names <- colnames (d)

    pts <- list ()
    for (i in seq_len (nrow (d))) { # The "from" vertices
        o <- v [match (from [i], v$id), ]
        pts [[i]] <- lapply (dlim, function (j) {
            pts_j <- pt_names [which (d [i, ] <= j)]
            pts_xy <- v [match (pts_j, v$id), c ("x", "y")]
            h0 <- grDevices::chull (pts_xy) - 1L # 0-indexed
            hull <- rcpp_concaveman (pts_xy, h0, concavity, length_threshold)
            res <- NULL
            if (length (hull) > 0) {
                # Then match back to `pts_j`:
                index <- geodist::geodist_min (hull, v [, c ("x", "y")])
                res <- v [c (index, index [1]), ]

                res$from <- o$id
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

    # flatten lists:
    pts <- do.call (rbind, lapply (pts, function (i) do.call (rbind, i)))

    rownames (pts) <- NULL

    return (pts)
}

# convert distance matrix with values equal to various isodistances into list of
# lists of points
dmat_to_pts <- function (d, from, v, dlim) {

    pt_names <- colnames (d)

    pts <- list ()
    for (i in seq_len (nrow (d))) { # The "from" vertices
        o <- v [match (from [i], v$id), ]
        pts [[i]] <- lapply (dlim, function (j) {
            pts_j <- pt_names [which (d [i, ] <= j)]
            res <- v [match (pts_j, v$id), ]
            res$from <- o$id
            res$dlim <- j
            res <- res [, c (
                "from",
                "dlim",
                "id",
                "x",
                "y"
            )]
            return (res)
        })
        names (pts [[i]]) <- paste (dlim)
    }
    names (pts) <- rownames (d)

    # flatten lists:
    pts <- do.call (rbind, lapply (pts, function (i) do.call (rbind, i)))

    rownames (pts) <- NULL

    return (pts)
}

check_concavity <- function (concavity) {
    if (!(is.numeric (concavity) || length (concavity) > 1)) {
        stop ("concavity must be numeric")
    }
    if (concavity < 0 || concavity > 1) {
        message ("concavity must be between 0 and 1; setting to default of 0")
        concavity <- 0
    }
    return (concavity)
}
