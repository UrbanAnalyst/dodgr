#' @title Calculate matrix of pair-wise distances between points.
#'
#' @description Calculates distances from input `data.frame` objects (`graph`),
#' which must minimally contain three columns of `from`, `to`, and `d` or
#' `dist`. If an additional column named `weight` or `wt` is present, shortest
#' paths are calculated according to values specified in that column, while
#' distances returned are calculated from the `d` or `dist` column. That is,
#' paths between any pair of points will be calculated according to the minimal
#' total sum of `weight` values (if present), while reported distances will be
#' total sums of `dist` values.
#'
#' Graphs derived from Open Street Map street networks, via the
#' \link{weight_streetnet} function, have columns labelled `d`, `d_weighted`,
#' `time`, and `time_weighted`. For these inputs, paths between origin and
#' destination points are always routed using `d_weighted` (or `t_weighted` for
#' times), while final distances are sums of values of `d` (or `t` for times)-
#' that is, of un-weighted distances or times - along those paths.
#'
#' The function is parallelized for efficient computation of distances between
#' multiple origin and destination points, as described in the `from` parameter
#' below.
#'
#' @param graph `data.frame` or equivalent object representing the network
#' graph (see Notes). For `dodgr` street networks, this may be a network derived
#' from either \pkg{sf} or \pkg{silicate} ("sc") data, generated with
#' \link{weight_streetnet}.
#'
#' The `from` and `to` columns of `graph` may be either single
#' columns of numeric or character values specifying the numbers or names of
#' graph vertices, or combinations to two columns specifying geographical
#' (longitude and latitude,) coordinates. In the latter case, almost any sensible
#' combination of names will be accepted (for example, `fromx, fromy`,
#' `from_x, from_y`, or `fr_lat, fr_lon`.)
#'
#' Note that longitude and latitude values are always interpreted in 'dodgr' to
#' be in EPSG:4326 / WSG84 coordinates. Any other kinds of coordinates should
#' first be reprojected to EPSG:4326 before submitting to any 'dodgr' routines.
#'
#' See further information in Details.
#'
#' @param from Vector or matrix of points **from** which route distances are to
#' be calculated, specified as one of the following:
#' \itemize{
#' \item Single character vector precisely matching node numbers or names
#' given in `graph$from` or `graph$to`.
#' \item Single vector of integer-ish values, in which case these will be
#' presumed to specify indices into \link{dodgr_vertices}, and NOT to
#' correspond to values in the 'from' or 'to' columns of the graph. See the
#' example below for a demonstration.
#' \item Matrix or equivalent of longitude and latitude coordinates, in which
#' case these will be matched on to the nearest coordinates of 'from' and 'to'
#' points in the graph.
#' }
#'
#' @param to Vector or matrix of points **to** which route distances are to be
#' calculated. If `to` is `NULL`, pairwise distances will be calculated from
#' all `from` points to all other nodes in `graph`. If both `from` and `to` are
#' `NULL`, pairwise distances are calculated between all nodes in `graph`.
#'
#' @param shortest If `FALSE`, calculate distances along the \emph{fastest}
#' rather than shortest routes. For street networks produced with
#' \link{weight_streetnet}, distances may also be calculated along the
#' \emph{fastest} routes with the `shortest = FALSE` option. Graphs must in
#' this case have columns of `time` and `time_weighted`. Note that the fastest
#' routes will only be approximate when derived from \pkg{sf}-format data
#' generated with the \pkg{osmdata} function `osmdata_sf()`, and will be much
#' more accurate when derived from `sc`-format data generated with
#' `osmdata_sc()`. See \link{weight_streetnet} for details.
#'
#' @param pairwise If `TRUE`, calculate distances only between the ordered
#' pairs of `from` and `to`.
#'
#' @param heap Type of heap to use in priority queue. Options include
#' Fibonacci Heap (default; `FHeap`), Binary Heap (`BHeap`),
#' `Trinomial Heap (`TriHeap`), Extended Trinomial Heap
#' (`TriHeapExt`, and 2-3 Heap (`Heap23`).
#'
#' @param parallel If `TRUE`, perform routing calculation in parallel.
#' Calculations in parallel ought very generally be advantageous. For small
#' graphs, calculating distances in parallel is likely to offer relatively
#' little gain in speed, but increases from parallel computation will generally
#' markedly increase with increasing graph sizes. By default, parallel
#' computation uses the maximal number of available cores or threads. This
#' number can be reduced by specifying a value via
#' `RcppParallel::setThreadOptions (numThreads = <desired_number>)`. Parallel
#' calculations are, however, not able to be interrupted (for example, by
#' `Ctrl-C`), and can only be stopped by killing the R process.
#'
#' @param quiet If `FALSE`, display progress messages on screen.
#'
#' @return square matrix of distances between nodes
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
#' # Example of "from" and "to" as integer-ish values, in which case they are
#' # interpreted to index into "dodgr_vertices()":
#' graph <- data.frame (
#'     from = c (1, 3, 2, 2, 3, 3, 4, 4),
#'     to = c (2, 1, 3, 4, 2, 4, 3, 1),
#'     d = c (1, 2, 1, 3, 2, 1, 2, 1)
#' )
#' dodgr_dists (graph, from = 1, to = 2)
#' # That then gives distance from "1" to "3" because the vertices are built
#' # sequentially along "graph$from":
#' dodgr_vertices (graph)
#' # And vertex$id [2] is "3"
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
#' # https://github.com/UrbanAnalyst/dodgr/issues/47
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
dodgr_dists <- function (graph,
                         from = NULL,
                         to = NULL,
                         shortest = TRUE,
                         pairwise = FALSE,
                         heap = "BHeap",
                         parallel = TRUE,
                         quiet = TRUE) {

    graph <- tbl_to_df (graph)

    hps <- get_heap (heap, graph)
    heap <- hps$heap
    graph <- hps$graph

    graph <- preprocess_spatial_cols (graph)
    gr_cols <- dodgr_graph_cols (graph)
    is_spatial <- is_graph_spatial (graph)
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

    if (parallel && heap == "TriHeapExt") {
        if (!quiet) {
            message (
                "Extended TriHeaps can not be calculated in parallel; ",
                "reverting to serial calculation"
            )
        }
        parallel <- FALSE
    }

    d <- calculate_distmat (
        graph,
        to_from_indices$vert_map,
        to_from_indices$from,
        to_from_indices$to,
        heap,
        is_spatial,
        parallel,
        pairwise
    )


    if (!quiet) {
        message ("done.")
    }

    return (d)
}

#' Calculate matrix of pair-wise distances between points.
#'
#' Alias for \link{dodgr_dists}
#' @inherit dodgr_dists
#' @family distances
#' @export
dodgr_distances <- function (graph,
                             from = NULL,
                             to = NULL,
                             shortest = TRUE,
                             pairwise = FALSE,
                             heap = "BHeap",
                             parallel = TRUE,
                             quiet = TRUE) {

    dodgr_dists (graph,
        from,
        to,
        shortest = shortest,
        pairwise = pairwise,
        heap = heap,
        parallel = parallel,
        quiet = quiet
    )
}

#' Get an index of `pts` matching `vert_map`.
#'
#' Also returns the corresonding names of those `pts`
#'
#' @return list of `index`, which is 0-based for C++, and corresponding
#' `id` values.
#' @noRd
get_index_id_cols <- function (graph,
                               gr_cols,
                               vert_map,
                               pts) {

    index <- -1
    id <- NULL
    if (!missing (pts)) {
        if (is.integer (pts) && is.vector (pts)) {
            index <- pts
        } else if (is.character (pts) ||
            is.numeric (pts) ||
            is.matrix (pts) ||
            is.data.frame (pts)) {
            index <- get_pts_index (graph, gr_cols, vert_map, pts)
            if ((is.matrix (pts) || is.data.frame (pts)) &&
                !any (duplicated (index))) {
                rownames (pts) <- vert_map$vert [index]
            }
        } else {
            stop (
                "routing points are of unknown form; must be either ",
                "character, matrix, or integer"
            )
        }

        if (length (pts == 2) && is.numeric (pts) &&
            ((any (grepl ("x", names (pts), ignore.case = TRUE)) &&
                any (grepl ("y", names (pts), ignore.case = TRUE))) ||
                (any (grepl ("lon", names (pts), ignore.case = TRUE) &&
                    (any (grepl ("lat", names (pts), ignore.case = TRUE))))))) {
            names (pts) <- NULL
        }

        id <- get_id_cols (pts)
        if (is.null (id)) {
            id <- vert_map$vert [index]
        } # index is 1-based
    }
    list (index = index, id = id)
}


#' Get the ID columns or rownames of from or to points
#'
#' @param pts The `from` or `to` args passed to `dodgr_dists`.
#' @return Character vector of names of points, if they exist in `pts`
#' @noRd
get_id_cols <- function (pts) {

    ids <- NULL
    if (any (grepl ("id", colnames (pts), ignore.case = TRUE))) {
        nmc <- which (grepl ("id", colnames (pts)))
        if (methods::is (pts, "data.frame")) {
            ids <- pts [[nmc]]
        } else if (is.matrix (pts)) {
            ids <- pts [, nmc, drop = TRUE]
        }
    } else if (is.vector (pts) && !is.null (names (pts))) {
        ids <- names (pts)
    } else if (!is.null (rownames (pts))) {
        ids <- rownames (pts)
    }
    return (ids)
}

#' Map unique vertex names to sequential numbers in matrix
#'
#' @noRd
make_vert_map <- function (graph,
                           gr_cols,
                           xy = FALSE) {

    # gr_cols are (edge_id, from, to, d, w, component, xfr, yfr, xto, yto)
    verts <- c (paste0 (graph [[gr_cols$from]]), paste0 (graph [[gr_cols$to]]))
    indx <- which (!duplicated (verts))
    if (!xy) {
        # Note id has to be 0-indexed:
        res <- data.frame (
            vert = paste0 (verts [indx]),
            id = seq (indx) - 1,
            stringsAsFactors = FALSE
        )
    } else {
        verts_x <- c (graph [[gr_cols$xfr]], graph [[gr_cols$xto]])
        verts_y <- c (graph [[gr_cols$yfr]], graph [[gr_cols$yto]])
        res <- data.frame (
            vert = paste0 (verts [indx]),
            id = seq (indx) - 1,
            x = verts_x [indx],
            y = verts_y [indx],
            stringsAsFactors = FALSE
        )
    }
    return (res)
}

#' Convert `from` or `to` args of \link{dodgr_dists} to indices into
#' \link{dodgr_vertices}.
#'
#' @param graph A dodgr graph
#' @param vert_map Two-column `data.frame` of unique vertices and
#' corresponding IDs, obtained from `make_vert_map`
#' @param gr_cols Returned from `dodgr_graph_cols()`
#' @param pts Either a vector of names, or a matrix or `data.frame` of
#' arbitrary geographical coordinates for which to get index into vertices of
#' graph.
#'
#' @noRd
get_pts_index <- function (graph,
                           gr_cols,
                           vert_map,
                           pts) {

    if (!(is.matrix (pts) || is.data.frame (pts))) {
        if (!is.numeric (pts)) {
            pts <- matrix (pts, ncol = 1)
        } else {
            nms <- names (pts)
            if (is.null (nms) && length (pts) > 1L) {
                nms <- c ("x", "y")
            }
            pts <- matrix (pts, nrow = 1) # vector of (x,y) vals
            colnames (pts) <- nms
        }
    }

    if (ncol (pts) == 1) {

        pts <- get_pts_index_vec (pts, vert_map)

    } else {

        pts <- get_pts_index_rect (pts, graph, gr_cols)
    }

    pts
}

get_pts_index_vec <- function (pts, vert_map) {

    pts <- pts [, 1]

    if (!is.numeric (pts)) {

        indx <- match (pts, vert_map$vert)

        if (any (is.na (indx))) {
            stop (paste0 (
                "from/to are not numeric yet can not be",
                " matched onto graph vertices"
            ))
        }
        pts <- indx
    }

    if (any (pts < 1 | pts > nrow (vert_map))) {
        stop (paste0 ("points exceed numbers of vertices"))
    }

    return (pts)
}

get_pts_index_rect <- function (pts, graph, gr_cols) {

    nms <- names (pts)
    if (is.null (nms)) {
        nms <- colnames (pts)
    }

    ix <- which (grepl ("x", nms, ignore.case = TRUE) |
        grepl ("lon", nms, ignore.case = TRUE))
    iy <- which (grepl ("y", nms, ignore.case = TRUE) |
        grepl ("lat", nms, ignore.case = TRUE))

    if (length (ix) != 1 || length (iy) != 1) {
        stop (paste0 (
            "Unable to determine geographical ",
            "coordinates in from/to"
        ))
    }

    index <- match (c ("xfr", "yfr", "xto", "yto"), names (gr_cols))
    if (any (is.na (gr_cols [index]))) {
        stop (paste0 (
            "Cannot determine geographical coordinates ",
            "against which to match pts"
        ))
    }

    if (is.data.frame (pts)) {
        names (pts) [ix] <- "x"
        names (pts) [iy] <- "y"
    } else {
        colnames (pts) [ix] <- "x"
        colnames (pts) [iy] <- "y"
    }

    # Result of rcpp_points_index is 0-indexed for C++
    pts <- rcpp_points_index_par (dodgr_vertices (graph), pts) + 1

    return (pts)
}

# nocov start

#' Download a street network when `graph` not passed to `dodgr_dists`,
#' by using the lists of from and to points.
#'
#' @param from Arg passed to `dodgr_dists`
#' @param to Arg passed to `dodgr_dists`
#' @param expand Factor by which street network is to be expanded beyond range
#' of `from` and `to` points.
#' @return Converted graph as `data.frame`
#' @noRd
graph_from_pts <- function (from,
                            to,
                            expand = 0.1,
                            wt_profile = "bicycle",
                            quiet = TRUE) {

    if (!quiet) {
        message (
            paste0 (
                "No graph submitted to dodgr_dists; ",
                "downloading street network ... "
            ),
            appendLF = FALSE
        )
    }

    pts <- NULL
    if (!missing (from)) {
        pts <- from
    }
    if (!missing (to)) {
        pts <- rbind (pts, to)
    }
    pts <- pts [which (!duplicated (pts)), ]
    graph <- dodgr_streetnet (pts = pts, expand = expand) %>%
        weight_streetnet (wt_profile = wt_profile)

    if (!quiet) {
        message ("done")
    }

    return (graph)
}

# nocov end

#' Flip from and two vertices of a graph.
#'
#' This is only called from `dodgr_distances`, and only on converted graph, so
#' just needs to switch from and to vertex columns.
#' @noRd
flip_graph <- function (graph) {

    grcols <- dodgr_graph_cols (graph)
    fr_temp <- graph [[grcols$from]]
    graph [[grcols$from]] <- graph [[grcols$to]]
    graph [[grcols$to]] <- fr_temp

    return (graph)
}

#' Call the actual C++ functions to calculate and return distance matrices
#' @noRd
calculate_distmat <- function (graph,
                               vert_map,
                               from_index,
                               to_index,
                               heap,
                               is_spatial,
                               parallel = TRUE,
                               pairwise = FALSE) {

    flip <- FALSE
    if (length (from_index$index) > length (to_index$index)) {
        flip <- TRUE
        graph <- flip_graph (graph)
        temp <- from_index
        from_index <- to_index
        to_index <- temp
    }

    if (parallel) {
        if (pairwise) {
            d <- rcpp_get_sp_dists_paired_par (
                graph,
                vert_map,
                from_index$index,
                to_index$index,
                heap,
                is_spatial
            )
        } else {
            d <- rcpp_get_sp_dists_par (
                graph,
                vert_map,
                from_index$index,
                to_index$index,
                heap,
                is_spatial
            )
        }
    } else {
        d <- rcpp_get_sp_dists (
            graph,
            vert_map,
            from_index$index,
            to_index$index,
            heap
        )
    }

    if (!pairwise) {
        if (!is.null (from_index$id)) {
            rownames (d) <- from_index$id
        } else {
            rownames (d) <- vert_map$vert
        }
        if (!is.null (to_index$id)) {
            colnames (d) <- to_index$id
        } else {
            colnames (d) <- vert_map$vert
        }

        if (get_turn_penalty (graph) > 0) {

            rownames (d) <- gsub ("\\_(start|end)$", "", rownames (d))
            colnames (d) <- gsub ("\\_(start|end)$", "", colnames (d))
        }

        if (flip) {
            d <- t (d)
        }
    }

    return (d)
}
