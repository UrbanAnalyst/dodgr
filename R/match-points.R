#' match_pts_to_verts
#'
#' Match spatial points to the vertices of a spatial graph. The
#' \link{match_pts_to_graph} function matches points to the nearest edge based
#' on geometric intersections; this function only matches to the nearest vertex
#' based on point-to-point distances.
#'
#' @param verts A `data.frame` of vertices obtained from
#' `dodgr_vertices(graph)`.
#' @param xy coordinates of points to be matched to the vertices, either as
#' matrix or \pkg{sf}-formatted `data.frame`.
#' @param connected Should points be matched to the same (largest) connected
#' component of graph? If `FALSE` and these points are to be used for a
#' `dodgr` routing routine (\link{dodgr_dists}, \link{dodgr_paths}, or
#' \link{dodgr_flows_aggregate}), then results may not be returned if points are
#' not part of the same connected component. On the other hand, forcing them to
#' be part of the same connected component may decrease the spatial accuracy of
#' matching.
#'
#' @return A vector index into verts
#' @family match
#' @export
#' @examples
#' net <- weight_streetnet (hampi, wt_profile = "foot")
#' verts <- dodgr_vertices (net)
#' # Then generate some random points to match to graph
#' npts <- 10
#' xy <- data.frame (
#'     x = min (verts$x) + runif (npts) * diff (range (verts$x)),
#'     y = min (verts$y) + runif (npts) * diff (range (verts$y))
#' )
#' pts <- match_pts_to_verts (verts, xy)
#' pts # an index into verts
#' pts <- verts$id [pts]
#' pts # names of those vertices
match_pts_to_verts <- function (verts, xy, connected = FALSE) {

    if (!all (c ("id", "x", "y") %in% names (verts))) {
        message (
            "First argument to match_pts_to_verts should be result of ",
            "dodgr_vertices;\npresuming you've submitted the network ",
            "itself and will now try extracting the vertices"
        )
        verts <- dodgr_vertices (verts)
    }

    indx <- seq (nrow (verts))
    if (connected) {
        vertsi <- verts [which (verts$component == 1), ]
        indx <- match (vertsi$id, verts$id)
    }

    xyi <- find_xy_col_simple (verts)
    verts <- data.frame (x = verts [indx, xyi [1]], y = verts [indx, xyi [2]])

    xy <- pre_process_xy (xy)

    # rcpp_points_index is 0-indexed, so ...
    indx [rcpp_points_index_par (verts, xy) + 1L]
}

#' match_points_to_verts
#'
#' Alias for \link{match_pts_to_verts}
#' @inherit match_pts_to_verts
#' @family match
#' @export
match_points_to_verts <- function (verts, xy, connected = FALSE) {

    match_pts_to_verts (verts, xy, connected = connected)
}

pre_process_xy <- function (xy) {

    if (!(is.matrix (xy) || is.data.frame (xy))) {
        stop ("xy must be a matrix or data.frame")
    }
    if (!is (xy, "sf")) {
        if (ncol (xy) != 2) {
            stop ("xy must have only two columns")
        }
    }

    if (is (xy, "tbl")) {
        xy <- data.frame (xy)
    }
    if (is (xy, "sf")) {
        if (!"geometry" %in% names (xy)) {
            stop ("xy has no sf geometry column")
        } # nocov
        if (!is (xy$geometry, "sfc_POINT")) {
            stop ("xy$geometry must be a collection of sfc_POINT objects")
        }
        xy <- unlist (lapply (xy$geometry, as.numeric)) %>%
            matrix (nrow = 2) %>%
            t ()
        xy <- data.frame (x = xy [, 1], y = xy [, 2])
    } else {
        xyi <- find_xy_col_simple (xy)
        xy <- data.frame (x = xy [, xyi [1]], y = xy [, xyi [2]])
    }

    return (xy)
}

#' match_pts_to_graph
#'
#' Match spatial points to the edges of a spatial graph, through finding the
#' edge with the closest perpendicular intersection. NOTE: Intersections are
#' calculated geometrically, and presume planar geometry. It is up to users of
#' projected geometrical data, such as those within a `dodgr_streetnet` object,
#' to ensure that either: (i) Data span an sufficiently small area that errors
#' from presuming planar geometry may be ignored; or (ii) Data are re-projected
#' to an equivalent planar geometry prior to calling this routine.
#'
#' @param graph A `dodgr` graph with spatial coordinates, such as a
#' `dodgr_streetnet` object.
#' @param distances If `TRUE`, return a 'data.frame' object with 'index' column
#' as described in return value; and additional 'dist' column with perpendicular
#' distance to nearest edge in graph. See description of return value for
#' details.
#' @inheritParams match_pts_to_verts
#'
#' @return For 'distances = FALSE' (default), a vector index matching the `xy`
#' coordinates to nearest edges. For bi-directional edges, only one match is
#' returned, and it is up to the user to identify and suitably process matching
#' edge pairs. For 'distances = TRUE', a 'data.frame' of two columns:
#' \itemize{
#' \item "index" The index of closest edges in "graph", as described above.
#' \item "d_signed" The perpendicular distance from ech point to the nearest
#' edge, with negative distances denoting points to the left of edges, and
#' positive distances denoting points to the right. Distances of zero denote
#' points lying precisely on the line of an edge (potentially including cases
#' where nearest point of bisection lies beyond the actual edge).
#' }
#' @family match
#' @export
#' @examples
#' graph <- weight_streetnet (hampi, wt_profile = "foot")
#' # Then generate some random points to match to graph
#' verts <- dodgr_vertices (graph)
#' npts <- 10
#' xy <- data.frame (
#'     x = min (verts$x) + runif (npts) * diff (range (verts$x)),
#'     y = min (verts$y) + runif (npts) * diff (range (verts$y))
#' )
#' edges <- match_pts_to_graph (graph, xy)
#' graph [edges, ] # The edges of the graph closest to `xy`
match_pts_to_graph <- function (graph, xy,
                                connected = FALSE, distances = FALSE) {

    if (!is_graph_spatial (graph)) {
        stop ("Points may only be matched to spatial graphs.")
    }

    if (connected) {
        graph <- graph [which (graph$component == 1), ]
    }

    xy <- pre_process_xy (xy)

    gr_cols <- dodgr_graph_cols (graph)
    gr_cols <- unlist (gr_cols [which (!is.na (gr_cols))])
    graph <- graph [, gr_cols]
    names (graph) <- names (gr_cols)

    res <- rcpp_points_to_edges_par (graph, xy)
    index <- seq (nrow (xy))

    # rcpp_points_index is 0-indexed, so ...
    graph_index <- as.integer (res [index]) + 1L

    if (distances) {
        ret <- data.frame (
            index = graph_index,
            d_signed = signed_intersection_dists (graph, xy, res)
        )
    } else {
        ret <- graph_index
    }

    return (ret)
}

#' match_points_to_graph
#'
#' Alias for \link{match_pts_to_graph}
#' @inherit match_pts_to_graph
#' @family match
#' @export
match_points_to_graph <- function (graph, xy, connected = FALSE) {

    match_pts_to_graph (graph, xy, connected = connected)
}

#' Get geodesic distances to intersection points for match_pts_to_graph.
#'
#' @param res Output of 'rcpp_points_index' function
#' @noRd
signed_intersection_dists <- function (graph, xy, res) {

    n <- nrow (xy)
    index <- seq (n)

    # rcpp_points_index is 0-indexed, so ...
    graph_index <- as.integer (res [index]) + 1L

    # coordinates not yet used here; see #103
    xy_intersect <- data.frame (
        x = res [index + nrow (xy)],
        y = res [index + 2L * nrow (xy)]
    )

    d <- geodist::geodist (
        xy,
        xy_intersect,
        paired = TRUE,
        measure = "geodesic"
    )

    # Then coordinates of graph edges for sign calculation
    xy_cols <- c ("xfr", "yfr", "xto", "yto")
    gxy <- graph [graph_index, xy_cols]

    which_side <- (gxy$xto - gxy$xfr) * (xy_intersect$y - gxy$yfr) -
        (gxy$yto - gxy$yfr) * (xy_intersect$x - gxy$xfr)
    which_side [which_side < 0.0] <- -1
    which_side [which_side > 0.0] <- 1

    return (d * which_side)
}

#' Insert new nodes into a graph, breaking edges at point of nearest
#' intersection.
#'
#' The "id" value of each edge to be divided through insertion of new points is
#' modified to produce two new "id" values with suffixes "_A" and "_B". This
#' routine presumes graphs to be `dodgr_streetnet` object, with geographical
#' coordinates.
#'
#' @inheritParams match_pts_to_graph
#' @return A modified version of `graph`, with additional edges formed by
#' breaking previous edges at nearest penpendicular intersections with the
#' points, `xy`.
#' @family match
#' @examples
#' graph <- weight_streetnet (hampi, wt_profile = "foot")
#' dim (graph)
#'
#' verts <- dodgr_vertices (graph)
#' set.seed (2)
#' npts <- 10
#' xy <- data.frame (
#'     x = min (verts$x) + runif (npts) * diff (range (verts$x)),
#'     y = min (verts$y) + runif (npts) * diff (range (verts$y))
#' )
#'
#' graph <- add_nodes_to_graph (graph, xy)
#' dim (graph) # more edges than original
#' @export
add_nodes_to_graph <- function (graph, xy) {

    index <- match_pts_to_graph (graph, xy)

    gr_cols <- dodgr_graph_cols (graph)
    gr_cols <- unlist (gr_cols [which (!is.na (gr_cols))])
    graph_std <- graph [, gr_cols] # standardise column names
    names (graph_std) <- names (gr_cols)

    # Expand index to include all potentially duplicated edges:
    index <- lapply (seq_along (index), function (i) {
        out <- which (
            (graph_std$from == graph_std$from [index [i]] &
                graph_std$to == graph_std$to [index [i]]) |
                (graph_std$from == graph_std$to [index [i]] &
                    graph_std$to == graph_std$from [index [i]])
        )
        cbind (rep (i, length (out)), out)
    })
    index <- data.frame (do.call (rbind, index))
    names (index) <- c ("n", "index")

    genhash <- function (len = 10) {
        paste0 (sample (c (0:9, letters, LETTERS), size = len), collapse = "")
    }

    edges_to_split <- graph_std [index$index, ]
    graph_to_add <- graph [index$index, ]

    graph_std <- graph_std [-index$index, ]
    graph <- graph [-index$index, ]

    edges_to_split$n <- index$n

    edges_split <- lapply (unique (index$n), function (i) {

        edges_to_split_i <- edges_to_split [which (edges_to_split$n == i), ]

        d_wt <- edges_to_split_i$d_weighted / edges_to_split_i$d
        t_wt <- edges_to_split_i$time_weighted / edges_to_split_i$time
        t_scale <- edges_to_split_i$time / edges_to_split_i$d

        new_edges_i <- lapply (seq (nrow (edges_to_split_i)), function (e) {

            edge_i <- rbind (edges_to_split_i [e, ], edges_to_split_i [e, ])
            edge_i$to [1] <- edge_i$from [2] <- genhash ()
            edge_i$xto [1] <- xy$x [i]
            edge_i$yto [1] <- xy$y [i]
            edge_i$xfr [2] <- xy$x [i]
            edge_i$yfr [2] <- xy$y [i]

            xy_i <- data.frame (
                x = c (edge_i [1, "xfr"], edge_i [1, "xto"], edge_i [2, "xto"]),
                y = c (edge_i [1, "yfr"], edge_i [1, "yto"], edge_i [2, "yto"])
            )
            dmat <- geodist::geodist (xy_i)
            edge_i$d [1] <- dmat [1, 2]
            edge_i$d [2] <- dmat [2, 3]

            edge_i$d_weighted <- edge_i$d * d_wt
            edge_i$time <- edge_i$d * t_scale
            edge_i$time_weighted <- edge_i$time * t_wt

            edge_i$edge_id <- paste0 (edge_i$edge_id, "_", LETTERS [e])

            return (edge_i)
        })

        return (do.call (rbind, new_edges_i))
    })

    edges_split <- do.call (rbind, edges_split)

    # Then match edges_split back on to original graph:
    graph_to_add <- graph_to_add [edges_split$n, ]
    gr_cols <- gr_cols [which (!is.na (gr_cols))]
    for (g in seq_along (gr_cols)) {
        graph_to_add [, gr_cols [g]] <- edges_split [[names (gr_cols) [g]]]
    }

    return (rbind (graph, graph_to_add))
}
