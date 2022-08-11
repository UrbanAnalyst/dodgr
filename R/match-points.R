#' match_pts_to_verts
#'
#' Match spatial points to the vertices of a spatial graph
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
#' @family misc
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
    if (!(is.matrix (xy) || is.data.frame (xy))) {
        stop ("xy must be a matrix or data.frame")
    }
    if (!is (xy, "sf")) {
        if (ncol (xy) != 2) {
            stop ("xy must have only two columns")
        }
    }

    indx <- seq (nrow (verts))
    if (connected) {
        vertsi <- verts [which (verts$component == 1), ]
        indx <- match (vertsi$id, verts$id)
    }

    xyi <- find_xy_col_simple (verts)
    verts <- data.frame (x = verts [indx, xyi [1]], y = verts [indx, xyi [2]])
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

    # rcpp_points_index is 0-indexed, so ...
    indx [rcpp_points_index_par (verts, xy) + 1L]
}

#' match_points_to_verts
#'
#' Alias for \link{match_pts_to_verts}
#' @inherit match_pts_to_verts
#' @family misc
#' @export
match_points_to_verts <- function (verts, xy, connected = FALSE) {

    match_pts_to_verts (verts, xy, connected = connected)
}
