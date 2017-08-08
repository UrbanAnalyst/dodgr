#' get_vertices
#'
#' Extract vertices of graph, including spatial coordinates if included
#'
#' @param graph A flat table of graph edges. Must contain columns labelled
#' \code{from} and \code{to}, or \code{start} and \code{stop}. May also contain
#' similarly labelled columns of spatial coordinates (for example
#' \code{from_x}) or \code{stop_lon}).
#' @return A \code{data.frame} of vertices with unique numbers (\code{n}).
#'
#' @note Values of \code{n} are 0-indexed
#'
#' @export
get_vertices <- function (graph)
{
    fr_col <- which (grepl ("fr", names (graph), ignore.case = TRUE) |
                     grepl ("sta", names (graph), ignore.case = TRUE))
    to_col <- which (grepl ("to", names (graph), ignore.case = TRUE) |
                     grepl ("sto", names (graph), ignore.case = TRUE))

    if (ncol (graph) > 4 &
        (any (grepl ("x", names (graph), ignore.case = TRUE)) |
         any (grepl ("y", names (graph), ignore.case = TRUE)) |
         any (grepl ("lon", names (graph), ignore.case = TRUE)) |
         any (grepl ("lat", names (graph), ignore.case = TRUE))))
    {
        # graph is spatial
        if (length (fr_col) != length (to_col))
            stop (paste0 ("from and to columns in graph appear ",
                          "to have different strutures"))
        if (length (fr_col) < 2 | length (to_col) < 2)
            stop (paste0 ("Graph appears to be spatial yet unable to ",
                          "extract coordinates."))

        if (length (fr_col) == 3)
        {
            frx_col <- find_xy_col (graph, fr_col, x = TRUE)
            fry_col <- find_xy_col (graph, fr_col, x = FALSE)
            frid_col <- fr_col [which (!fr_col %in%
                                       c (frx_col, fry_col))]
            fr_col <- c (frx_col, fry_col)
            xy_fr_id <- graph [, frid_col]
            if (!is.character (xy_fr_id))
                xy_fr_id <- paste0 (xy_fr_id)

            tox_col <- find_xy_col (graph, to_col, x = TRUE)
            toy_col <- find_xy_col (graph, to_col, x = FALSE)
            toid_col <- to_col [which (!to_col %in%
                                       c (tox_col, toy_col))]
            to_col <- c (tox_col, toy_col)
            xy_to_id <- graph [, toid_col]
            if (!is.character (xy_to_id))
                xy_to_id <- paste0 (xy_to_id)
        } else # len == 2, so must be only x-y
        {
            xy_fr_id <- paste0 (graph [, fr_col [1]],
                                graph [, fr_col [2]])
            xy_to_id <- paste0 (graph [, to_col [1]],
                                graph [, to_col [2]])
        }

        xy_fr <- graph [, fr_col]
        xy_to <- graph [, to_col]
        if (!(all (apply (xy_fr, 2, is.numeric)) |
              all (apply (xy_to, 2, is.numeric))))
            stop (paste0 ("graph appears to have non-numeric ",
                          "longitudes and latitudes"))

        verts <- data.frame (id = c (xy_fr_id, xy_to_id),
                             x = c (graph [, fr_col [1]],
                                    graph [, to_col [1]]),
                             y = c (graph [, fr_col [2]],
                                    graph [, to_col [2]]),
                             stringsAsFactors = FALSE)
    } else # non-spatial graph
    {
        if (!(length (fr_col) == 1 & length (to_col) == 1))
            stop ("Graph appears to be non-spatial, yet unable to ",
                  "determine vertex columns")

        verts <- data.frame (from = graph [, fr_col],
                             to = graph [, to_col])
    }
    indx <- which (!duplicated (verts))
    verts <- verts [indx, ]
    verts$n <- seq (nrow (verts)) - 1

    return (verts)
}

#' match_pts_to_graph
#'
#' Match spatial points to a spatial graph which contains vertex coordindates
#'
#' @param verts A \code{data.frame} of vertices obtained from
#' \code{get_vertices(graph)}.
#' @param xy coordinates of points to be matched to the vertices
#'
#' @return A vector index into verts
#' @export
match_pts_to_graph <- function (verts, xy)
{
    if (!(is.matrix (xy) | is.data.frame (xy)))
        stop ("xy must be a matrix or data.frame")
    if (ncol (xy) != 2)
        stop ("xy must have only two columns")

    nms <- names (verts)
    if (is.null (nms))
        nms <- colnames (verts)
    if (!is.null (nms))
    {
        ix <- which (grepl ("x", nms, ignore.case = TRUE) |
                     grepl ("lon", nms, ignore.case = TRUE))
        iy <- which (grepl ("y", nms, ignore.case = TRUE) |
                     grepl ("lat", nms, ignore.case = TRUE))
    } else
    {
        message ("xy has no named columns; assuming order is x then y")
        ix <- 1
        iy <- 2
    }

    verts <- data.frame (x = verts [, ix], y = verts [, iy])

    rcpp_points_index (verts, xy)
}
