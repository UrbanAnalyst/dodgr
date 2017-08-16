#' dodgr_streetnet
#'
#' Use the \code{osmdata} package to extract the street network for a given
#' location. For routing between a given set of points (passed as \code{pts}),
#' the \code{bbox} argument may be ommitted, in which case a bounding box will
#' be constructed by expanding the range of \code{pts} by the relative amount of
#' \code{expand}.
#'
#' @param bbox Bounding box as vector or matrix of coordinates, or location
#' name. Passed to \code{osmdata::getbb}.
#' @param pts List of points presumably containing spatial coordinates
#' @param expand Relative factor by which street network should extend beyond
#' limits defined by pts (only if \code{bbox} not given).
#' @return A Simple Features (\code{sf}) object with coordinates of all lines in
#' the street network.
#'
#' @export
dodgr_streetnet <- function (bbox, pts, expand = 0.05)
{
    if (!missing (bbox))
        bbox <- osmdata::getbb (bbox)
    else if (!missing (pts))
    {
        nms <- names (pts)
        if (is.null (nms))
            nms <- colnames (pts)
        colx <- which (grepl ("x", nms, ignore.case = TRUE) |
                       grepl ("lon", nms, ignore.case = TRUE))
        coly <- which (grepl ("y", nms, ignore.case = TRUE) |
                       grepl ("lat", nms, ignore.case = TRUE))
        if (! (length (colx) == 1 | length (coly) == 1))
            stop ("Can not unambiguous determine coordinates in graph")

        x <- range (pts [, colx])
        x <- x + c (-expand, expand) * diff (x)
        y <- range (pts [, coly])
        y <- y + c (-expand, expand) * diff (y)

        bbox <- c (x [1], y [1], x [2], y [2])
    } else
        stop ('Either bbox or pts must be specified.')

    dat <- osmdata::opq (bbox) %>%
        osmdata::add_osm_feature (key = "highway") %>%
        osmdata::osmdata_sf ()

    return (dat$osm_lines)
}

#' weight_streetnet
#'
#' Weight (or re-weight) an \code{sf}-formatted OSM street network according to
#' a named routino profile, selected from (foot, horse, wheelchair, bicycle,
#' moped, motorcycle, motorcar, goods, hgv, psv).
#'
#' @param graph Street network extracted with \code{get_stretnet}
#' @param wt_profile Name of weighting profile
#'
#' @return A \code{data.frame} of edges representing the street network, along
#' with a column of graph component numbers.
#'
#' @export
weight_streetnet <- function (graph, wt_profile = "bicycle")
{
    if (!is (graph, "sf"))
        stop ('graph must be class "sf"')
    if (!all (c ("geometry", "highway", "osm_id") %in% names (graph)))
        stop ('graph must be class "sf" and have highway and geometry columns')

    prf_names <- c ("foot", "horse", "wheelchair", "bicycle", "moped",
                    "motorcycle", "motorcar", "goods", "hgv", "psv")
    wt_profile <- match.arg (tolower (wt_profile), prf_names)
    profiles <- dodgr::weighting_profiles
    wt_profile <- profiles [profiles$name == wt_profile, ]
    wt_profile$value <- wt_profile$value / 100

    dat <- rcpp_sf_as_network (graph, pr = wt_profile)
    data.frame (edge_id = seq (nrow (dat [[1]])),
                from_id = as.character (dat [[2]] [, 1]),
                from_lon = dat [[1]] [, 1],
                from_lat = dat [[1]] [, 2],
                to_id = as.character (dat [[2]] [, 2]),
                to_lon = dat [[1]] [, 3],
                to_lat = dat [[1]] [, 4],
                d = dat [[1]] [, 5],
                d_weighted = dat [[1]] [, 6],
                highway = as.character (dat [[2]] [, 3]),
                stringsAsFactors = FALSE
                )
}
