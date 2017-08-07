#' get_streetnet
#'
#' @param pts List of points presumably containing spatial coordinates
#' @param expand Relative factor by which street network should extend beyond
#' limits defined by pts
#' @return square matrix of distances between nodes
#'
#' @export
get_streetnet <- function (pts, expand = 0.05)
{
    colx <- which (grepl ("x", names (pts), ignore.case = TRUE) |
                   grepl ("lon", names (pts), ignore.case = TRUE))
    coly <- which (grepl ("y", names (pts), ignore.case = TRUE) |
                   grepl ("lat", names (pts), ignore.case = TRUE))
    if (! (length (colx) == 1 | length (coly) == 1))
        stop ("Can not unambiguous determine coordinates in graph")

    x <- range (pts [, colx])
    x <- x + c (-expand, expand) * diff (x)
    y <- range (pts [, coly])
    y <- y + c (-expand, expand) * diff (y)

    dat <- osmdata::opq (c (x [1], y [1], x [2], y [2])) %>%
        osmdata::add_osm_feature (key = "highway") %>%
        osmdata::osmdata_sf ()

    return (dat$osm_lines)
}

#' weight_streetnet
#'
#' Weight an OSM street network according to a named routino profile, selected
#' from (foot, horse, wheelchair, bicycle, moped, motorcycle, motorcar, goods,
#' hgv, psv).
#'
#' @param net Street network extracted with \code{get_stretnet}
#' @param profile Name of weighting profile
#'
#' @return A \code{data.frame} of edges representing the street network.
#'
#' @export
weight_streetnet <- function (net, profile = "bicycle")
{
    prf_names <- c ("foot", "horse", "wheelchair", "bicycle", "moped",
                    "motorcycle", "motorcar", "goods", "hgv", "psv")
    profile = match.arg (tolower (profile), prf_names)
    profiles <- dodgr::weighting_profiles
    pr <- profiles [profiles$name == profile, ]

    dat <- rcpp_lines_as_network (net, pr = pr)
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
