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
#' @examples
#' \dontrun{
#' streetnet <- dodgr_streetnet ("hampi india", expand = 0)
#' # convert to form needed for \code{dodgr} functions:
#' graph <- weight_streetnet (streetnet)
#' nrow (graph) # 5,742 edges
#' # Alternative ways of extracting street networks by using a small selection of
#' # graph vertices to define bounding box:
#' verts <- dodgr_vertices (graph)
#' verts <- verts [sample (nrow (verts), size = 200), ]
#' streetnet <- dodgr_streetnet (pts = verts, expand = 0)
#' graph <- weight_streetnet (streetnet)
#' nrow (graph)
#' # This will generally have many more rows because most street networks include
#' # streets that extend considerably beyond the specified bounding box.
#' }
dodgr_streetnet <- function (bbox, pts, expand = 0.05)
{
    if (!missing (bbox))
    {
        bbox <- osmdata::getbb (bbox)
        bbox [1, ] <- bbox [1, ] + c (-expand, expand) * diff (bbox [1, ])
        bbox [2, ] <- bbox [2, ] + c (-expand, expand) * diff (bbox [2, ])
    }
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
#' @param sf_lines A street network represented as \code{sf} \code{LINESTRING}
#' objects, typically extracted with \code{get_stretnet}
#' @param wt_profile Name of weighting profile
#'
#' @return A \code{data.frame} of edges representing the street network, along
#' with a column of graph component numbers.
#'
#' @export
#' @examples
#' net <- weight_streetnet (hampi) # internal sf-formatted street network
#' class(net) # data.frame
#' dim(net) # 6096  11; 6096 streets
weight_streetnet <- function (sf_lines, wt_profile = "bicycle")
{
    if (!is (sf_lines, "sf"))
        stop ('sf_lines must be class "sf"')
    if (!all (c ("geometry", "highway", "osm_id") %in% names (sf_lines)))
        stop (paste0 ('sf_lines must be class "sf" and ',
                      'have highway and geometry columns'))

    prf_names <- c ("foot", "horse", "wheelchair", "bicycle", "moped",
                    "motorcycle", "motorcar", "goods", "hgv", "psv")
    wt_profile <- match.arg (tolower (wt_profile), prf_names)
    profiles <- dodgr::weighting_profiles
    wt_profile <- profiles [profiles$name == wt_profile, ]
    wt_profile$value <- wt_profile$value / 100

    dat <- rcpp_sf_as_network (sf_lines, pr = wt_profile)
    graph <- data.frame (edge_id = seq (nrow (dat [[1]])),
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

    # get component numbers for each edge
    graph$component <- dodgr_components (graph)$component

    return (graph)
}
