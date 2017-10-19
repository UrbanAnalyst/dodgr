#' dodgr_streetnet
#'
#' Use the \code{osmdata} package to extract the street network for a given
#' location. For routing between a given set of points (passed as \code{pts}),
#' the \code{bbox} argument may be omitted, in which case a bounding box will
#' be constructed by expanding the range of \code{pts} by the relative amount of
#' \code{expand}.
#'
#' @param bbox Bounding box as vector or matrix of coordinates, or location
#' name. Passed to \code{osmdata::getbb}.
#' @param pts List of points presumably containing spatial coordinates
#' @param expand Relative factor by which street network should extend beyond
#' limits defined by pts (only if \code{bbox} not given).
#' @param quiet If \code{FALSE}, display progress messages
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
dodgr_streetnet <- function (bbox, pts, expand = 0.05, quiet = TRUE)
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

    net <- osmdata::opq (bbox) %>%
                osmdata::add_osm_feature (key = "highway") %>%
                osmdata::osmdata_sf (quiet = quiet) %>%
                extract2 ("osm_lines")
    if (nrow (net) == 0)
        stop ("Street network unable to be downloaded")

    return (net)
}

#' weight_streetnet
#'
#' Weight (or re-weight) an \code{sf}-formatted OSM street network according to
#' a named routino profile, selected from (foot, horse, wheelchair, bicycle,
#' moped, motorcycle, motorcar, goods, hgv, psv).
#'
#' @param sf_lines A street network represented as \code{sf} \code{LINESTRING}
#' objects, typically extracted with \code{get_streetnet}
#' @param wt_profile Name of weighting profile, or vector of values with names
#' corresponding to names in \code{type_col}
#' @param type_col Specify column of the \code{sf} \code{data.frame} object
#' which designates different types of highways to be used for weighting
#' (default works with \code{osmdata} objects).
#' @param id_col Specify column of the code{sf} \code{data.frame} object which
#' provides unique identifiers for each highway (default works with
#' \code{osmdata} objects).
#'
#' @return A \code{data.frame} of edges representing the street network, along
#' with a column of graph component numbers.
#'
#' @export
#' @examples
#' # hampi is included with package as an 'osmdata' sf-formatted street network
#' net <- weight_streetnet (hampi)
#' class(net) # data.frame
#' dim(net) # 6096  11; 6096 streets
#' # os_roads_bristol is also included as an sf data.frame, but in a different
#' # format requiring identification of columns and specification of custom
#' # weighting scheme.
#' colnm <- "formOfWay"
#' wts <- c (0.1, 0.2, 0.8, 1)
#' names (wts) <- unique (os_roads_bristol [[colnm]])
#' net <- weight_streetnet (os_roads_bristol, wt_profile = wts,
#'                          type_col = colnm, id_col = "identifier")
#' dim (net) # 406 11; 406 streets
weight_streetnet <- function (sf_lines, wt_profile = "bicycle",
                              type_col = "highway", id_col = "osm_id")
{
    if (!is (sf_lines, "sf"))
        stop ('sf_lines must be class "sf"')
    if (!all (c ("geometry", type_col, id_col) %in% names (sf_lines)))
        stop (paste0 ('sf_lines must be class "sf" and ',
                      'have highway and geometry columns'))

    if (type_col != "highway")
        names (sf_lines) [which (names (sf_lines) == type_col)] <- "highway"
    if (id_col != "osm_id")
        names (sf_lines) [which (names (sf_lines) == id_col)] <- "osm_id"

    if (is.null (names (sf_lines$geometry)))
        names (sf_lines$geometry) <- sf_lines$osm_id

    if (is.character (wt_profile))
    {
        prf_names <- c ("foot", "horse", "wheelchair", "bicycle", "moped",
                        "motorcycle", "motorcar", "goods", "hgv", "psv")
        wt_profile <- match.arg (tolower (wt_profile), prf_names)
        profiles <- dodgr::weighting_profiles
        wt_profile <- profiles [profiles$name == wt_profile, ]
        wt_profile$value <- wt_profile$value / 100
    } else if (is.numeric (wt_profile) & !is.null (names (wt_profile)))
    {
        nms <- names (wt_profile)
        wt_profile <- data.frame (name = "custom",
                                  way = nms,
                                  value = wt_profile)
    } else
        stop ("Custom named profiles must be vectors with named values")


    dat <- rcpp_sf_as_network (sf_lines, pr = wt_profile)
    graph <- data.frame (edge_id = seq (nrow (dat$character_values)),
                         from_id = as.character (dat$character_values [, 1]),
                         from_lon = dat$numeric_values [, 1],
                         from_lat = dat$numeric_values [, 2],
                         to_id = as.character (dat$character_values [, 2]),
                         to_lon = dat$numeric_values [, 3],
                         to_lat = dat$numeric_values [, 4],
                         d = dat$numeric_values [, 5],
                         d_weighted = dat$numeric_values [, 6],
                         highway = as.character (dat$character_values [, 3]),
                         way_id = as.character (dat$character_values [, 4]),
                         stringsAsFactors = FALSE
                         )

    # get component numbers for each edge
    class (graph) <- c (class (graph), "dodgr_streetnet")
    graph <- dodgr_components (graph)

    return (graph)
}
