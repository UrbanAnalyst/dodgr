# nocov start

#' dodgr_streetnet
#'
#' Use the `osmdata` package to extract the street network for a given
#' location. For routing between a given set of points (passed as `pts`),
#' the `bbox` argument may be omitted, in which case a bounding box will
#' be constructed by expanding the range of `pts` by the relative amount of
#' `expand`.
#'
#' @param bbox Bounding box as vector or matrix of coordinates, or location
#' name. Passed to `osmdata::getbb`.
#' @param pts List of points presumably containing spatial coordinates
#' @param expand Relative factor by which street network should extend beyond
#' limits defined by pts (only if `bbox` not given).
#' @param quiet If `FALSE`, display progress messages
#' @return A Simple Features (`sf`) object with coordinates of all lines in
#' the street network.
#'
#' @note Calls to this function may return "General overpass server error" with
#' a note that "Query timed out." The overpass served used to access the data
#' has a sophisticated queueing system which prioritises requests that are
#' likely to require little time. These timeout errors can thus generally *not*
#' be circumvented by changing "timeout" options on the HTTP requests, and
#' should rather be interpreted to indicate that a request is too large, and may
#' need to be refined, or somehow broken up into smaller queries.
#'
#' @family extraction
#' @export
#' @examples
#' \dontrun{
#' streetnet <- dodgr_streetnet ("hampi india", expand = 0)
#' # convert to form needed for `dodgr` functions:
#' graph <- weight_streetnet (streetnet)
#' nrow (graph) # around 5,900 edges
#' # Alternative ways of extracting street networks by using a small selection
#' # of graph vertices to define bounding box:
#' verts <- dodgr_vertices (graph)
#' verts <- verts [sample (nrow (verts), size = 200), ]
#' streetnet <- dodgr_streetnet (pts = verts, expand = 0)
#' graph <- weight_streetnet (streetnet)
#' nrow (graph)
#' # This will generally have many more rows because most street networks
#' # include streets that extend considerably beyond the specified bounding box.
#'
#' # bbox can also be a polygon:
#' bb <- osmdata::getbb ("gent belgium") # rectangular bbox
#' nrow (dodgr_streetnet (bbox = bb)) # around 30,000
#' bb <- osmdata::getbb ("gent belgium", format_out = "polygon")
#' nrow (dodgr_streetnet (bbox = bb)) # around 17,000
#' # The latter has fewer rows because only edges within polygon are returned
#'
#' # Example with access restrictions
#' bbox <- c (-122.2935, 47.62663, -122.28, 47.63289)
#' x <- dodgr_streetnet_sc (bbox)
#' net <- weight_streetnet (x, keep_cols = "access", turn_penalty = TRUE)
#' # has many streets with "access" = "private"; these can be removed like this:
#' net2 <- net [which (!net$access != "private"), ]
#' # or modified in some other way such as strongly penalizing use of those
#' # streets:
#' index <- which (net$access == "private")
#' net$time_weighted [index] <- net$time_weighted [index] * 100
#' }
dodgr_streetnet <- function (bbox,
                             pts = NULL,
                             expand = 0.05,
                             quiet = TRUE) {

    bb <- process_bbox (bbox, pts, expand)

    # osm_poly2line merges all street polygons with the line ones
    net <- osmdata::opq (bb$bbox) %>%
        osmdata::add_osm_feature (key = "highway") %>%
        osmdata::osmdata_sf (quiet = quiet) %>%
        osmdata::osm_poly2line ()

    if (nrow (net$osm_lines) == 0) {
        stop ("Street network unable to be downloaded")
    }

    if (!is.null (bb$bbox_poly)) {
        net <- osmdata::trim_osmdata (net, bb$bbox_poly)
    }

    return (net$osm_lines)
}

#' dodgr_streetnet_sc
#'
#' Use the `osmdata` package to extract the street network for a given
#' location and return it in `SC`-format. For routing between a given set of
#' points (passed as `pts`), the `bbox` argument may be omitted, in which case a
#' bounding box will be constructed by expanding the range of `pts` by the
#' relative amount of `expand`.
#'
#' @family extraction
#' @inherit dodgr_streetnet
#' @export
dodgr_streetnet_sc <- function (bbox,
                                pts = NULL,
                                expand = 0.05,
                                quiet = TRUE) {

    bb <- process_bbox (bbox, pts, expand)

    # This extracts only some of the most common restrictions; see full list at
    # https://wiki.openstreetmap.org/wiki/Restrictions
    fts <- c (
        "\"highway\"",
        "\"restriction\"",
        "\"access\"",
        "\"bicycle\"",
        "\"foot\"",
        "\"motorcar\"",
        "\"motor_vehicle\"",
        "\"vehicle\"",
        "\"toll\""
    )

    osmdata::opq (bb$bbox) %>%
        osmdata::add_osm_features (features = fts) %>%
        osmdata::osmdata_sc (quiet = quiet)
}

# nocov end

process_bbox <- function (bbox,
                          pts = NULL,
                          expand = 0.05) {

    bbox_poly <- NULL

    if (!missing (bbox)) {

        if (is.character (bbox)) {

            bbox <- osmdata::getbb (bbox) # nocov

        } else if (is.list (bbox)) {

            if (!all (vapply (bbox, is.numeric, logical (1)))) {
                stop (
                    "bbox is a list, so items must be numeric ",
                    "(as in osmdata::getbb (..., format_out = 'polygon'))"
                )
            }
            if (length (bbox) > 1) {
                message ("selecting the first polygon from bbox")
            } # nocov

            bbox_poly <- bbox [[1]]
            colnames (bbox_poly) <- c ("x", "y")
            bbox <- apply (bbox [[1]], 2, range)
            colnames (bbox) <- c ("x", "y")
            rownames (bbox) <- c ("min", "max")

        } else if (is.numeric (bbox)) {

            if (!inherits (bbox, "matrix")) {

                if (length (bbox) != 4) {
                    stop ("bbox must have four numeric values")
                }
                bbox <- vec_to_bbox (bbox)

            } else if (nrow (bbox) > 2) {

                bbox_poly <- bbox
                colnames (bbox_poly) <- c ("x", "y")
                bbox <- apply (bbox, 2, range)
            }

            # if x/y are rows, then transpose to columns:
            ptn <- "^x|x$|^y|y$|^lon|lon$|long$|longitude|^lat|lat$|latitude"
            if (all (grepl (ptn, rownames (bbox)))) {
                bbox <- t (bbox)
            }

            rownames (bbox) <- c ("min", "max")
            colnames (bbox) <- c ("x", "y")
        }

        if (identical (rownames (bbox), c ("x", "y")) ||
            identical (colnames (bbox), c ("min", "max"))) {
            bbox <- t (bbox)
        }

        bbox [, 1] <- bbox [, 1] + c (-expand, expand) * diff (bbox [, 1])
        bbox [, 2] <- bbox [, 2] + c (-expand, expand) * diff (bbox [, 2])

    } else if (!is.null (pts)) {

        nms <- names (pts)
        if (is.null (nms)) {
            nms <- colnames (pts)
        }

        colx <- which (grepl ("x", nms, ignore.case = TRUE) |
            grepl ("lon", nms, ignore.case = TRUE))
        coly <- which (grepl ("y", nms, ignore.case = TRUE) |
            grepl ("lat", nms, ignore.case = TRUE))

        if (!(length (colx) == 1 || length (coly) == 1)) {
            stop ("Can not unambiguously determine coordinates in graph")
        }

        x <- range (pts [, colx])
        x <- x + c (-expand, expand) * diff (x)
        y <- range (pts [, coly])
        y <- y + c (-expand, expand) * diff (y)

        bbox <- cbind (
            c (x [1], x [2]),
            c (y [1], y [2])
        )
        rownames (bbox) <- c ("min", "max")
        colnames (bbox) <- c ("x", "y")

    } else {
        stop ("Either bbox or pts must be specified.")
    }

    # columns are then [x, y]:
    if (bbox [1, 1] < -180) {
        bbox [1, 1] <- bbox [1, 2] + 360
    }
    if (bbox [2, 1] > 180) {
        bbox [2, 1] <- bbox [2, 2] - 360
    }
    if (bbox [1, 2] < -90) {
        bbox [1, 2] <- -90
    }
    if (bbox [2, 2] > 90) {
        bbox [2, 2] <- 90
    }

    list (bbox = bbox, bbox_poly = bbox_poly)
}

vec_to_bbox <- function (bbox) {

    if (is.null (names (bbox))) {

        # presume (lon, lat, lon, lat)
        bbox <- cbind (
            sort (bbox [c (1, 3)]),
            sort (bbox [c (2, 4)])
        )
        colnames (bbox) <- c ("x", "y")
        rownames (bbox) <- c ("min", "max")
    } else {

        minptn <- c ("^min", "min$", "0$")
        chk <- vapply (
            minptn, function (i) {
                any (grepl (i, names (bbox)))
            },
            logical (1)
        )
        if (length (which (chk)) != 1L) {
            stop (
                "names of bbox elements should clearly label ",
                "min & max longitude and latitude"
            )
        }

        minptn <- minptn [which (chk)]
        mincols <- grep (minptn, names (bbox))

        if (minptn == "min$") {

            maxcols <- grep ("max$", names (bbox))
            xcols <- grep ("^x|^lon", names (bbox))
            ycols <- grep ("^y|^lat", names (bbox))

        } else if (minptn == "^min") {

            maxcols <- grep ("^max", names (bbox))
            xcols <- grep ("x$|lon$", names (bbox))
            ycols <- grep ("y$|lat$", names (bbox))

        } else if (minptn == "0$") {

            maxcols <- grep ("1$", names (bbox))
            xcols <- grep ("^x|^lon", names (bbox))
            ycols <- grep ("^y|^lat", names (bbox))
        }

        if (!(length (mincols) == 2 &&
            length (maxcols) == 2 &&
            length (xcols) == 2 &&
            length (ycols) == 2)) {

            stop (
                "names of bbox elements should clearly label ",
                "min & max longitude and latitude"
            )
        }

        lonmin <- xcols [which (xcols %in% mincols)]
        latmin <- ycols [which (ycols %in% mincols)]
        lonmax <- xcols [which (xcols %in% maxcols)]
        latmax <- ycols [which (ycols %in% maxcols)]

        bbox <- rbind (
            c (bbox [lonmin], bbox [lonmax]),
            c (bbox [latmin], bbox [latmax])
        )
    }

    return (bbox)
}
