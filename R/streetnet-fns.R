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

    # osm_poly2line merges all street polygons with the line ones
    net <- osmdata::opq (bbox) %>%
                osmdata::add_osm_feature (key = "highway") %>%
                osmdata::osmdata_sf (quiet = quiet) %>%
                osmdata::osm_poly2line () %>%
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
#' objects, typically extracted with \code{dodgr_streetnet}
#' @param wt_profile Name of weighting profile, or vector of values with names
#' corresponding to names in \code{type_col}
#' @param type_col Specify column of the \code{sf} \code{data.frame} object
#' which designates different types of highways to be used for weighting
#' (default works with \code{osmdata} objects).
#' @param id_col Specify column of the code{sf} \code{data.frame} object which
#' provides unique identifiers for each highway (default works with
#' \code{osmdata} objects).
#' @param keep_cols Vectors of columns from \code{sf_lines} to be kept in the
#' resultant \code{dodgr} network; vector can be either names or indices of
#' desired columns.
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
                              type_col = "highway", id_col = "osm_id",
                              keep_cols = NULL)
{
    if (!is (sf_lines, "sf"))
        stop ('sf_lines must be class "sf"')
    #if (!all (c ("geometry", type_col, id_col) %in% names (sf_lines)))
    #    stop (paste0 ('sf_lines must be class "sf" and ',
    #                  'have highway and geometry columns'))
    if (!"geometry" %in% names (sf_lines))
        stop (paste0 ('sf_lines must be class "sf" and have a geometry column'))

    if (type_col != "highway")
        names (sf_lines) [which (names (sf_lines) == type_col)] <- "highway"
    if (id_col != "osm_id")
        names (sf_lines) [which (names (sf_lines) == id_col)] <- "osm_id"

    if (!"highway" %in% names (sf_lines))
        stop ("Please specify type_col to be used for weighting streetnet")
    if (!"osm_id" %in% names (sf_lines))
        stop ("Please specifiy id_col to be used to identify streetnet rows")

    if (is.null (names (sf_lines$geometry)))
        names (sf_lines$geometry) <- sf_lines$osm_id

    if (is.character (wt_profile))
    {
        prf_names <- c ("foot", "horse", "wheelchair", "bicycle", "moped",
                        "motorcycle", "motorcar", "goods", "hgv", "psv")
        wt_profile <- match.arg (tolower (wt_profile), prf_names)
        profiles <- dodgr::weighting_profiles
        wt_profile <- profiles [profiles$name == wt_profile, ]
    } else if (is.numeric (wt_profile))
    {
        nms <- names (wt_profile)
        if (is.null (nms))
            nms <- NA
        wt_profile <- data.frame (name = "custom",
                                  way = nms,
                                  value = wt_profile,
                                  stringsAsFactors = FALSE)
    } else if (is.data.frame (wt_profile))
    {
        # assert that is has the standard structure
        if (ncol (wt_profile) != 3 |
            !identical (names (wt_profile), c ("name", "way", "value")))
            stop ("Weighting profiles must have three columsn of ",
                  "(name, way, value); see 'weighting_profiles' for examples")
    } else
        stop ("Custom named profiles must be vectors with named values")

    if (nrow (wt_profile) > 1 & all (wt_profile$name != "custom"))
        sf_lines <- remap_way_types (sf_lines, wt_profile)

    dat <- rcpp_sf_as_network (sf_lines, pr = wt_profile)
    graph <- data.frame (geom_num = dat$numeric_values [, 1] + 1, # 1-indexed!
                         edge_id = seq (nrow (dat$character_values)),
                         from_id = as.character (dat$character_values [, 1]),
                         from_lon = dat$numeric_values [, 2],
                         from_lat = dat$numeric_values [, 3],
                         to_id = as.character (dat$character_values [, 2]),
                         to_lon = dat$numeric_values [, 4],
                         to_lat = dat$numeric_values [, 5],
                         d = dat$numeric_values [, 6],
                         d_weighted = dat$numeric_values [, 7],
                         highway = as.character (dat$character_values [, 3]),
                         way_id = as.character (dat$character_values [, 4]),
                         stringsAsFactors = FALSE
                         )
    # rcpp_sf_as_network now flags non-routable ways with -1, so:
    graph$d_weighted [graph$d_weighted < 0] <- max (graph$d_weighted) * 1e6
    if (all (graph$highway == ""))
        graph$highway <- NULL
    if (all (graph$way_id == ""))
        graph$way_id <- NULL

    # If original geometries did not have rownames (meaning it's not from
    # osmdata), then reassign unique vertex from/to IDs based on coordinates
    if (is.null (rownames (as.matrix (sf_lines$geometry [[1]]))))
    {
        # Find unique vertices by coorindates along:
        xy <- data.frame (x = c (graph$from_lon, graph$to_lon),
                          y = c (graph$from_lat, graph$to_lat))
        indx <- which (!duplicated (xy))
        #ids <- seq (indx) - 1 # 0-indexed
        xy_indx <- xy [indx, ]
        xy_indx$indx <- seq (nrow (xy_indx))

        # Then match coordinates to those unique values and replace IDs. Note
        # that `merge` does not preserve row order, see:
        # https://www.r-statistics.com/2012/01/merging-two-data-frame-objects-while-preserving-the-rows-order/
        xy_from <- data.frame (x = graph$from_lon, y = graph$from_lat,
                               ord = seq (nrow (graph)))
        xy_from <- merge (xy_from, xy_indx) # re-sorts rows
        graph$from_id <- as.character (xy_from$indx [order (xy_from$ord)])
        xy_to <- data.frame (x = graph$to_lon, y = graph$to_lat,
                             ord = seq (nrow (graph)))
        xy_to <- merge (xy_to, xy_indx)
        graph$to_id <- as.character (xy_to$indx [order (xy_to$ord)])
    }

    # get component numbers for each edge
    class (graph) <- c (class (graph), "dodgr_streetnet")
    graph <- dodgr_components (graph)

    # And finally, re-insert keep_cols:
    if (length (keep_cols) > 0)
    {
        keep_names <- NULL
        if (is.character (keep_cols))
        {
            keep_names <- keep_cols
            keep_cols <- match (keep_cols, names (sf_lines))
        } else if (is.numeric (keep_cols))
        {
            keep_names <- names (sf_lines) [keep_cols]
        } else
        {
            stop ("keep_cols must be either character or numeric")
        }
        indx <- match (graph$geom_num, seq (sf_lines$geometry))
        for (k in seq (keep_names))
            graph [[keep_names [k] ]] <- sf_lines [indx, keep_cols [k],
                                                   drop = TRUE]
    }

    return (graph)
}

# re-map any OSM 'highway' types with pmatch to standard types
remap_way_types <- function (sf_lines, wt_profile)
{
    way_types <- unique (as.character (sf_lines$highway))
    dodgr_types <- unique (wt_profile$way)
    # clearer to code as a for loop
    for (i in seq (way_types))
    {
        if (!way_types [i] %in% dodgr_types)
        {
            pos <- which (pmatch (dodgr_types, way_types  [i]) > 0)
            if (length (pos) > 0)
                sf_lines$highway [sf_lines$highway == way_types [i]] <-
                    dodgr_types [pos]
        }
    }

    # re-map some common types
    indx <- which (sf_lines$highway %in% c ("pedestrian", "footway",
                                            "bridleway"))
    if (length (indx) > 0)
        sf_lines$highway [indx] <- "path"

    way_types <- unique (as.character (sf_lines$highway))
    not_in_wt_prof <- way_types [which (!way_types %in% dodgr_types)]
    if (length (not_in_wt_prof) > 0)
        message ("The following highway types are present in data yet ",
                 "lack corresponding weight_profile values: ",
                 paste0 (not_in_wt_prof, sep = ", "))

    # remove not_in_wt_prof types. Note that this subsetting strips most sf
    # attributes, so is not sf-compliant, but that's okay here because this
    # object is only passed to internal C++ routines,
    indx <- which (!sf_lines$highway %in% not_in_wt_prof)
    sf_lines <- sf_lines [indx, ]

    return (sf_lines)
}

#' dodgr_to_sfc
#'
#' Convert a \code{dodgr} graph as \code{data.frame} of all edges into a Simple
#' Features (\pkg{sf}) object through aggregating edges into \code{LINESTRING}
#' objects representing longest sequences between all junction nodes. The
#' resultant objects will generally contain more \code{LINESTRING} objects than
#' the original \pkg{sf} object, because the former will be bisected at every
#' junction point.
#'
#' @param net A \pkg{dodgr} network
#' @return A list containing (1) A \code{data.frame} of data associated with the
#' `sf` geometries; and (ii) A Simple Features Collection (\code{sfc}) list of
#' \code{LINESTRING} objects.
#'
#' @note The output of this function corresponds to the edges obtained from
#' \code{dodgr_contract_graph}. An \pkg{sf} \code{data.frame} may be created by
#' appending any data from the latter to the \code{sfc} output of this function -
#' see \pkg{sf} for details.
#'
#' @export
#' @examples
#' hw <- weight_streetnet (hampi)
#' xy <- dodgr_to_sfc (hw)
#' dim (hw) # 5.845 edges
#' length (xy$geom) # 682 aggregated linestrings aggregated from those edges
#' nrow (hampi) # compared to 191 linestrings in original sf object
#' dim (xy$dat) # same number of rows as there are geometries
dodgr_to_sfc <- function (net)
{
    # force sequential IDs. TODO: Allow non-sequential by replacing indices in
    # src/dodgr_to_sf::get_edge_to_vert_maps with maps to sequential indices.
    net$edge_id <- seq (nrow (net))

    gc <- dodgr_contract_graph (net)
    geoms <- rcpp_aggregate_to_sf (net, gc$graph, gc$edge_map)

    # Then match data of `net` potentially including way_id, back on to the
    # geometries:
    #edge_ids <- gc$graph$edge_id [match (names (geoms), gc$graph$edge_id)]
    #indx1 <- which (edge_ids %in% gc$edge_map$edge_new)
    #indx2 <- seq (edge_ids) [!seq (edge_ids) %in% indx1]
    #edge_ids <- c (gc$edge_map$edge_old [indx1], edge_ids [indx2])
    #index <- match (edge_ids, net$edge_id)
    #dat <- net [index, ]
    #dat$from_id <- dat$from_lat <- dat$from_lon <- NULL
    #dat$to_id <- dat$to_lat <- dat$to_lon <- NULL
    #dat$d <- dat$d_weighted <- dat$edge_id <- NULL

    geoms <- geoms [match (gc$graph$edge_id, names (geoms))]

    return (list (dat = gc$graph, geoms = geoms))
}

#' dodgr_to_sf
#' @inherit dodgr_to_sfc
#' @export
dodgr_to_sf <- function (net)
{
    .Deprecated ("dodgr_to_sfc")
    dodgr_to_sfc (net)
}
