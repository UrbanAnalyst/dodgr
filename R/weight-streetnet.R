#' weight_streetnet
#'
#' Weight (or re-weight) an `sf`-formatted OSM street network according to
#' a named routino profile, selected from (foot, horse, wheelchair, bicycle,
#' moped, motorcycle, motorcar, goods, hgv, psv).
#'
#' @param x A street network represented as `sf` `LINESTRING`
#' objects, typically extracted with `dodgr_streetnet`
#' @param wt_profile Name of weighting profile, or vector of values with names
#' corresponding to names in `type_col` (see Details)
#' @param type_col Specify column of the `sf` `data.frame` object
#' which designates different types of highways to be used for weighting
#' (default works with `osmdata` objects).
#' @param id_col Specify column of the code{sf} `data.frame` object which
#' provides unique identifiers for each highway (default works with
#' `osmdata` objects).
#' @param keep_cols Vectors of columns from `sf_lines` to be kept in the
#' resultant `dodgr` network; vector can be either names or indices of
#' desired columns.
#'
#' @return A `data.frame` of edges representing the street network, with
#' distances in kilometres, along with a column of graph component numbers.
#'
#' @note Names for the `wt_profile` parameter are taken from
#' \link{weighting_profiles}, which is a `data.frame` of weights for
#' different modes of transport. Current modes included there are "bicycle",
#' "foot", "goods", "hgv", "horse", "moped", "motorcar", "motorcycle", "psv",
#' and "wheelchair". Railway routing can be implemented with
#' the separate function \link{weight_railway}.
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
#'
#' # An example for a generic (non-OSM) highway, represented as the
#' # `routes_fast` object of the \pkg{stplanr} package, which is a
#' # SpatialLinesDataFrame.
#' \dontrun{
#' library (stplanr)
#' # merge all of the 'routes_fast' lines into a single network
#' r <- overline (routes_fast, attrib = "length", buff_dist = 1)
#' r <- sf::st_as_sf (r)
#' # We need to specify both a `type` and `id` column for the
#' # \link{weight_streetnet} function.
#' r$type <- 1
#' r$id <- seq (nrow (r))
#' graph <- weight_streetnet (r, type_col = "type", id_col = "id",
#'                            wt_profile = 1)
#' }
weight_streetnet <- function (x, wt_profile = "bicycle",
                              type_col = "highway", id_col = "osm_id",
                              keep_cols = NULL)
{
    UseMethod ("weight_streetnet")
}

#' @name weight_streetnet
#' @export
weight_streetnet.default <- function (x, wt_profile = "bicycle",
                              type_col = "highway", id_col = "osm_id",
                              keep_cols = NULL)
{
    stop ("Unknown class")
}

# ********************************************************************
# *************************     sf class     ************************* 
# ********************************************************************

#' @name weight_streetnet
#' @export
weight_streetnet.sf <- function (x, wt_profile = "bicycle",
                              type_col = "highway", id_col = "osm_id",
                              keep_cols = NULL)
{
    geom_column <- get_sf_geom_col (x)
    attr (x, "sf_column") <- geom_column

    if (type_col != "highway")
        names (x) [which (names (x) == type_col)] <- "highway"
    if (id_col != "osm_id")
        names (x) [which (names (x) == id_col)] <- "osm_id"

    if (!"highway" %in% names (x) & !is.numeric (wt_profile))
        stop ("Please specify type_col to be used for weighting streetnet")
    if (!"osm_id" %in% names (x))
    {
        idcol <- grep ("id", names (x), ignore.case = TRUE)
        if (length (idcol) == 1)
        {
            message ("Using column ", names (x) [idcol],
                     " as ID column for edges; please specify explicitly if",
                     " a different column should be used.")
            names (x) [idcol] <- "osm_id"
        } else if (length (idcol) > 1)
        {
            stop ("Multiple potential ID columns: [",
                  paste0 (names (x) [idcol], collapse = " "),
                  "]; please explicitly specify one of these.")
        } else if (length (idcol) == 1)
        {
            message ("x appears to have no ID column;",
                     "sequential edge numbers will be used.")
            x$osm_id <- seq (nrow (x))
        }
    }

    if (is.null (names (x [geom_column])))
        names (x [geom_column]) <- x$osm_id
    # Then rename geom_column to "geometry" for the C++ routine
    names (x) [match (geom_column, names (x))] <- "geometry"
    attr (x, "sf_column") <- "geometry"

    if (is.character (wt_profile))
    {
        if (grepl ("rail", wt_profile, ignore.case = TRUE))
        {
            stop ("Please use the weight_railway function for railway routing.")
        } else
        {
            prf_names <- c ("foot", "horse", "wheelchair", "bicycle", "moped",
                            "motorcycle", "motorcar", "goods", "hgv", "psv")
            wt_profile <- match.arg (tolower (wt_profile), prf_names)
            profiles <- dodgr::weighting_profiles
            wt_profile <- profiles [profiles$name == wt_profile, ]
        }
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
        x <- remap_way_types (x, wt_profile)

    dat <- rcpp_sf_as_network (x, pr = wt_profile)
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
    graph$d_weighted [graph$d_weighted < 0] <- .Machine$double.xmax
    if (all (graph$highway == ""))
        graph$highway <- NULL
    if (all (graph$way_id == ""))
        graph$way_id <- NULL

    # If original geometries did not have rownames (meaning it's not from
    # osmdata), then reassign unique vertex from/to IDs based on coordinates
    if (is.null (rownames (as.matrix (x$geometry [[1]]))))
        graph <- rownames_from_xy (graph)

    # get component numbers for each edge
    class (graph) <- c (class (graph), "dodgr_streetnet")
    graph <- dodgr_components (graph)

    # And finally, re-insert keep_cols:
    if (length (keep_cols) > 0)
        graph <- reinsert_keep_cols (x, graph, keep_cols)

    graph$d_weighted [graph$d_weighted == .Machine$double.xmax] <- NA

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
    indx <- which (sf_lines$highway %in% c ("pedestrian", "footway"))
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

# Return the name of the sf geometry column, which this routines permits to be
# either anything that greps "geom" (so "geom", "geoms", "geometry"), or else
# just plain "g". See Issue#66.
get_sf_geom_col <- function (graph)
{
    gcol <- grep ("geom", names (graph))
    if (length (gcol) > 1)
    {
        gnames <- c ("geometry", "geom", "geoms")
        mg <- match (gnames, names (graph))
        if (length (which (!is.na (mg))) == 1)
            gcol <- mg [which (!is.na (mg))]
        else
            stop ("Unable to determine geometry column from [",
                  paste0 (names (graph) [gcol], collapse = ", "), "]")
    } else if (length (gcol) == 0)
    {
        gcol <- match ("g", names (graph))
        if (is.na (gcol) | length (gcol) != 1)
            stop ("Unable to determine geometry column")
    }

    return (names (graph) [gcol])
}

rownames_from_xy <- function (graph)
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
    xy_from <- data.frame (x = graph$from_lon,
                           y = graph$from_lat,
                           ord = seq (nrow (graph)))
    xy_from <- merge (xy_from, xy_indx) # re-sorts rows
    graph$from_id <- as.character (xy_from$indx [order (xy_from$ord)])
    xy_to <- data.frame (x = graph$to_lon,
                         y = graph$to_lat,
                         ord = seq (nrow (graph)))
    xy_to <- merge (xy_to, xy_indx)
    graph$to_id <- as.character (xy_to$indx [order (xy_to$ord)])

    return (graph)
}

reinsert_keep_cols <- function (sf_lines, graph, keep_cols)
{
    keep_names <- NULL
    if (is.character (keep_cols))
    {
        keep_names <- keep_cols
        keep_cols <- match (keep_cols, names (sf_lines))
        # NA is no keep_cols match 
    } else if (is.numeric (keep_cols))
    {
        keep_names <- names (sf_lines) [keep_cols]
    } else
    {
        stop ("keep_cols must be either character or numeric")
    }
    indx <- which (is.na (keep_cols))
    if (length (indx) > 0)
        message ("Data has no columns named ",
                 paste0 (keep_names, collapse = ", "))
    keep_cols <- keep_cols [!is.na (keep_cols)]
    if (length (keep_cols) > 0)
    {
        indx <- match (graph$geom_num, seq (sf_lines$geometry))
        if (!is.na (keep_cols))
            for (k in seq (keep_names))
                graph [[keep_names [k] ]] <- sf_lines [indx, keep_cols [k],
                                                       drop = TRUE]
    }

    return (graph)
}

# ********************************************************************
# *************************     sc class     ************************* 
# ********************************************************************

#' @name weight_streetnet
#' @export
weight_streetnet.sc <- weight_streetnet.SC <- function (x, wt_profile = "bicycle",
                                                        type_col = "highway",
                                                        id_col = "osm_id",
                                                        keep_cols = NULL)
{
    requireNamespace ("geodist")
    requireNamespace ("dplyr")
    check_sc (x)

    extract_sc_edges_xy (x) %>%
        sc_edge_dist () %>%
        extract_sc_edges_highways (x) %>%
        weight_sc_edges (wt_profile) %>%
        sc_lanes_surface (wt_profile) %>%
        sc_edge_time (wt_profile, x)
}

has_elevation <- function (x)
{
    "z_" %in% names (x$vertex)
}

check_sc <- function (x)
{
    if (!"osmdata_sc" %in% class (x))
        stop ("weight_streetnet currently only works for 'sc'-class objects ",
              "extracted with osmdata::osmdata_sc.")
}

# First step of edge extraction: join x and y coordinates
extract_sc_edges_xy <- function (x)
{
    rename0 <- c (.vx0_x = "x_", .vx0_y = "y_", .vx0_z = "z_")
    rename1 <- c (.vx1_x = "x_", .vx1_y = "y_", .vx1_z = "z_")
    if (!has_elevation (x))
    {
        rename0 <- rename0 [1:2]
        rename1 <- rename1 [1:2]
    }

    dplyr::left_join (x$edge, x$vertex, by = c (".vx0" = "vertex_")) %>%
        dplyr::rename (!!rename0) %>%
        dplyr::left_join (x$vertex, by = c (".vx1" = "vertex_")) %>%
        dplyr::rename (!!rename1)
}

sc_edge_dist <- function (edges)
{
    # no visible binding notes:
    .vx0_z <- .vx1_z <- NULL

    edges$d <- geodist::geodist (edges [, c (".vx0_x", ".vx0_y")],
                                 edges [, c (".vx1_x", ".vx1_y")], paired = TRUE)
    if (".vx0_z" %in% names (edges) & ".vx1_z" %in% names (edges))
        edges <- dplyr::mutate (edges, "dz" = .vx1_z - .vx0_z) %>%
            dplyr::select (-c(.vx0_z, .vx1_z))
    return (edges)
}

extract_sc_edges_highways <- function (edges, x)
{
    # no visible binding notes:
    native_ <- key <- `:=` <- value <- NULL

    edges <- dplyr::left_join (edges, x$object_link_edge, by = "edge_") %>%
        dplyr::select (-native_)
    keep_types <- c ("highway", "oneway", "oneway:bicycle", "lanes",
                     "maxspeed", "surface")
    for (k in keep_types)
    {
        objs <- dplyr::filter (x$object, key == k)
        edges <- dplyr::left_join (edges, objs, by = "object_") %>%
            dplyr::rename (!!dplyr::quo_name (k) := value) %>%
            dplyr::select (-key)
    }
    # oneway:bicycle doesn't enquote properly, so:
    i <- grep ("bicycle", names (edges))
    names (edges) [i] <- "oneway_bicycle"

    # re-map the oneway values to boolean
    edges$oneway [!edges$oneway %in% c ("no", "yes")] <- "no"
    edges$oneway <- ifelse (edges$oneway == "no", FALSE, TRUE)
    edges$oneway_bicycle [!edges$oneway_bicycle %in% c ("no", "yes")] <- "no"
    edges$oneway_bicycle <- ifelse (edges$oneway_bicycle == "no", FALSE, TRUE)

    # TODO: Do this in convert_graph function
    # replace NA distances and times
    #m <- .Machine$double.xmax
    #edges$d_weighted [is.na (edges$d_weighted)] <- m
    #edges$time [is.na (edges$time)] <- m

    return (edges)
}

weight_sc_edges <- function (edges, wt_profile)
{
    # no visible binding notes:
    value <- d <- NULL

    wp <- dodgr::weighting_profiles
    wp <- wp [wp$name == wt_profile, c ("way", "value")]
    dplyr::left_join (edges, wp, by = c ("highway" = "way")) %>%
        dplyr::filter (!is.na (value)) %>%
        dplyr::mutate (d_weighted = ifelse (value == 0, NA, d / value)) %>%
        dplyr::select (-value)
}

# adjust weighted distances according to numbers of lanes and surfaces
# NOTE: This is currently hard-coded for active transport only, and will not
# work for anything else
sc_lanes_surface <- function (edges, wt_profile)
{
    # no visible binding notes:
    lanes <- maxspeed <- surface <- NULL

    if (wt_profile %in% c ("foot", "bicycle"))
    {
        lns <- c (4, 5, 6, 7, 8)
        wts <- c (0.05, 0.05, 0.1, 0.1, 0.2)
        for (i in seq (lns))
        {
            index <- which (edges$lanes == lns [i])
            if (i == length (lns))
                index <- which (edges$lanes >= lns [i])
            edges$d_weighted [index] <- edges$d_weighted [index] * (1 + wts [i])
        }

        srf <- c ("dirt", "grass", "natural")
        wts <- c (0.1, 0.1, 0.1)
        for (i in seq (srf))
        {
            index <- which (edges$surface == srf [i])
            edges$d_weighted [index] <- edges$d_weighted [index] * (1 + wts [i])
        }
    }
    edges <- dplyr::select (edges, -c(lanes, maxspeed, surface))

    return (edges)
}

# Convert weighted distances in metres to time in seconds
# NOTE: This is currently hard-coded for foot and bicycle only, and will not work for
# anything else. 
sc_edge_time <- function (edges, wt_profile, x)
{
    # no visible binding messages:
    object_ <- NULL

    if (wt_profile == "foot")
    {
        # Uses 
        # [Naismith's Rule](https://en.wikipedia.org/wiki/Naismith%27s_rule)
        edges$time = 3600 * (edges$d_weighted / 1000) / 5 # 5km per hour
        if ("dz" %in% names (edges))
        {
            index <- which (edges$dz > 0)
            edges$time [index] <- edges$time [index] + edges$dz [index] / 10
            edges$dz <- NULL
        }

        # add waiting times at traffic lights
        wp <- dodgr::weighting_profiles
        wait <- wp$value [wp$name == paste0 ("lights_", wt_profile)]

        crossings <- traffic_light_objs_ped (x) # way IDs
        objs <- x$object %>% dplyr::filter (object_ %in% crossings$crossings)
        oles <- x$object_link_edge %>% dplyr::filter (object_ %in% objs$object_)
        index <- match (oles$edge_, edges$edge_)
        edges$time [index] <- edges$time [index] + wait
    } else if (wt_profile == "bicycle")
    {
        # http://theclimbingcyclist.com/gradients-and-cycling-how-much-harder-are-steeper-climbs/
        # http://cycleseven.org/effect-of-hills-on-cycling-effort
        # The latter argues for a linear relationship with a reduction in speed
        # of "about 11% for every 1% change in steepness". For 0.01 to translate
        # to 0.11, it needs to be multiplied by 0.11 / 0.01, or 11
        edges$time = 3600 * (edges$d_weighted / 1000) / 12 # 12km per hour
        if ("dz" %in% names (edges))
        {
            index <- which (edges$dz > 0)
            edges$time [index] <- edges$time [index] * 
                (1 + 11 * edges$dz [index] / edges$d [index])
            edges$dz <- NULL
        }
        # ... TODO: Downhill
        # http://www.sportsci.org/jour/9804/dps.html
        # downhill cycling speed ~ sqrt (slope)
    }
    return (edges)
}

# ********************************************************************
# **********************     weight railway     **********************
# ********************************************************************

#' weight_railway
#'
#' Weight (or re-weight) an `sf`-formatted OSM street network for routing
#' along railways.
#'
#' @param sf_lines A street network represented as `sf` `LINESTRING`
#' objects, typically extracted with `dodgr_streetnet`
#' @param type_col Specify column of the `sf` `data.frame` object
#' which designates different types of railways to be used for weighting
#' (default works with `osmdata` objects).
#' @param id_col Specify column of the code{sf} `data.frame` object which
#' provides unique identifiers for each railway (default works with
#' `osmdata` objects).
#' @param keep_cols Vectors of columns from `sf_lines` to be kept in the
#' resultant `dodgr` network; vector can be either names or indices of
#' desired columns.
#' @param excluded Types of railways to exclude from routing.
#'
#' @return A `data.frame` of edges representing the rail network, along
#' with a column of graph component numbers.
#'
#' @note Default railway weighting is by distance. Other weighting schemes, such
#' as by maximum speed, can be implemented simply by modifying the
#' `d_weighted` column returned by this function accordingly.
#'
#' @examples
#' \dontrun{
#' # sample railway extraction with the 'osmdata' package
#' library (osmdata)
#' dat <- opq ("shinjuku") %>%
#'     add_osm_feature (key = "railway") %>%
#'     osmdata_sf (quiet = FALSE)
#' graph <- weight_railway (dat$osm_lines)
#' }
#' @export
weight_railway <- function (sf_lines, type_col = "railway", id_col = "osm_id",
                            keep_cols = c ("maxspeed"),
                            excluded = c ("abandoned", "disused",
                                          "proposed", "razed"))
{
    if (!is (sf_lines, "sf"))
        stop ('sf_lines must be class "sf"')
    geom_column <- get_sf_geom_col (sf_lines)
    attr (sf_lines, "sf_column") <- geom_column

    if (type_col != "railway")
        names (sf_lines) [which (names (sf_lines) == type_col)] <- "railway"
    if (id_col != "osm_id")
        names (sf_lines) [which (names (sf_lines) == id_col)] <- "osm_id"

    if (!"railway" %in% names (sf_lines))
        stop ("Please specify type_col to be used for weighting railway")
    if (!"osm_id" %in% names (sf_lines))
        stop ("Please specifiy id_col to be used to identify railway rows")

    if (is.null (names (sf_lines$geometry)))
        names (sf_lines$geometry) <- sf_lines$osm_id


    sf_lines <- sf_lines [which (!(sf_lines$railway %in% excluded |
                                   is.na (sf_lines$railway))), ]
    # routing is based on matching the given profile to the "highway" field of
    # sf_lines, so:
    sf_lines$highway <- sf_lines$railway

    wt_profile <- data.frame (name = "custom",
                              way = unique (sf_lines$highway),
                              value = 1,
                              stringsAsFactors = FALSE)

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
    graph$d_weighted [graph$d_weighted < 0] <- .Machine$double.xmax
    if (all (graph$highway == ""))
        graph$highway <- NULL
    if (all (graph$way_id == ""))
        graph$way_id <- NULL

    # If original geometries did not have rownames (meaning it's not from
    # osmdata), then reassign unique vertex from/to IDs based on coordinates
    if (is.null (rownames (as.matrix (sf_lines$geometry [[1]]))))
        graph <- rownames_from_xy (graph)

    # get component numbers for each edge
    class (graph) <- c (class (graph), "dodgr_streetnet")
    graph <- dodgr_components (graph)

    # And finally, re-insert keep_cols:
    if (length (keep_cols) > 0)
        graph <- reinsert_keep_cols (sf_lines, graph, keep_cols)

    return (graph)
}
