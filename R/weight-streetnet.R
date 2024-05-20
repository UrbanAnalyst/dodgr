#' Weight a street network according to a specified weighting profile.
#'
#' Weight (or re-weight) an \pkg{sf} or `SC` (`silicate`)-formatted OSM street
#' network according to a named profile, selected from (foot, horse, wheelchair,
#' bicycle, moped, motorcycle, motorcar, goods, hgv, psv), or a cusstomized
#' version dervied from those.
#'
#' @param x A street network represented either as `sf` `LINESTRING`
#' objects, typically extracted with \link{dodgr_streetnet}, or as an `SC`
#' (`silicate`) object typically extracted with the \link{dodgr_streetnet_sc}.
#' @param wt_profile Name of weighting profile, or `data.frame` specifying
#' custom values (see Details)
#' @param wt_profile_file Name of locally-stored, `.json`-formatted version of
#' `dodgr::weighting_profiles`, created with \link{write_dodgr_wt_profile}, and
#' modified as desired.
#' @param type_col Specify column of the `sf` `data.frame` object
#' which designates different types of highways to be used for weighting
#' (default works with `osmdata` objects).
#' @param id_col For `sf`-formatted data only: Specify column of the \pkg{sf}
#' `data.frame` object which provides unique identifiers for each highway
#' (default works with `osmdata` objects).
#' @param keep_cols Vectors of columns from `x` to be kept in the resultant
#' `dodgr` network; vector can be either names, regex-patterns,  or indices of
#' desired columns (see notes).
#' @param turn_penalty Including time penalty on edges for turning across
#' oncoming traffic at intersections (see Note).
#' @param left_side Does traffic travel on the left side of the road (`TRUE`) or
#' the right side (`FALSE`)? - only has effect on turn angle calculations for
#' edge times.
#'
#' @return A `data.frame` of edges representing the street network, with
#' distances in metres and times in seconds, along with a column of graph
#' component numbers. Times for \pkg{sf}-formatted street networks are only
#' approximate, and do not take into account traffic lights, turn angles, or
#' elevation changes. Times for \pkg{sc}-formatted street networks take into
#' account all of these factors, with elevation changes automatically taken into
#' account for networks generated with the \pkg{osmdata} function
#' `osm_elevation()`.
#'
#' @note Names for the `wt_profile` parameter are taken from
#' \link{weighting_profiles}, which is a list including a `data.frame` also
#' called `weighting_profiles` of weights for different modes of transport.
#' Values for `wt_profile` are taken from current modes included there, which
#' are "bicycle", "foot", "goods", "hgv", "horse", "moped", "motorcar",
#' "motorcycle", "psv", and "wheelchair". Railway routing can be implemented
#' with the separate function \link{weight_railway}. Alternatively, the entire
#' `weighting_profile` structures can be written to a local `.json`-formatted
#' file with \link{write_dodgr_wt_profile}, the values edited as desired, and
#' the name of this file passed as the `wt_profile_file` parameter.
#'
#' @note Realistic routing include factors such as access restrictions, turn
#' penalties, and effects of incline, can only be implemented when the objects
#' passed to `weight_streetnet` are of \pkg{sc} ("silicate") format, generated
#' with \link{dodgr_streetnet_sc}. Restrictions applied to ways (in Open
#' Streetmap Terminology) may be controlled by ensuring specific columns are
#' retained in the `dodgr` network with the `keep_cols` argument. For example,
#' restrictions on access are generally specified by specifying a value for the
#' key of "access". Include "access" in `keep_cols` will ensure these values are
#' retained in the `dodgr` version, from which ways with specified values can
#' easily be removed or modified, as demonstrated in the examples.
#'
#' The additional Open Street Map (OSM) keys which can be used to specify
#' restrictions are which are automatically extracted with
#' \link{dodgr_streetnet_sc}, and so may be added to the `keep_cols` argument,
#' include:
#' \itemize{
#' \item "access"
#' \item "bicycle"
#' \item "foot"
#' \item "highway"
#' \item "motorcar"
#' \item "motor_vehicle"
#' \item "restriction"
#' \item "toll"
#' \item "vehicle"
#' }
#'
#' Restrictions and time-penalties on turns can be implemented from such
#' objects by setting `turn_penalty = TRUE`. Setting `turn_penalty = TRUE` will
#' honour turn restrictions specified in Open Street Map (unless the "penalties"
#' table of \link{weighting_profiles} has `restrictions = FALSE` for a specified
#' `wt_profile`). Resultant graphs are fundamentally different from the default
#' for distance-based routing. These graphs may be used directly in the
#' \link{dodgr_dists} function. Use in any other functions requires additional
#' information obtained in a file in the temporary directory of the current R
#' session with a name starting with "dodgr_junctions_", and including the
#' value of `attr(graph, "hash")`. If graphs with turn penalties are to be used
#' in subsequent R sessions, this "dodgr_junctions_" file will need to be moved
#' to a more permanent storage location, and then replaced in the temporary
#' directory of any subsequent R sessions.
#'
#' @note The resultant graph includes only those edges for which the given
#' weighting profile specifies finite edge weights. Any edges of types not
#' present in a given weighting profile are automatically removed from the
#' weighted streetnet.
#'
#' @note If the resultant graph is to be contracted via
#' \link{dodgr_contract_graph}, **and** if the columns of the graph have been,
#' or will be, modified, then automatic caching must be switched off with
#' \link{dodgr_cache_off}. If not, the \link{dodgr_contract_graph} function will
#' return the automatically cached version, which is the contracted version of
#' the full graph prior to any modification of columns.
#'
#' @seealso \link{write_dodgr_wt_profile}, \link{dodgr_times}
#'
#' @family extraction
#' @export
#' @examples
#' # hampi is included with package as an 'osmdata' sf-formatted street network
#' net <- weight_streetnet (hampi)
#' class (net) # data.frame
#' dim (net) # 6096  11; 6096 streets
#' # os_roads_bristol is also included as an sf data.frame, but in a different
#' # format requiring identification of columns and specification of custom
#' # weighting scheme.
#' colnm <- "formOfWay"
#' wts <- data.frame (
#'     name = "custom",
#'     way = unique (os_roads_bristol [[colnm]]),
#'     value = c (0.1, 0.2, 0.8, 1)
#' )
#' net <- weight_streetnet (
#'     os_roads_bristol,
#'     wt_profile = wts,
#'     type_col = colnm, id_col = "identifier"
#' )
#' dim (net) # 406 11; 406 streets
#'
#' # An example for a generic (non-OSM) highway, represented as the
#' # `routes_fast` object of the \pkg{stplanr} package, which is a
#' # SpatialLinesDataFrame.
#' \dontrun{
#' library (stplanr)
#' # merge all of the 'routes_fast' lines into a single network
#' r <- overline (routes_fast, attrib = "length", buff_dist = 1)
#' r <- sf::st_as_sf (r, crs = 4326)
#' # We need to specify both a `type` and `id` column for the
#' # \link{weight_streetnet} function.
#' r$type <- 1
#' r$id <- seq (nrow (r))
#' graph <- weight_streetnet (
#'     r,
#'     type_col = "type",
#'     id_col = "id",
#'     wt_profile = 1
#' )
#' }
weight_streetnet <- function (x,
                              wt_profile = "bicycle",
                              wt_profile_file = NULL,
                              turn_penalty = FALSE,
                              type_col = "highway",
                              id_col = "osm_id",
                              keep_cols = NULL,
                              left_side = FALSE) {

    UseMethod ("weight_streetnet")
}

#' @name weight_streetnet
#' @family extraction
#' @export
weight_streetnet.default <- function (x,
                                      wt_profile = "bicycle",
                                      wt_profile_file = NULL,
                                      turn_penalty = FALSE,
                                      type_col = "highway",
                                      id_col = "osm_id",
                                      keep_cols = NULL,
                                      left_side = FALSE) {

    stop ("Unknown class")
}

# ********************************************************************
# ***********************   generic variables   ***********************
# ********************************************************************

way_types_to_keep <- c (
    "bicycle",
    "foot",
    "highway",
    "oneway",
    "oneway.bicycle",
    "oneway:bicycle",
    "lanes",
    "maxspeed",
    "junction"
)

# ********************************************************************
# *************************     sf class     *************************
# ********************************************************************

#' @name weight_streetnet
#' @family extraction
#' @export
weight_streetnet.sf <- function (x,
                                 wt_profile = "bicycle",
                                 wt_profile_file = NULL,
                                 turn_penalty = FALSE,
                                 type_col = "highway",
                                 id_col = "osm_id",
                                 keep_cols = NULL,
                                 left_side = FALSE) {

    if (turn_penalty) {
        stop (
            "Turn-penalty calculations only currently implemented for ",
            "street network data generated with the ",
            "`osmdata::osmdata_sc()` function."
        )
    }
    geom_column <- get_sf_geom_col (x)
    attr (x, "sf_column") <- geom_column

    x <- change_col_names (x, type_col, "highway")
    x <- change_col_names (x, id_col, "osm_id")
    x <- check_highway_osmid (x, wt_profile)

    if (is.null (names (x [[geom_column]]))) {
        names (x [[geom_column]]) <- x$osm_id
    }
    # Then rename geom_column to "geometry" for the C++ routine
    names (x) [match (geom_column, names (x))] <- "geometry"
    attr (x, "sf_column") <- "geometry"

    wp <- get_wt_profile (x, wt_profile, wt_profile_file)
    # convert oneway and oneway*bicycle values to boolean (fn is in
    # wt_streetnet-times.R):
    x <- convert_hw_types_to_bool (x, wt_profile)

    wt_profile <- wp$wt_profile
    wt_profile_name <- wp$wt_profile_name
    if (nrow (wt_profile) > 1 && all (wt_profile$name != "custom")) {
        x <- remap_way_types (x, wt_profile)
    }

    dat <- rcpp_sf_as_network (x, pr = wt_profile)
    graph <- data.frame (
        geom_num = dat$numeric_values [, 1] + 1, # 1-indexed!
        edge_id = seq_len (nrow (dat$character_values)),
        from_id = as.character (dat$character_values [, 1]),
        from_lon = dat$numeric_values [, 2],
        from_lat = dat$numeric_values [, 3],
        to_id = as.character (dat$character_values [, 2]),
        to_lon = dat$numeric_values [, 4],
        to_lat = dat$numeric_values [, 5],
        stringsAsFactors = FALSE
    )

    graph$d <- geodist::geodist (graph [, c ("from_lon", "from_lat")],
        graph [, c ("to_lon", "to_lat")],
        paired = TRUE,
        measure = "geodesic"
    )
    graph$d_weighted <- graph$d * dat$numeric_values [, 6]

    graph$highway <- as.character (dat$character_values [, 3])
    graph$way_id <- as.character (dat$character_values [, 4])

    # rcpp_sf_as_network flags non-routable ways with -1, so:
    graph$d_weighted [graph$d_weighted < 0] <- NA
    if (all (graph$highway == "")) {
        graph$highway <- NULL
    }
    if (all (graph$way_id == "")) {
        graph$way_id <- NULL
    } # nocov

    # If original geometries did not have rownames (meaning it's not from
    # osmdata), then reassign unique vertex from/to IDs based on coordinates
    if (is.null (rownames (as.matrix (x$geometry [[1]])))) {
        graph <- rownames_from_xy (graph)
    }

    graph <- dodgr_components (graph)

    if (!is.null (wt_profile_name)) {
        if (wt_profile_name == "bicycle") {
            if (is.integer (keep_cols)) {
                keep_cols <- names (x) [keep_cols]
            }
            keep_cols <- unique (c (keep_cols, c (
                "bicycle",
                "cycleway",
                "cycleway:left",
                "cycleway:right"
            )))

        }
    }
    if (length (keep_cols) > 0) {
        graph <- reinsert_keep_cols (x, graph, keep_cols)
    }

    graph <- add_extra_sf_columns (graph, x)
    if (!is.null (wt_profile_name)) {
        graph <- set_maxspeed (graph, wt_profile_name, wt_profile_file) %>%
            weight_by_num_lanes (wt_profile_name) %>%
            calc_edge_time (wt_profile_name)
    }

    gr_cols <- dodgr_graph_cols (graph)
    graph <- graph [which (!is.na (graph [[gr_cols$d_weighted]])), ]

    class (graph) <- c (class (graph), "dodgr_streetnet")
    attr (graph, "turn_penalty") <- FALSE

    hash <- get_hash (graph, contracted = FALSE, force = TRUE)
    attr (graph, "hash") <- hash
    if (is_dodgr_cache_on ()) {
        attr (graph, "px") <- cache_graph (graph, gr_cols$edge_id)
    }

    return (graph)
}

# changed type_col and id_col to expected values of "highway" and "osm_id"
change_col_names <- function (x, colvar, expected) {

    if (colvar != expected) {
        names (x) [which (names (x) == colvar)] <- expected
    }
    return (x)
}

check_highway_osmid <- function (x, wt_profile) {

    if (!"highway" %in% names (x) && !is.numeric (wt_profile)) {
        stop (
            "Please specify type_col to be used for ", # nocov
            "weighting streetnet"
        )
    } # nocov
    if (!"osm_id" %in% names (x)) {
        idcol <- grep ("^id|id$", names (x), ignore.case = TRUE)
        if (length (idcol) == 1) {
            message (
                "Using column ", names (x) [idcol],
                " as ID column for edges; please specify explicitly if",
                " a different column should be used."
            )
            names (x) [idcol] <- "osm_id"
        } else if (length (idcol) > 1) {
            stop (
                "Multiple potential ID columns: [",
                paste0 (names (x) [idcol], collapse = " "),
                "]; please explicitly specify one of these."
            )
        } else if (length (idcol) == 0) {
            message (
                "x appears to have no ID column; ",
                "sequential edge numbers will be used."
            )
            x$osm_id <- seq_len (nrow (x))
        }
    }

    return (x)
}

get_wt_profile <- function (x, wt_profile, wt_profile_file) {

    if (!is.data.frame (wt_profile) && length (wt_profile) > 1) {
        stop ("wt_profile can only be one element")
    }

    wt_profile_name <- NULL
    if (is.character (wt_profile)) {
        if (grepl ("rail", wt_profile, ignore.case = TRUE)) {
            stop ("Please use the weight_railway function for railway routing.")
        } else {
            wt_profile_name <- wt_profile
            wt_profile <- get_profile (wt_profile_name, wt_profile_file)
        }
    } else if (is.numeric (wt_profile)) {
        nms <- names (wt_profile)
        if (is.null (nms)) {
            nms <- NA
        }
        wt_profile <- data.frame (
            name = "custom",
            way = nms,
            value = wt_profile,
            stringsAsFactors = FALSE
        )
    } else if (is.data.frame (wt_profile)) {
        # assert that is has the standard structure
        if (!all (c ("name", "way", "value") %in% names (wt_profile))) {
            stop (
                "Weighting profiles must have three columsn of ",
                "(name, way, value); see 'weighting_profiles' for examples"
            )
        }
    } else {
        stop ("Custom named profiles must be vectors with named values")
    }

    list (wt_profile = wt_profile, wt_profile_name = wt_profile_name)
}

# re-map any OSM 'highway' types with pmatch to standard types
remap_way_types <- function (sf_lines, wt_profile) {

    way_types <- unique (as.character (sf_lines$highway))
    dodgr_types <- unique (wt_profile$way)
    # clearer to code as a for loop
    for (i in seq (way_types)) {
        if (!way_types [i] %in% dodgr_types) {
            pos <- which (pmatch (dodgr_types, way_types [i]) > 0)
            if (length (pos) > 0) {
                sf_lines$highway [sf_lines$highway == way_types [i]] <-
                    dodgr_types [pos]
            }
        }
    }

    # re-map some common types
    indx <- grep ("pedestrian|footway", sf_lines$highway)
    if (length (indx) > 0) {
        sf_lines$highway [indx] <- "path"
    }

    way_types <- unique (as.character (sf_lines$highway))
    not_in_wt_prof <- way_types [which (!way_types %in% dodgr_types)]
    if (length (not_in_wt_prof) > 0) {
        message (
            "The following highway types are present in data yet ",
            "lack corresponding weight_profile values: ",
            paste0 (not_in_wt_prof, sep = ", ")
        )
    }

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
get_sf_geom_col <- function (graph) {

    gcol <- grep ("geom", names (graph))
    if (length (gcol) > 1) {
        gnames <- c ("geometry", "geom", "geoms")
        mg <- match (gnames, names (graph))
        if (length (which (!is.na (mg))) == 1) {
            gcol <- mg [which (!is.na (mg))]
        } else {
            stop (
                "Unable to determine geometry column from [",
                paste0 (names (graph) [gcol], collapse = ", "), "]"
            )
        }
    } else if (length (gcol) == 0) {
        gcol <- match ("g", names (graph))
        if (is.na (gcol) || length (gcol) != 1) {
            stop ("Unable to determine geometry column")
        }
    }

    return (names (graph) [gcol])
}

rownames_from_xy <- function (graph) {

    xyf <- data.frame (
        x = graph$from_lon,
        y = graph$from_lat
    )
    xyt <- data.frame (
        x = graph$to_lon,
        y = graph$to_lat
    )

    ids <- rcpp_unique_rownames (xyf, xyt, precision = 10)
    graph$from_id <- ids$from_id
    graph$to_id <- ids$to_id

    return (graph)
}

reinsert_keep_cols <- function (sf_lines, graph, keep_cols) {

    keep_names <- NULL
    if (is.character (keep_cols)) {
        keep_cols <- lapply (keep_cols, function (i) match (i, names (sf_lines)))
        keep_cols <- sort (unique (unlist (keep_cols)))
        keep_names <- names (sf_lines) [keep_cols]
        # NA is no keep_cols match
    } else if (is.numeric (keep_cols)) {
        if (min (keep_cols) < 1 || max (keep_cols) > nrow (sf_lines)) {
            stop (
                "Numeric keep_cols must index into columns of 'sf' input",
                call. = FALSE
            )
            keep_names <- names (sf_lines) [keep_cols]
        }
    } else {
        stop ("keep_cols must be either character or numeric", .call = FALSE)
    }
    index <- which (!is.na (keep_cols))
    keep_cols <- keep_cols [index]
    keep_names <- keep_names [index]
    if (length (keep_cols) > 0) {
        indx <- match (graph$geom_num, seq (sf_lines$geometry))
        for (k in seq_along (keep_cols)) {
            graph [[keep_names [k]]] <-
                sf_lines [indx, keep_cols [k], drop = TRUE]
        }
    }

    return (graph)
}

add_extra_sf_columns <- function (graph, x) {

    if (!"way_id" %in% names (graph)) { # only works for OSM data
        return (graph)
    } # nocov

    hi <- match ("highway", names (graph))
    if (is.na (hi)) {
        hi <- ncol (graph)
        index2 <- NULL
    } else if (hi == ncol (graph)) {
        index2 <- NULL # nocov
    } else {
        index2 <- (hi + 1):ncol (graph)
    }

    keep_types <- c ("lanes", "maxspeed", "surface")
    keep_df <- array (NA_character_,
        dim = c (nrow (graph), length (keep_types))
    )
    nms <- c (names (graph) [1:hi], keep_types, names (graph) [index2])
    graph <- cbind (
        graph [, 1:hi],
        data.frame (keep_df, stringsAsFactors = FALSE),
        graph [, index2]
    )
    names (graph) <- nms

    row_index <- match (graph$way_id, x$osm_id)
    col_index_x <- match (keep_types, names (x))
    keep_types <- keep_types [which (!is.na (col_index_x))]
    col_index_x <- col_index_x [which (!is.na (col_index_x))]
    col_index_graph <- match (keep_types, names (graph))

    x [[attr (x, "sf_column")]] <- NULL
    x <- data.frame (x, stringsAsFactors = FALSE)
    # that still sometimes produces factors, so:
    for (i in seq_len (ncol (x))) {
        x [, i] <- paste0 (x [, i])
    }
    graph [, col_index_graph] <- x [row_index, col_index_x]

    return (graph)
}

# ********************************************************************
# *************************     sc class     *************************
# ********************************************************************
#
# most functions are defined in weight-streetnet-times.R

#' @name weight_streetnet
#' @family extraction
#' @export
weight_streetnet.sc <- function (x,
                                 wt_profile = "bicycle",
                                 wt_profile_file = NULL,
                                 turn_penalty = FALSE,
                                 type_col = "highway",
                                 id_col = "osm_id",
                                 keep_cols = NULL,
                                 left_side = FALSE) {

    requireNamespace ("geodist")
    requireNamespace ("dplyr")
    check_sc (x)

    x$vertex <- x$vertex [which (!duplicated (x$vertex)), ]

    if (wt_profile == "bicycle") {
        if (is.integer (keep_cols)) {
            stop (
                "keep_cols for 'sc' networks must be names of columns, not indices",
                call. = FALSE
            )
        }
        keep_cols <- unique (c (keep_cols, c (
            "bicycle",
            "cycleway",
            "cycleway:left",
            "cycleway:right"
        )))
    }
    keep_cols <- unique (c (way_types_to_keep, keep_cols))

    graph <- extract_sc_edges_xy (x) %>% # vert, edge IDs + coordinates
        sc_edge_dist () %>% # append dist
        extract_sc_edges_highways (
            x,
            wt_profile,
            wt_profile_file,
            keep_cols
        ) %>% # hw key-val pairs
        weight_sc_edges (
            wt_profile,
            wt_profile_file
        ) %>% # add d_weighted col
        set_maxspeed (
            wt_profile,
            wt_profile_file
        ) %>% # modify d_weighted
        weight_by_num_lanes (wt_profile) %>%
        calc_edge_time (wt_profile) %>% # add time
        sc_traffic_lights (
            x,
            wt_profile,
            wt_profile_file
        ) %>% # modify time
        rm_duplicated_edges () %>%
        sc_duplicate_edges (wt_profile)

    cl <- class (graph)
    graph <- dodgr_components (graph) # strips tbl class
    class (graph) <- cl

    gr_cols <- dodgr_graph_cols (graph)
    graph <- graph [which (!is.na (graph [[gr_cols$d_weighted]])), ]

    attr (graph, "turn_penalty") <- 0

    if (turn_penalty) {
        attr (graph, "turn_penalty") <-
            get_turn_penalties (wt_profile, wt_profile_file)$turn
        attr (graph, "wt_profile") <- wt_profile
        attr (graph, "wt_profile_file") <- wt_profile_file
        attr (graph, "left_side") <- left_side

        restrictions <- extract_turn_restictions (x)
        attr (graph, "turn_restrictions_no") <- restrictions$rw_no
        attr (graph, "turn_restrictions_only") <- restrictions$rw_only
    }

    gr_cols <- dodgr_graph_cols (graph)
    graph <- graph [which (!is.na (graph [[gr_cols$d_weighted]]) |
        !is.na (graph [[gr_cols$time]])), ]

    class (graph) <- c (
        class (graph),
        "dodgr_streetnet",
        "dodgr_streetnet_sc"
    )

    attr (graph, "hash") <-
        get_hash (graph, contracted = FALSE, force = TRUE)

    if (is_dodgr_cache_on ()) {
        attr (graph, "px") <- cache_graph (graph, gr_cols$edge_id)
    }

    return (graph)
}

#' @name weight_streetnet
#' @family extraction
#' @export
weight_streetnet.SC <- function (x,
                                 wt_profile = "bicycle",
                                 wt_profile_file = NULL,
                                 turn_penalty = FALSE,
                                 type_col = "highway",
                                 id_col = "osm_id",
                                 keep_cols = NULL,
                                 left_side = FALSE) {

    weight_streetnet.sc (
        x,
        wt_profile = wt_profile,
        wt_profile_file = wt_profile_file,
        turn_penalty = turn_penalty,
        type_col = type_col,
        id_col = id_col,
        keep_cols = keep_cols,
        left_side = left_side
    )
}
