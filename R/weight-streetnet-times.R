# definitions used in weight_streetnet.sc, including functions for time-based
# network weighting.

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

sc_edge_dist <- function (graph)
{
    # no visible binding notes:
    .vx0_z <- .vx1_z <- NULL

    graph$d <- geodist::geodist (graph [, c (".vx0_x", ".vx0_y")],
                                 graph [, c (".vx1_x", ".vx1_y")], paired = TRUE)
    if (".vx0_z" %in% names (graph) & ".vx1_z" %in% names (graph))
        graph <- dplyr::mutate (graph, "dz" = .vx1_z - .vx0_z) %>%
            dplyr::select (-c(.vx0_z, .vx1_z))
    return (graph)
}

extract_sc_edges_highways <- function (graph, x, wt_profile,
    keep_types = way_types_to_keep)
{
    # no visible binding notes:
    native_ <- key <- `:=` <- value <- NULL

    surface <- dodgr::weighting_profiles$surface_speeds
    surface <- surface [surface$name == wt_profile, ]
    if (nrow (surface) > 0)
    {
        keep_types <- c (keep_types, unique (surface$key))
    }

    graph <- dplyr::left_join (graph, x$object_link_edge, by = "edge_") %>%
        dplyr::select (-native_)
    for (k in keep_types)
    {
        objs <- dplyr::filter (x$object, key == k)
        graph <- dplyr::left_join (graph, objs, by = "object_") %>%
            dplyr::rename (!!dplyr::quo_name (k) := value) %>%
            dplyr::select (-key)
    }

    graph <- convert_hw_types_to_bool (graph)

    # TODO: Do this in convert_graph function
    # replace NA distances and times
    #m <- .Machine$double.xmax
    #graph$d_weighted [is.na (graph$d_weighted)] <- m
    #graph$time [is.na (graph$time)] <- m

    return (graph)
}

convert_hw_types_to_bool <- function (graph)
{
    # oneway:bicycle doesn't enquote properly, so:
    names (graph) [grep ("bicycle", names (graph))] <- "oneway_bicycle"

    # re-map the oneway values to boolean
    graph$oneway [!graph$oneway %in% c ("no", "yes")] <- "no"
    graph$oneway <- ifelse (graph$oneway == "no", FALSE, TRUE)
    graph$oneway_bicycle [!graph$oneway_bicycle %in% c ("no", "yes")] <- "no"
    graph$oneway_bicycle <- ifelse (graph$oneway_bicycle == "no", FALSE, TRUE)
    return (graph)
}

weight_sc_edges <- function (graph, wt_profile)
{
    # no visible binding notes:
    value <- d <- NULL

    wp <- dodgr::weighting_profiles$weighting_profiles
    wp <- wp [wp$name == wt_profile, c ("way", "value")]
    dplyr::left_join (graph, wp, by = c ("highway" = "way")) %>%
        dplyr::filter (!is.na (value)) %>%
        dplyr::mutate (d_weighted = ifelse (value == 0, NA, d / value)) %>%
        dplyr::select (-value)
}

# adjust weighted distances according to numbers of lanes and surfaces
# NOTE: This is currently hard-coded for active transport only, and will not
# work for anything else
wt_lanes_surface <- function (graph, wt_profile)
{
    # no visible binding notes:
    lanes <- maxspeed <- surface <- NULL

    if (wt_profile %in% c ("foot", "bicycle"))
    {
        # increase weighted distance according to numbers of lanes
        lns <- c (4, 5, 6, 7, 8)
        wts <- c (0.05, 0.05, 0.1, 0.1, 0.2)
        for (i in seq (lns))
        {
            index <- which (graph$lanes == lns [i])
            if (i == length (lns))
                index <- which (graph$lanes >= lns [i])
            graph$d_weighted [index] <- graph$d_weighted [index] * (1 + wts [i])
        }
    }

    if ("lanes" %in% names (graph)) graph$lanes <- NULL
    if ("maxspeed" %in% names (graph)) graph$maxspeed <- NULL
    #graph <- dplyr::select (graph, -c(lanes, maxspeed))

    return (graph)
}

# Convert weighted distances in metres to time in seconds
sc_edge_time <- function (graph, wt_profile, x)
{
    # no visible binding messages:
    object_ <- NULL

    graph$time <- NA_real_

    # General times based on speed profile
    w <- dodgr::weighting_profiles$weighting_profiles
    w <- w [w$name == wt_profile, ]
    # speeds are in km/h, but distances are in m. speeds in km/s = s / 3600,
    # and speeds in m/s = s * 1000 / 3600
    for (i in seq (nrow (w)))
    {
        speed_m_per_s <- w$speeds [i] * 1000 / 3600 # m/h -> m/s
        index <- which (graph$highway == w$way [i])
        graph$time [index] <- graph$d [index] / speed_m_per_s
    }

    if (wt_profile %in% c ("foot", "bicycle"))
    {
        if ("dz" %in% names (graph))
            graph <- times_by_incline (graph, wt_profile)

        speeds <- dodgr::weighting_profiles$surface_speeds
        speeds <- speeds [speeds$name == wt_profile, ]
        for (i in seq (nrow (speeds)))
        {
            index <- which (graph [[speeds$key [i] ]] ==
                            graph [[speeds$value [i] ]])
            if (length (index) > 0)
                graph$time [index] <- graph$time [index] / speeds$speed [i]
        }
    }

    return (graph)
}

times_by_incline <- function (graph, wt_profile)
{
    if (wt_profile == "foot")
    {
        # Uses 
        # [Naismith's Rule](https://en.wikipedia.org/wiki/Naismith%27s_rule)
        graph$time = 3600 * (graph$d_weighted / 1000) / 5 # 5km per hour
        if ("dz" %in% names (graph))
        {
            index <- which (graph$dz > 0)
            graph$time [index] <- graph$time [index] + graph$dz [index] / 10
            graph$dz <- NULL
        }

    } else if (wt_profile == "bicycle")
    {
        # http://theclimbingcyclist.com/gradients-and-cycling-how-much-harder-are-steeper-climbs/
        # http://cycleseven.org/effect-of-hills-on-cycling-effort
        # The latter argues for a linear relationship with a reduction in speed
        # of "about 11% for every 1% change in steepness". For 0.01 to translate
        # to 0.11, it needs to be multiplied by 0.11 / 0.01, or 11
        graph$time = 3600 * (graph$d_weighted / 1000) / 12 # 12km per hour
        if ("dz" %in% names (graph))
        {
            index <- which (graph$dz > 0)
            graph$time [index] <- graph$time [index] * 
                (1 + 11 * graph$dz [index] / graph$d [index])
            graph$dz <- NULL
        }
        # ... TODO: Downhill
        # http://www.sportsci.org/jour/9804/dps.html
        # downhill cycling speed ~ sqrt (slope)
    }
}

sc_traffic_lights <- function (graph, wt_profile, x)
{
    # no visible binding NOTES:
    object_ <- NULL

    w <- dodgr::weighting_profiles$penalties
    wait <- w$traffic_lights [w$name == wt_profile]

    # first for intersections marked as crossings
    crossings <- traffic_light_objs_ped (x) # way IDs
    objs <- x$object %>% dplyr::filter (object_ %in% crossings$crossings)
    oles <- x$object_link_edge %>% dplyr::filter (object_ %in% objs$object_)
    index <- match (oles$edge_, graph$edge_)
    graph$time [index] <- graph$time [index] + wait

    # then all others with nodes simply marked as traffic lights - match
    # those to *start* nodes and simply add the waiting time
    nodes <- traffic_signal_nodes (x)
    index2 <- which (graph$.vx0 %in% nodes &
                     !graph$.vx0 %in% graph$.vx0 [index])
    graph$time [index2] <- graph$time [index2] + wait
    
    return (graph)
}

rm_duplicated_edges <- function (graph)
{
    index <- cbind (which (duplicated (graph [, c (".vx0", ".vx1")])),
                    which (duplicated (graph [, c (".vx0", ".vx1")], fromLast = TRUE)))
    removes <- apply (index, 1, function (i)
                      ifelse (graph$time [i [1] ] > graph$time [i [2] ],
                              i [1], i [2]))
    graph [!seq (nrow (graph)) %in% removes, ]
}

# up to that point, all edges are non-duplicated, and so need to be duplicated
# for non-oneway
sc_duplicate_edges <- function (x, wt_profile)
{
    oneway_modes <-  c ("motorcycle", "motorcar", "goods", "hgv", "psv")

    index <- seq (nrow (x))
    if (wt_profile %in% c( "bicycle", "moped"))
        index <- which (!x$oneway_bicycle)
    else if (wt_profile %in% oneway_modes)
        index <- which (!x$oneway)

    xnew <- x [index, ]
    xnew <- swap_cols (xnew, ".vx0", ".vx1")
    xnew <- swap_cols (xnew, ".vx0_x", ".vx1_x")
    xnew <- swap_cols (xnew, ".vx0_y", ".vx1_y")
    xnew$edge_ <- rcpp_gen_hash (nrow (xnew), 10)
    rbind (x, xnew)
}

swap_cols <- function (x, cola, colb)
{
    temp <- x [[cola]]
    x [[cola]] <- x [[colb]]
    x [[colb]] <- temp
    return (x)
}

get_turn_penalty <- function (wt_profile)
{
    p <- dodgr::weighting_profiles$penalties
    p$turn [p$name == wt_profile]
}


# traffic lights for pedestrians
# https://wiki.openstreetmap.org/wiki/Tag:highway%3Dtraffic_signals#Complex_intersections

# return silicate "object" instances -> OSM ways IDs asosicated with given sets
# of key-val pairs
get_key_val_pair <- function (x, kv)
{
    # no visible binding notes:
    key <- value <- object_ <- NULL

    xo <- lapply (kv, function (i)
                  dplyr::filter (x$object, key == i [1], value == i [2]) %>%
                      dplyr::select (object_) %>%
                      dplyr::pull (object_))
    xo <- table (do.call (c, xo))

    res <- NULL
    if (any (xo == length (kv)))
        res <- names (xo) [which (xo == length (kv))]

    return (res)
}

get_key_val_pair_node <- function (x, kv)
{
    # no visible binding notes:
    key <- value <- vertex_ <- NULL

    if (is.null (x$nodes))
        return (NULL)

    xo <- lapply (kv, function (i)
                  dplyr::filter (x$nodes, key == i [1], value == i [2]) %>%
                      dplyr::select (vertex_) %>%
                      dplyr::pull (vertex_))
    unique (unlist (xo))
}

# Get all OSM way IDs associated with traffic lights from osmdata_sc object x
traffic_light_objs_ped <- function (x)
{
    # 1. Traffic signal without intersection (e.g. before bridge), no pedestrian
    # crossing
    x1 <- get_key_val_pair (x, list (c ("highway", "traffic_signals"),
                                     c ("crossing", "no")))

    # 2. Pedestrian crossing without intersection
    x2a <- get_key_val_pair (x, list (c ("highway", "crossing"),
                                      c ("crossing", "traffic_signals")))
    x2b <- get_key_val_pair (x, list (c ("highway", "traffic_signals"),
                                      c ("crossing", "traffic_signals")))
    x2_forw <- get_key_val_pair (x, list (c ("highway", "traffic_signals"),
                                          c ("crossing", "no"),
                                  c ("traffic_signals:direction", "forward")))
    x2_back <- get_key_val_pair (x, list (c ("highway", "traffic_signals"),
                                          c ("crossing", "no"),
                                  c ("traffic_signals:direction", "backward")))

    # 3. Simple Intersection
    x3a <- get_key_val_pair (x, list (c ("highway", "traffic_signals"),
                                      c ("crossing", "traffic_signals")))
    x3b <- get_key_val_pair (x, list (c ("highway", "traffic_signals"),
                                      c ("crossing", "no")))
    x3c <- get_key_val_pair (x, list (c ("highway", "crossing"),
                                      c ("crossing", "traffic_signals")))

    # 4. Intersection of divided and undivided highway with no lights
    crossings <- get_key_val_pair (x, list (c ("highway", "footway"),
                                            c ("footway", "crossing")))

    xout <- unique (c (x1, x2a, x2b, x3a, x3b, x3c))
    list ("both" = xout,
          "forward" = x2_forw,
          "back" = x2_back,
          "crossings" = crossings)
}

# Get all OSM node IDs that are traffic lights from osmdata_sc object x
traffic_signal_nodes <- function (x)
{
    x1 <- get_key_val_pair_node (x, list (c ("highway", "traffic_signals")))
    x2 <- get_key_val_pair_node (x, list (c ("highway", "crossing"),
                                          c ("crossing", "traffic_signals")))
    unique (c (x1, x2))
}

join_junctions_to_graph <- function (graph, wt_profile, left_side = FALSE)
{
    turn_penalty <- get_turn_penalty (wt_profile)
    resbind <- NULL
    if (turn_penalty > 0)
    {
        res <- rcpp_route_times (graph, left_side, turn_penalty)

        index <- which (graph$.vx0 %in% res$junction_vertices)
        v_start <- graph$.vx0 [index]
        graph$.vx0 [index] <- paste0 (graph$.vx0 [index], "_start")
        index <- which (graph$.vx1 %in% res$junction_vertices)
        v_end <- graph$.vx1 [index]
        graph$.vx1 [index] <- paste0 (graph$.vx1 [index], "_end")

        # pad out extra columns of res to match any extra in original graph
        resbind <- data.frame (array (NA, dim = c (nrow (res$graph), ncol (graph))))
        names (resbind) <- names (graph)
        resbind [, which (names (graph) %in% names (res$graph))] <- res$graph
    }

    rbind (graph, resbind)
}
