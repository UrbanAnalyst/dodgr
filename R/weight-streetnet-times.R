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

extract_sc_edges_highways <- function (edges, x, wt_profile,
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

    edges <- dplyr::left_join (edges, x$object_link_edge, by = "edge_") %>%
        dplyr::select (-native_)
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

    wp <- dodgr::weighting_profiles$weighting_profiles
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
        # increase weighted distance according to numbers of lanes
        lns <- c (4, 5, 6, 7, 8)
        wts <- c (0.05, 0.05, 0.1, 0.1, 0.2)
        for (i in seq (lns))
        {
            index <- which (edges$lanes == lns [i])
            if (i == length (lns))
                index <- which (edges$lanes >= lns [i])
            edges$d_weighted [index] <- edges$d_weighted [index] * (1 + wts [i])
        }
    }

    if ("lanes" %in% names (edges)) edges$lanes <- NULL
    if ("maxspeed" %in% names (edges)) edges$maxspeed <- NULL
    #edges <- dplyr::select (edges, -c(lanes, maxspeed))

    return (edges)
}

# Convert weighted distances in metres to time in seconds
sc_edge_time <- function (edges, wt_profile, x)
{
    # no visible binding messages:
    object_ <- NULL

    edges$time <- NA_real_

    # General times based on speed profile
    w <- dodgr::weighting_profiles$weighting_profiles
    w <- w [w$name == wt_profile, ]
    # speeds are in km/h, but distances are in m. speeds in km/s = s / 3600,
    # and speeds in m/s = s * 1000 / 3600
    for (i in seq (nrow (w)))
    {
        speed_m_per_s <- w$speeds [i] * 1000 / 3600 # m/h -> m/s
        index <- which (edges$highway == w$way [i])
        edges$time [index] <- edges$d [index] / speed_m_per_s
    }

    if (wt_profile %in% c ("foot", "bicycle"))
    {
        if ("dz" %in% names (edges))
            edges <- times_by_incline (edges, wt_profile)

        speeds <- dodgr::weighting_profiles$surface_speeds
        speeds <- speeds [speeds$name == wt_profile, ]
        for (i in seq (nrow (speeds)))
        {
            index <- which (edges [[speeds$key [i] ]] ==
                            edges [[speeds$value [i] ]])
            if (length (index) > 0)
                edges$time [index] <- edges$time [index] / speeds$speed [i]
        }
    }

    return (edges)
}

times_by_incline <- function (edges, wt_profile)
{
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
}

sc_traffic_lights <- function (edges, wt_profile, x)
{
    # no visible binding NOTES:
    object_ <- NULL

    w <- dodgr::weighting_profiles$penalties
    wait <- w$traffic_lights [w$name == wt_profile]

    # first for intersections marked as crossings
    crossings <- traffic_light_objs_ped (x) # way IDs
    objs <- x$object %>% dplyr::filter (object_ %in% crossings$crossings)
    oles <- x$object_link_edge %>% dplyr::filter (object_ %in% objs$object_)
    index <- match (oles$edge_, edges$edge_)
    edges$time [index] <- edges$time [index] + wait

    # then all others with nodes simply marked as traffic lights - match
    # those to *start* nodes and simply add the waiting time
    nodes <- traffic_signal_nodes (x)
    index2 <- which (edges$.vx0 %in% nodes &
                     !edges$.vx0 %in% edges$.vx0 [index])
    edges$time [index2] <- edges$time [index2] + wait
    
    return (edges)
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
