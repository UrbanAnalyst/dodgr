# definitions used in weight_streetnet.sc, including functions for time-based
# network weighting.

has_elevation <- function (x) {

    "z_" %in% names (x$vertex)
}

check_sc <- function (x) {

    if (!"osmdata_sc" %in% class (x))
        stop ("weight_streetnet currently only works for 'sc'-class objects ",
              "extracted with osmdata::osmdata_sc.")
}

# First step of edge extraction: join x and y coordinates
extract_sc_edges_xy <- function (x) {

    rename0 <- c (.vx0_x = "x_", .vx0_y = "y_", .vx0_z = "z_")
    rename1 <- c (.vx1_x = "x_", .vx1_y = "y_", .vx1_z = "z_")
    if (!has_elevation (x)) {
        rename0 <- rename0 [1:2]
        rename1 <- rename1 [1:2]
    }

    dplyr::left_join (x$edge, x$vertex, by = c (".vx0" = "vertex_")) %>%
        dplyr::rename (!!rename0) %>%
        dplyr::left_join (x$vertex, by = c (".vx1" = "vertex_")) %>%
        dplyr::rename (!!rename1)
}

sc_edge_dist <- function (graph) {

    # no visible binding notes:
    .vx0_z <- .vx1_z <- NULL

    xy0 <- as.data.frame (graph [, c (".vx0_x", ".vx0_y")])
    xy1 <- as.data.frame (graph [, c (".vx1_x", ".vx1_y")])
    graph$d <- geodist::geodist (xy0, xy1, paired = TRUE)
    if (".vx0_z" %in% names (graph) & ".vx1_z" %in% names (graph))
        graph <- dplyr::mutate (graph, "dz" = .vx1_z - .vx0_z) %>%
            dplyr::select (-c(.vx0_z, .vx1_z))
    return (graph)
}

extract_sc_edges_highways <- function (graph, x, wt_profile, wt_profile_file,
                                       keep_cols) {

    # no visible binding notes:
    native_ <- key <- `:=` <- value <- NULL # nolint

    surface <- get_surface_speeds (wt_profile, wt_profile_file)
    if (nrow (surface) > 0) {
        keep_cols <- c (keep_cols, unique (surface$key))
    }

    graph <- dplyr::left_join (graph, x$object_link_edge, by = "edge_") %>%
        dplyr::select (-native_)
    for (k in keep_cols) {
        objs <- dplyr::filter (x$object, key == k)
        graph <- dplyr::left_join (graph, objs, by = "object_") %>%
            dplyr::rename (!!dplyr::quo_name (k) := value) %>%
            dplyr::select (-key)
    }

    convert_hw_types_to_bool (graph, wt_profile)
}

convert_hw_types_to_bool <- function (graph, wt_profile) {

    if (!"oneway" %in% names (graph))
        return (graph)
    if (is.logical (graph$oneway))
        return (graph)

    if (!(is.character (wt_profile) | is.data.frame (wt_profile)))
        return (graph)

    if (!is.character (wt_profile))
        wt_profile <- unique (wt_profile$name)
    b <- grep ("oneway.*bicycle|bicycle.*oneway", names(graph))
    if ("oneway" %in% names (graph) |
        (length (b) == 1 & wt_profile == "bicycle")) {
        index <- which (!graph$oneway %in% c ("no", "yes"))
        if (length (index) > 0)
            graph$oneway [index] <- "no"
        graph$oneway <- ifelse (graph$oneway == "no", FALSE, TRUE)

        if (length (b) == 1) {
            # oneway:bicycle doesn't enquote properly, so:
            names (graph) [b] <- "oneway_bicycle"

            index <- which (!graph$oneway_bicycle %in% c ("no", "yes"))
            if (length (index) > 0)
                graph$oneway_bicycle [index] <- "no"
            graph$oneway_bicycle <-
                ifelse (graph$oneway_bicycle == "no", FALSE, TRUE)

            if (wt_profile == "bicycle") {
                graph$oneway <- graph$oneway_bicycle
                graph$oneway_bicycle <- NULL
            }
        }
    }
    return (graph)
}

weight_sc_edges <- function (graph, wt_profile, wt_profile_file) {

    # no visible binding notes:
    value <- d <- NULL

    wp <- get_profile (wt_profile, wt_profile_file)
    wp <- wp [, c ("way", "value")]

    dplyr::left_join (graph, wp, by = c ("highway" = "way")) %>%
        dplyr::filter (!is.na (value)) %>%
        dplyr::mutate (d_weighted = ifelse (value == 0, NA, d / value)) %>%
        dplyr::select (-value)
}

# Set maximum speed for each edge.
set_maxspeed <- function (graph, wt_profile, wt_profile_file) {

    if (!"maxspeed" %in% names (graph))
        graph$maxspeed <- NA_real_ # nocov
    if (!"highway" %in% names (graph))
        return (graph) # nocov

    maxspeed <- rep (NA_real_, nrow (graph))
    index <- grep ("mph", graph$maxspeed)
    maxspeed [index] <- as.numeric (gsub ("[^[:digit:]. ]", "",
                                            graph$maxspeed [index]))
    maxspeed [index] <- maxspeed [index] * 1.609344

    index <- seq (nrow (graph)) [!(seq (nrow (graph)) %in% index)]
    maxspeed_char <- graph$maxspeed [index] # character string
    # some maxspeeds have two values, where the 1st is generally the "default"
    # value. This gsub extracts those only:
    maxspeed_char <- gsub ("[[:punct:]].*$", "", maxspeed_char)
    # some (mostly Austria and Germany) have "maxspeed:walk" for living streets.
    # This has no numeric value, but is replaced here with 10km/h
    maxspeed_char <- gsub ("walk", "10", maxspeed_char)
    maxspeed_char <- gsub ("none", NA, maxspeed_char)
    index2 <- which (!(is.na (maxspeed_char) |
                       maxspeed_char == "" |
                       maxspeed_char == "NA"))
    index2 <- index2 [which (!grepl ("[[:alpha:]]", maxspeed_char [index2]))]

    maxspeed_numeric <- rep (NA_real_, length (index))
    maxspeed_numeric [index2] <- as.numeric (maxspeed_char [index2])
    maxspeed [index] <- maxspeed_numeric

    graph$maxspeed <- maxspeed
    # Those are the OSM values, which must then be combined with values
    # determined from the specified profile. The lowest value is ultimately
    # chosen.
    wp <- get_profile (wt_profile, wt_profile_file)

    wp_index <- match (graph$highway, wp$way)
    graph_index <- which (!is.na (wp_index))
    wp_index <- wp_index [graph_index]
    maxspeed <- cbind (graph$maxspeed, rep (NA, nrow (graph)))
    maxspeed [graph_index, 2] <- wp$max_speed [wp_index]
    graph$maxspeed <- apply (maxspeed, 1, function (i)
                             ifelse (all (is.na (i)),
                                     NA_real_,
                                     min (i, na.rm = TRUE)))

    na_highways <- wp$way [which (is.na (wp$max_speed))]
    graph$maxspeed [graph$highway %in% na_highways] <- NA_real_
    # Also set weighted distance for all these to NA:
    #gr_cols <- dodgr_graph_cols (graph)
    #graph [[gr_cols$d_weighted]] [graph$highway %in% na_highways] <- NA_real_

    if (wt_profile %in% c ("horse", "wheelchair") |
        !"surface" %in% names (graph))
        return (graph)

    # And then repeat for max speeds according to surface profiles
    s <- get_surface_speeds (wt_profile, wt_profile_file)
    s <- s [s$name == wt_profile, c ("key", "value", "max_speed")]
    surf_vals <- unique (graph$surface [graph$surface != "NA"])
    surf_speeds <- s$max_speed [match (surf_vals, s$value)]
    surf_vals <- surf_vals [!is.na (surf_speeds)]
    surf_speeds <- surf_speeds [!is.na (surf_speeds)]

    surf_index <- match (graph$surface, surf_vals)
    graph_index <- which (!is.na (surf_index))
    surf_index <- surf_index [graph_index]
    maxspeed <- cbind (as.numeric (graph$maxspeed),
                       rep (NA_real_, nrow (graph)))
    maxspeed [graph_index, 2] <- surf_speeds [surf_index]
    graph$maxspeed <- apply (maxspeed, 1, function (i)
                             ifelse (all (is.na (i)),
                                     NA_real_,
                                     min (i, na.rm = TRUE)))

    graph$surface <- NULL

    return (graph)
}

# adjust weighted distances according to numbers of lanes
weight_by_num_lanes <- function (graph, wt_profile) {

    # only weight these profiles:
    profile_names <- c ("foot", "bicycle", "wheelchair", "horse")
    if (!(wt_profile %in% profile_names | "lanes" %in% names (graph)))
        return (graph) # nocov

    lns <- c (4, 5, 6, 7, 8)
    wts <- c (0.05, 0.05, 0.1, 0.1, 0.2)
    for (i in seq (lns)) {
        index <- which (graph$lanes == lns [i])
        if (i == length (lns))
            index <- which (graph$lanes >= lns [i])
        graph$d_weighted [index] <- graph$d_weighted [index] * (1 + wts [i])
    }

    graph$lanes <- NULL

    return (graph)
}

# Convert distances in metres to time in seconds. Up to this point, distances
# have been weighted for type of way (via
# weighting_profiles$weighting_profiles), and there is a maxspeed column
# reflecting profile values plus effect of different surfaces.
# The time is distance scaled by maxspeed, and time_weighted is d_weighted
# scaled by maxspeed
calc_edge_time <- function (graph, wt_profile) {

    gr_cols <- dodgr_graph_cols (graph)
    speed_m_per_s <- graph$maxspeed * 1000 / 3600 # maxspeeds are km/hr
    graph$time <- graph [[gr_cols$d]] / speed_m_per_s
    graph$time_weighted <- graph [[gr_cols$d_weighted]] / speed_m_per_s

    if ("dz" %in% names (graph) &
        wt_profile %in% c ("foot", "bicycle")) {
            graph <- times_by_incline (graph, wt_profile)
    }
    graph$maxspeed <- NULL

    return (graph)
}

# increase both real and weighted times according to elevation increases:
times_by_incline <- function (graph, wt_profile) {

    if (wt_profile == "foot") {
        # Uses
        # [Naismith's Rule](https://en.wikipedia.org/wiki/Naismith%27s_rule)
        if ("dz" %in% names (graph)) {
            index <- which (graph$dz > 0)
            graph$time [index] <- graph$time [index] + graph$dz [index] / 10
            graph$time_weighted [index] <- graph$time_weighted [index] +
                graph$dz [index] / 10
        }

    } else if (wt_profile == "bicycle") {
        # http://theclimbingcyclist.com/gradients-and-cycling-how-much-harder-are-steeper-climbs/ # nolint
        # http://cycleseven.org/effect-of-hills-on-cycling-effort
        # The latter argues for a linear relationship with a reduction in speed
        # of "about 11% for every 1% change in steepness". For 0.01 to translate
        # to 0.11, it needs to be multiplied by 0.11 / 0.01, or 11
        if ("dz" %in% names (graph)) {
            index <- which (graph$dz > 0)
            graph$time [index] <- graph$time [index] *
                (1 + 11 * graph$dz [index] / graph$d [index])
            graph$time_weighted [index] <- graph$time_weighted [index] *
                (1 + 11 * graph$dz [index] / graph$d [index])
        }
        # ... TODO: Downhill
        # http://www.sportsci.org/jour/9804/dps.html
        # downhill cycling speed ~ sqrt (slope)
    }
    return (graph)
}

sc_traffic_lights <- function (graph, x, wt_profile, wt_profile_file) {

    # no visible binding NOTES:
    object_ <- NULL

    wait <- get_turn_penalties (wt_profile, wt_profile_file)$traffic_lights
    if (length (wait) == 0) wait <- 0

    # first for intersections marked as crossings
    crossings <- traffic_light_objs (x) # way IDs
    objs <- x$object %>% dplyr::filter (object_ %in% crossings$crossings)
    oles <- x$object_link_edge %>% dplyr::filter (object_ %in% objs$object_)
    # Then the actual nodes with the traffic lights
    nodes <- traffic_signal_nodes (x)
    # Increment waiting times for edges ending at those nodes
    index <- which (graph$edge_ %in% oles$edge_ &
                    graph$.vx1 %in% nodes)
    graph$time [index] <- graph$time [index] + wait

    # then all others with nodes simply marked as traffic lights - match
    # those to *start* nodes and simply add the waiting time
    index2 <- which (graph$.vx0 %in% nodes &
                     !graph$.vx0 %in% graph$.vx0 [index])
    graph$time [index2] <- graph$time [index2] + wait

    return (graph)
}

rm_duplicated_edges <- function (graph) {

    gr_cols <- dodgr_graph_cols (graph)
    ft <- graph [, c (gr_cols$from, gr_cols$to)]
    index <- cbind (which (duplicated (ft)),
                    which (duplicated (ft, fromLast = TRUE)))

    index <- cbind (which (duplicated (graph [, c (".vx0", ".vx1")])),
                    which (duplicated (graph [, c (".vx0", ".vx1")],
                                       fromLast = TRUE)))
    removes <- apply (index, 1, function (i)
                      ifelse (graph$time [i [1] ] > graph$time [i [2] ], # nolint
                              i [1], i [2]))
    graph [!seq (nrow (graph)) %in% removes, ]
}

# up to that point, all edges are non-duplicated, and so need to be duplicated
# for non-oneway
sc_duplicate_edges <- function (x, wt_profile) {

    oneway_modes <-  c ("bicycle", "moped", "motorcycle", "motorcar", "goods",
                        "hgv", "psv")

    index <- seq (nrow (x))
    if (wt_profile %in% oneway_modes)
        index <- which (!x$oneway)

    xnew <- x [index, ]
    xnew <- swap_cols (xnew, ".vx0", ".vx1")
    xnew <- swap_cols (xnew, ".vx0_x", ".vx1_x")
    xnew <- swap_cols (xnew, ".vx0_y", ".vx1_y")
    xnew$edge_ <- rcpp_gen_hash (nrow (xnew), 10)

    res <- rbind (x, xnew)
    res$oneway <- NULL

    return (res)
}

swap_cols <- function (x, cola, colb) {

    temp <- x [[cola]]
    x [[cola]] <- x [[colb]]
    x [[colb]] <- temp
    return (x)
}



# traffic lights for pedestrians
# https://wiki.openstreetmap.org/wiki/Tag:highway%3Dtraffic_signals#Complex_intersections # nolint

# return silicate "object" instances -> OSM ways IDs asosicated with given sets
# of key-val pairs
get_key_val_pair <- function (x, kv) {

    # no visible binding notes:
    key <- value <- object_ <- NULL

    xo <- lapply (kv, function (i)
                  dplyr::filter (x$object, key == i [1], value == i [2]) %>%
                  dplyr::select (object_) %>%
                  dplyr::pull (object_))
    xo <- table (do.call (c, xo))

    res <- NULL
    if (any (xo == length (kv)))
        res <- names (xo) [which (xo == length (kv))] # nocov - not tested

    return (res)
}

get_key_val_pair_node <- function (x, kv) {

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
traffic_light_objs <- function (x) {

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
traffic_signal_nodes <- function (x) {

    x1 <- get_key_val_pair_node (x, list (c ("highway", "traffic_signals")))
    x2 <- get_key_val_pair_node (x, list (c ("highway", "crossing"),
                                          c ("crossing", "traffic_signals")))
    unique (c (x1, x2))
}

join_junctions_to_graph <- function (graph, wt_profile, wt_profile_file,
                                     left_side = FALSE) {

    turn_penalty <- get_turn_penalties (wt_profile, wt_profile_file)$turn
    resbind <- edge_map <- NULL

    if (turn_penalty > 0) {

        res <- rcpp_route_times (graph, left_side, turn_penalty)
        edge_map <- data.frame ("edge" = res$graph$edge_,
                                "e_in" = res$graph$old_edge_in,
                                "e_out" = res$graph$old_edge_out,
                                stringsAsFactors = FALSE)
        res$graph$old_edge_in <- res$graph$old_edge_out <- NULL

        index <- which (graph$.vx0 %in% res$junction_vertices)
        #v_start <- graph$.vx0 [index]
        graph$.vx0 [index] <- paste0 (graph$.vx0 [index], "_start")
        index <- which (graph$.vx1 %in% res$junction_vertices)
        #v_end <- graph$.vx1 [index]
        graph$.vx1 [index] <- paste0 (graph$.vx1 [index], "_end")

        # pad out extra columns of res to match any extra in original graph
        resbind <- data.frame (array (NA, dim = c (nrow (res$graph),
                                                   ncol (graph))))
        names (resbind) <- names (graph)
        resbind [, which (names (graph) %in% names (res$graph))] <- res$graph
        graph <- rbind (graph, resbind)
    }
    list (graph = graph, edge_map = edge_map)
}

#' Remove turn restrictions
#'
#' @param x The original `sc` object which ,when generated from
#' `dodgr_streetnet_sc`, includes turn restriction data
#' @param graph The processed by not yet turn-contracted graph
#' @param res The result of `join_junctions_to_graph`, with turn-contracted
#' `graph` and `edge_map` components.
#' @noRd
remove_turn_restrictions <- function (x, graph, res) {

    rels <- x$relation_properties # x from restrictions query above!!
    restriction_rels <- rels [rels$key == "restriction", ]
    index <- which (x$relation_members$relation_ %in%
                    restriction_rels$relation_)
    restriction_ways <- x$relation_members [index, ]

    rr_no <- restriction_rels [grep ("^no\\_", restriction_rels$value), ]
    rr_only <- restriction_rels [grep ("^only\\_", restriction_rels$value), ]
    rw_no <- restriction_ways [restriction_ways$relation_ %in%
                               rr_no$relation_, ]
    rw_only <- restriction_ways [restriction_ways$relation_ %in%
                                 rr_only$relation_, ]

    r_to_df <- function (r) {
        r <- lapply (split (r, f = factor (r$relation_)),
                     function (i) c (i$relation_ [1],
                                     i$member [2],
                                     i$member [c (1, 3)]))
        r <- data.frame (do.call (rbind, r))
        names (r) <- c ("relation", "node", "from", "to")
        return (stats::na.omit (r))
    }
    rw_no <- r_to_df (rw_no)
    rw_only <- r_to_df (rw_only)

    index0 <- match (rw_no$node, graph$.vx1) # in-edges
    index1 <- match (rw_no$node, graph$.vx0) # out-edges
    in_edges <- graph$edge_ [index0 [which (!is.na (index0))]]
    out_edges <- graph$edge_ [index1 [which (!is.na (index1))]]
    index <- which (res$edge_map$e_in %in% in_edges &
                    res$edge_map$e_out %in% out_edges)
    no_turn_edges <- res$edge_map$edge [index]

    index0 <- match (rw_only$node, graph$.vx1) # in-edges
    index1 <- match (rw_only$node, graph$.vx0) # out-edges
    in_edges <- graph$edge_ [index0 [which (!is.na (index0))]]
    out_edges <- graph$edge_ [index1 [which (!is.na (index1))]]
    # index of turns to edges other than "only" turn edges, so also to edges which
    # are to be excluded:
    index <- which (res$edge_map$e_in %in% in_edges &
                    !res$edge_map$e_out %in% out_edges)
    no_turn_edges <- unique (c (no_turn_edges, res$edge_map$edge [index]))

    res$graph <- res$graph [which (!res$graph$edge_ %in% no_turn_edges), ]
    res$edge_map <- res$edge_map [which (!res$edge_map$edge %in% no_turn_edges), ]

    return (res)
}
