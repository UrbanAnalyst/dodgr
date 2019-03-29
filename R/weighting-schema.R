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
    get_key_val_pair_node (x, list (c ("highway", "traffic_signals")))
}
