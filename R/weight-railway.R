#' Weight a network for routing along railways.
#'
#' Weight (or re-weight) an `sf`-formatted OSM street network for routing
#' along railways.
#'
#' @param x A street network represented either as `sf` `LINESTRING`
#' objects, typically extracted with \link{dodgr_streetnet}.
#' @param type_col Specify column of the `sf` `data.frame` object
#' which designates different types of railways to be used for weighting
#' (default works with `osmdata` objects).
#' @param id_col Specify column of the \pkg{sf} `data.frame` object which
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
#' @family extraction
#' @export
#' @examples
#' \dontrun{
#' # sample railway extraction with the 'osmdata' package
#' library (osmdata)
#' dat <- opq ("shinjuku") %>%
#'     add_osm_feature (key = "railway") %>%
#'     osmdata_sf (quiet = FALSE)
#' graph <- weight_railway (dat$osm_lines)
#' }
weight_railway <- function (x,
                            type_col = "railway",
                            id_col = "osm_id",
                            keep_cols = c ("maxspeed"),
                            excluded = c (
                                "abandoned",
                                "disused",
                                "proposed",
                                "razed"
                            )) {

    if (!is (x, "sf")) {
        stop ('x must be class "sf"')
    }
    geom_column <- get_sf_geom_col (x)
    attr (x, "sf_column") <- geom_column

    if (type_col != "railway") {
        names (x) [which (names (x) == type_col)] <- "railway"
    }
    if (id_col != "osm_id") {
        names (x) [which (names (x) == id_col)] <- "osm_id"
    } # nocov

    if (!"railway" %in% names (x)) {
        stop ("Please specify type_col to be used for weighting railway")
    }
    if (!"osm_id" %in% names (x)) {
        stop (
            "Please specifiy id_col to be used to identify ", # nocov
            "railway rows"
        )
    } # nocov

    requireNamespace ("geodist")

    if (is.null (names (x$geometry))) {
        names (x$geometry) <- x$osm_id
    } # nocov


    x <- x [which (!(x$railway %in% excluded | is.na (x$railway))), ]
    # routing is based on matching the given profile to the "highway" field of
    # x, so:
    x$highway <- x$railway

    wt_profile <- data.frame (
        name = "custom",
        way = unique (x$highway),
        value = 1,
        stringsAsFactors = FALSE
    )

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

    # rcpp_sf_as_network now flags non-routable ways with -1, so:
    graph$d_weighted [graph$d_weighted < 0] <- .Machine$double.xmax
    if (all (graph$highway == "")) {
        graph$highway <- NULL
    } # nocov
    if (all (graph$way_id == "")) {
        graph$way_id <- NULL
    } # nocov

    # If original geometries did not have rownames (meaning it's not from
    # osmdata), then reassign unique vertex from/to IDs based on coordinates
    if (is.null (rownames (as.matrix (x$geometry [[1]])))) {
        graph <- rownames_from_xy (graph)
    } # nocov

    # get component numbers for each edge
    class (graph) <- c (class (graph), "dodgr_streetnet")
    graph <- dodgr_components (graph)

    # And finally, re-insert keep_cols:
    if (length (keep_cols) > 0) {
        graph <- reinsert_keep_cols (x, graph, keep_cols)
    }

    return (graph)
}
