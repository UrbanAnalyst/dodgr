#' Convert a `dodgr` graph into an equivalent \pkg{sf} object.
#'
#' Works by aggregating edges into `LINESTRING` objects representing longest
#' sequences between all junction nodes. The resultant objects will generally
#' contain more `LINESTRING` objects than the original \pkg{sf} object, because
#' the former will be bisected at every junction point.
#'
#' @param graph A `dodgr` graph
#' @return Equivalent object of class \pkg{sf}.
#'
#' @note Requires the \pkg{sf} package to be installed.
#'
#' @family conversion
#' @export
#' @examples
#' hw <- weight_streetnet (hampi)
#' nrow (hw) # 5,729 edges
#' xy <- dodgr_to_sf (hw)
#' dim (xy) # 764 edges; 14 attributes
dodgr_to_sf <- function (graph) {

    if (methods::is (graph, "dodgr_contracted")) {
        warning (
            "'dodgr_to_sf' is intended to be called only on ",
            "non-contracted graphs.\nCalling on a contracted ",
            "graph will result in loss of information."
        )
    }

    graph <- tbl_to_df (graph)
    attr (graph, "hash") <- get_hash (graph, contracted = FALSE, force = TRUE)

    requireNamespace ("sf")
    res <- dodgr_to_sfc (graph)
    sf::st_sf (res$dat, geometry = res$geometry, crs = 4326)
}

#' Convert a `dodgr` graph into an equivalent `sf::sfc` object.
#'
#' Convert a `dodgr` graph into a `list` composed of
#' two objects: `dat`, a `data.frame`; and
#' `geometry`, an `sfc` object from the (\pkg{sf}) package.
#' Works by aggregating edges into `LINESTRING`
#' objects representing longest sequences between all junction nodes. The
#' resultant objects will generally contain more `LINESTRING` objects than
#' the original \pkg{sf} object, because the former will be bisected at every
#' junction point.
#'
#' @param graph A `dodgr` graph
#' @return A list containing (1) A `data.frame` of data associated with the
#' `sf` geometries; and (ii) A Simple Features Collection (`sfc`) list of
#' `LINESTRING` objects.
#'
#' @note The output of this function corresponds to the edges obtained from
#' `dodgr_contract_graph`. This function does not require the \pkg{sf} package
#' to be installed; the corresponding function that creates a full \pkg{sf}
#' object - \link{dodgr_to_sf} does requires \pkg{sf} to be installed.
#'
#' @family conversion
#' @export
#' @examples
#' hw <- weight_streetnet (hampi)
#' nrow (hw)
#' xy <- dodgr_to_sfc (hw)
#' dim (hw) # 5.845 edges
#' length (xy$geometry) # more linestrings aggregated from those edges
#' nrow (hampi) # than the 191 linestrings in original sf object
#' dim (xy$dat) # same number of rows as there are geometries
#' # The dodgr_to_sf function then just implements this final conversion:
#' # sf::st_sf (xy$dat, geometry = xy$geometry, crs = 4326)
dodgr_to_sfc <- function (graph) {

    graph <- tbl_to_df (graph)

    # force sequential IDs. TODO: Allow non-sequential by replacing indices in
    # src/dodgr_to_sf::get_edge_to_vert_maps with maps to sequential indices.
    gr_cols <- dodgr_graph_cols (graph)
    graph [[gr_cols$edge_id]] <- seq_len (nrow (graph))

    # hard-code column names for Rcpp routine:
    names (graph) [gr_cols$edge_id] <- "edge_id"
    names (graph) [gr_cols$from] <- "from_id"
    names (graph) [gr_cols$to] <- "to_id"

    names (graph) [gr_cols$xfr] <- "from_lon"
    names (graph) [gr_cols$yfr] <- "from_lat"
    names (graph) [gr_cols$xto] <- "to_lon"
    names (graph) [gr_cols$yto] <- "to_lat"

    gc <- dodgr_contract_graph (graph, nocache = TRUE)
    edge_map <- get_edge_map (gc)
    if (nrow (edge_map) == 0L) {
        stop (
            "Graph has already been contracted, and can not be ",
            "converted to 'sf'",
            call. = FALSE
        )
    }
    geometry <- rcpp_aggregate_to_sf (graph, gc, edge_map)

    # geometry has full WKT representation of CRS via osmdata, but still
    # triggers "old-style crs" message from `sf`. This uses `sf` to replace the
    # crs with identical values which don't trigger that message.
    out <- suppressMessages (
        sf::st_crs (geometry) <- 4326
    )

    # Then match data of `graph` potentially including way_id, back on to the
    # geometries:
    # edge_ids <- gc$graph$edge_id [match (names (geometry), gc$graph$edge_id)]
    # indx1 <- which (edge_ids %in% gc$edge_map$edge_new)
    # indx2 <- seq (edge_ids) [!seq (edge_ids) %in% indx1]
    # edge_ids <- c (gc$edge_map$edge_old [indx1], edge_ids [indx2])
    # index <- match (edge_ids, graph$edge_id)
    # dat <- graph [index, ]
    # dat$from_id <- dat$from_lat <- dat$from_lon <- NULL
    # dat$to_id <- dat$to_lat <- dat$to_lon <- NULL
    # dat$d <- dat$d_weighted <- dat$edge_id <- NULL

    geometry <- geometry [match (gc [[gr_cols$edge_id]], names (geometry))]

    return (list (dat = gc, geometry = geometry))
}

#' Convert a `dodgr` graph to an \pkg{igraph}.
#'
#' @param graph A `dodgr` graph
#' @param weight_column The column of the `dodgr` network to use as the edge
#' weights in the `igraph` representation.
#'
#' @return The `igraph` equivalent of the input. Note that this will \emph{not}
#' be a dual-weighted graph.
#'
#' @seealso \link{igraph_to_dodgr}
#'
#' @family conversion
#' @export
#' @examples
#' graph <- weight_streetnet (hampi)
#' graphi <- dodgr_to_igraph (graph)
dodgr_to_igraph <- function (graph, weight_column = "d") {

    requireNamespace ("igraph")
    graph <- tbl_to_df (graph)
    if (!weight_column %in% names (graph)) {
        stop ("graph contains no column named '", weight_column, "'")
    }

    gr_cols <- dodgr_graph_cols (graph)
    if (is.na (gr_cols$from) || is.na (gr_cols$to)) {
        scols <- find_spatial_cols (graph)
        graph$from_id <- scols$xy_id$xy_fr_id
        graph$to_id <- scols$xy_id$xy_to_id
        gr_cols <- dodgr_graph_cols (graph)
    }

    v <- dodgr_vertices (graph)
    graph <- graph [, do.call (c, gr_cols [!is.na (gr_cols)])]
    names (graph) <- names (gr_cols) [!is.na (gr_cols)]
    gr_cols <- dodgr_graph_cols (graph)
    names (graph) [which (names (gr_cols) == weight_column)] <- "weight"

    # remove edge_id if it exists
    if (!is.na (gr_cols$edge_id)) {
        graph [[gr_cols$edge_id]] <- NULL
    }

    igraph::graph_from_data_frame (graph, directed = TRUE, vertices = v)
}

#' Convert a \pkg{igraph} network to an equivalent `dodgr` representation.
#'
#' @param graph An \pkg{igraph} network
#'
#' @return The `dodgr` equivalent of the input.
#'
#' @seealso \link{dodgr_to_igraph}
#'
#' @family conversion
#' @export
#' @examples
#' graph <- weight_streetnet (hampi)
#' graphi <- dodgr_to_igraph (graph)
#' graph2 <- igraph_to_dodgr (graphi)
#' identical (graph2, graph) # FALSE
igraph_to_dodgr <- function (graph) {

    requireNamespace ("igraph")
    ei <- igraph::edge_attr (graph)
    vi <- igraph::vertex_attr (graph)
    index <- grep ("^lon|^lat|lon$|lat$|^x|^y|x$|y$", names (ei))
    if (length (index) > 0) {
        ei [index] <- NULL
    }
    index <- which (names (vi) %in% names (ei))
    if (length (index) > 0) {
        vi [index] <- NULL
    }
    vi <- data.frame (do.call (cbind, vi), stringsAsFactors = FALSE)

    res <- data.frame (
        cbind (
            igraph::as_edgelist (graph),
            do.call (cbind, ei)
        ),
        stringsAsFactors = FALSE
    )
    names (res) [1:2] <- c ("from_id", "to_id")

    nms <- paste0 (names (vi), "_from") [-1]
    res [, nms] <- vi [match (res$from_id, vi$name), 2:ncol (vi)]

    nms <- paste0 (names (vi), "_to") [-1]
    res [, nms] <- vi [match (res$to_id, vi$name), 2:ncol (vi)]

    for (i in 3:ncol (res)) {
        res [, i] <- convert_col (res, i)
    }

    res <- cbind ("edge_id" = seq_len (nrow (res)), res)

    return (res)
}

# Convert columns from character to numeric or integer where those given
# identical round-trip values (char -> numeric -> char, for example).
convert_col <- function (x, n = 3) {
    xn <- xnd <- x [, n]
    cl0 <- class (xn)
    storage.mode (xnd) <- "numeric"
    xni <- round (xnd)
    storage.mode (xni) <- "integer"
    if (identical (methods::as (xni, cl0), xn)) {
        xn <- xni
    } else if (identical (methods::as (xnd, cl0), xn)) {
        xn <- xnd
    }
    return (xn)
}


#' Convert a `dodgr` graph to an \pkg{tidygraph}.
#'
#' @param graph A `dodgr` graph
#'
#' @return The `tidygraph` equivalent of the input
#' @family conversion
#' @export
#' @examples
#' graph <- weight_streetnet (hampi)
#' grapht <- dodgr_to_tidygraph (graph)
dodgr_to_tidygraph <- function (graph) {

    requireNamespace ("igraph")
    requireNamespace ("tidygraph")

    dodgr_to_igraph (graph) %>%
        tidygraph::as_tbl_graph ()
}
