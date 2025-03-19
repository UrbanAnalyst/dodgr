# Miscellaneous **non-exported** graph functions for spatial graphs

#' preprocess spatial columns of graph
#'
#' @noRd
preprocess_spatial_cols <- function (graph) {

    gr_cols <- dodgr_graph_cols (graph)
    if (is.na (gr_cols$from) || is.na (gr_cols$to)) {
        scols <- find_spatial_cols (graph)
        graph$from_id <- scols$xy_id$xy_fr_id
        graph$to_id <- scols$xy_id$xy_to_id
    }
    return (graph)
}

#' is_graph_spatial
#'
#' Is the graph spatial or not?
#' @param graph A `data.frame` of edges
#' @return `TRUE` is `graph` is spatial, otherwise `FALSE`
#' @noRd
is_graph_spatial <- function (graph) {

    ncol (graph) > 4 &
        (any (grepl (
            "^x|x$",
            names (graph) [find_fr_col (graph)],
            ignore.case = TRUE
        )) |
            any (grepl (
                "^y|y$",
                names (graph) [find_to_col (graph)],
                ignore.case = TRUE
            )) |
            any (grepl ("lon", names (graph), ignore.case = TRUE)) |
            any (grepl ("lat", names (graph), ignore.case = TRUE)))
}

#' find_spatial_cols
#'
#' @return `fr_col` and `to_col` as vectors of 2 values of `x`
#' then `y` coordinates
#'
#' @noRd
find_spatial_cols <- function (graph) {

    graph <- tbl_to_df (graph)

    fr_col <- find_fr_col (graph)
    to_col <- find_to_col (graph)

    if (length (fr_col) < 2 || length (to_col) < 2) {
        stop (paste0 (
            "Graph appears to be spatial yet unable to ",
            "extract coordinates."
        ))
    }

    if (length (fr_col) == 3) {
        frx_col <- find_xy_col (graph, fr_col, x = TRUE)
        fry_col <- find_xy_col (graph, fr_col, x = FALSE)
        frid_col <- fr_col [which (!fr_col %in% c (frx_col, fry_col))]
        fr_col <- c (frx_col, fry_col)
        xy_fr_id <- graph [, frid_col]
        if (!is.character (xy_fr_id)) {
            xy_fr_id <- paste0 (xy_fr_id)
        }
    } else { # len == 2, so must be only x-y
        if (length (grep ("lon|lat|x|y", names (graph) [fr_col])) != 2) {
            stop ("Unable to determine coordinate columns of graph")
        } # nocov
        xy_fr_id <- paste0 (
            graph [, fr_col [1]], "-",
            graph [, fr_col [2]]
        )
    }

    if (length (to_col) == 3) {
        tox_col <- find_xy_col (graph, to_col, x = TRUE)
        toy_col <- find_xy_col (graph, to_col, x = FALSE)
        toid_col <- to_col [which (!to_col %in% c (tox_col, toy_col))]
        to_col <- c (tox_col, toy_col)
        xy_to_id <- graph [, toid_col]
        if (!is.character (xy_to_id)) {
            xy_to_id <- paste0 (xy_to_id)
        }
    } else { # len == 2, so must be only x-y
        if (length (grep ("lon|lat|x|y", names (graph) [to_col])) != 2) {
            stop ("Unable to determine coordinate columns of graph")
        } # nocov
        xy_to_id <- paste0 (
            graph [, to_col [1]], "-",
            graph [, to_col [2]]
        )
    }

    list (
        fr_col = fr_col,
        to_col = to_col,
        xy_id = data.frame (
            xy_fr_id = xy_fr_id,
            xy_to_id = xy_to_id,
            stringsAsFactors = FALSE
        )
    )
}

max_spatial_dist <- function (graph) {

    sp_cols <- find_spatial_cols (graph)
    graph_xy <- graph [, c (sp_cols$fr_col, sp_cols$to_col)]
    # Rename because 'sc' have different names:
    names (graph_xy) <- c ("from_lon", "from_lat", "to_lon", "to_lat")

    rx <- range (c (
        range (graph_xy$from_lon),
        range (graph_xy$to_lon)
    ))
    ry <- range (c (
        range (graph_xy$from_lat),
        range (graph_xy$to_lat)
    ))
    xylims <- c (rx [1], ry [1], rx [2], ry [2])

    suppressMessages (
        geodist::geodist (
            xylims [1:2],
            xylims [3:4],
            measure = "haversine"
        ) [1, 1]
    )
}
