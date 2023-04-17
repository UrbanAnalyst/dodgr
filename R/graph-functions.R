#' Identify the essential columns of the graph table (data.frame, tibble,
#' whatever) to be analysed in the C++ routines.
#'
#' @param graph A `data.frame` containing the edges of the graph
#' @return A list of column numbers of `edge_id`, `from`,
#' `to`, `d`, `w`, `time`, `xfr`, `yfr`, `xto`, `yto`, and `component`, some of
#' which may be NA.
#'
#' @noRd
dodgr_graph_cols <- function (graph) {

    nms <- names (graph)
    component <- grep ("comp", nms) %>% null_to_na ()
    if (methods::is (graph, "dodgr_streetnet") &&
        !methods::is (graph, "dodgr_streetnet_sc") &&
        ncol (graph) >= 11) {
        # columns are always identically structured
        edge_id <- which (nms == "edge_id") %>% null_to_na ()
        fr_col <- which (nms %in% c ("from", "from_id")) %>% null_to_na ()
        to_col <- which (nms %in% c ("to", "to_id")) %>% null_to_na ()
        d_col <- which (nms == "d")
        w_col <- which (nms == "d_weighted")

        xfr <- which (nms %in% c ("xfr", "from_lon"))
        if (length (xfr) == 0) xfr <- NA
        yfr <- which (nms %in% c ("yfr", "from_lat"))
        if (length (yfr) == 0) yfr <- NA
        xto <- which (nms %in% c ("xto", "to_lon"))
        if (length (xto) == 0) xto <- NA
        yto <- which (nms %in% c ("yto", "to_lat"))
        if (length (yto) == 0) yto <- NA
    } else {
        edge_id <- grep ("edge_id|edge_$", nms) %>% null_to_na ()

        d_col <- find_d_col (graph)
        w_col <- find_w_col (graph)
        # sc ensures this never happens, so not covered
        if (length (w_col) == 0) {
            w_col <- d_col # nocov
        }

        fr_col <- find_fr_id_col (graph)
        to_col <- find_to_id_col (graph)

        xfr <- yfr <- xto <- yto <- NA
        # TODO: Modify for other complex but non-spatial types of graph
        if (is_graph_spatial (graph)) {
            spcols <- find_spatial_cols (graph)
            graph <- tbl_to_df (graph)

            fr_is_num <- vapply (spcols$fr_col, function (i) {
                is.numeric (graph [[i]])
            }, logical (1))
            to_is_num <- vapply (spcols$to_col, function (i) {
                is.numeric (graph [[i]])
            }, logical (1))
            if (!(all (fr_is_num) && all (to_is_num))) {
                stop (paste0 (
                    "graph appears to have non-numeric ",
                    "longitudes and latitudes"
                ))
            }

            xfr <- spcols$fr_col [1]
            yfr <- spcols$fr_col [2]
            xto <- spcols$to_col [1]
            yto <- spcols$to_col [2]
        } else {
            if (length (fr_col) != 1 && length (to_col) != 1) {
                stop ("Unable to determine from and to columns in graph")
            } # nolint # nocov
        }
    }

    time_col <- grep ("time", nms)
    if (length (time_col) != 1) {
        time_col <- grep ("time$", nms)
        if (length (time_col) != 1) {
            time_col <- NA
        }
    }
    timew_col <- grep ("time_w|timew|tw", nms)
    if (length (timew_col) != 1) {
        timew_col <- grep ("time_w|timew|^tw", nms)
        if (length (timew_col) != 1) {
            timew_col <- NA
        }
    }

    ret <- c (
        edge_id,
        fr_col,
        to_col,
        d_col,
        w_col,
        time_col,
        timew_col,
        xfr,
        yfr,
        xto,
        yto,
        component
    )
    names (ret) <- c (
        "edge_id",
        "from",
        "to",
        "d",
        "d_weighted",
        "time",
        "time_weighted",
        "xfr",
        "yfr",
        "xto",
        "yto",
        "component"
    )
    class (ret) <- c (class (ret), "graph_columns")

    # This is passed to many C++ routines, in which case it needs to be
    # converted to a vector (`do.call (c, gr_cols)`), and the R-style 1-indexeso
    # need to be converted to equivalent 0-indexed forms
    return (as.list (ret))
}

#' Convert graph to a standard form suitable for submission to C++ routines
#'
#' @noRd
convert_graph <- function (graph, gr_cols) {

    tp <- NULL
    if ("turn_penalty" %in% names (attributes (graph))) {
        tp <- attr (graph, "turn_penalty")
    }

    keep_cols <- c (
        "edge_id", "from", "to", "d", "d_weighted",
        "time", "time_weighted"
    )
    index <- do.call (c, gr_cols [keep_cols])
    index <- index [!is.na (index)]
    graph <- graph [, index]
    names (graph) <- names (index)

    if ("edge_id" %in% names (graph)) {
        graph$edge_id <- convert_to_char (graph$edge_id)
    }
    graph$from <- convert_to_char (graph$from)
    graph$to <- convert_to_char (graph$to)

    if (!"time_weighted" %in% names (graph)) {
        graph$time_weighted <- graph$time
    }

    if (!is.null (tp)) {
        attr (graph, "turn_penalty") <- tp
    }

    return (graph)
}

convert_to_char <- function (x) {

    if (!is.character (x)) x <- paste0 (x)
    return (x)
}

tbl_to_df <- function (graph) {

    if (methods::is (graph, "tbl")) {
        classes <- class (graph) [!grepl ("tbl", class (graph))]
        graph <- as.data.frame (graph)
        class (graph) <- classes
    }
    return (graph)
}



#' Extract vertices of graph, including spatial coordinates if included.
#'
#' @param graph A flat table of graph edges. Must contain columns labelled
#' `from` and `to`, or `start` and `stop`. May also contain
#' similarly labelled columns of spatial coordinates (for example
#' `from_x`) or `stop_lon`).
#' @return A `data.frame` of vertices with unique numbers (`n`).
#'
#' @note Values of `n` are 0-indexed
#'
#' @family misc
#' @export
#' @examples
#' graph <- weight_streetnet (hampi)
#' v <- dodgr_vertices (graph)
dodgr_vertices <- function (graph) {

    # vertices are calculated as background process, so wait if that's not
    # finished.
    if ("px" %in% names (attributes (graph))) {
        px <- attr (graph, "px")
        while (px$is_alive ()) {
            px$wait ()
        }
    }

    hash <- ifelse (methods::is (graph, "dodgr_contracted"),
        "hashc", "hash"
    )
    hash <- attr (graph, hash)
    # make sure rows of graph have not been changed
    gr_cols <- dodgr_graph_cols (graph)
    hashe <- digest::digest (graph [[gr_cols$edge_id]])
    if (!identical (hashe, hash)) {
        hash <- NULL
    }
    if (!is.null (hash)) {
        hashe_ref <- attr (graph, "hashe")
        hashe_ref <- ifelse (is.null (hashe_ref), "", hashe_ref)
        hashe <- digest::digest (graph [[gr_cols$edge_id]])
        if (hashe != hashe_ref) {
            hash <- hashe
        }
    } else {
        if (is.na (gr_cols$edge_id)) {
            hash <- "" # nocov
        } else {
            hash <- get_hash (graph, contracted = FALSE, force = TRUE)
        }
    }

    verts_up_to_date <- FALSE

    fname <- fs::path (fs::path_temp (), paste0 ("dodgr_verts_", hash, ".Rds"))
    if (hash != "" && fs::file_exists (fname)) {
        verts <- readRDS (fname)
        # Then double-check to ensure those vertices match current graph:
        gr_cols <- dodgr_graph_cols (graph)
        graph_verts <- unique (c (graph [[gr_cols$from]], graph [[gr_cols$to]]))
        verts_up_to_date <-
            (all (graph_verts %in% verts$id) && all (verts$id %in% graph_verts))
    }

    if (!verts_up_to_date) {
        verts <- dodgr_vertices_internal (graph)
        saveRDS (verts, fname)
    }

    return (verts)
}

dodgr_vertices_internal <- function (graph) {

    graph <- tbl_to_df (graph)

    gr_cols <- dodgr_graph_cols (graph)
    # cols are (edge_id, from, to, d, w, component, xfr, yfr, xto, yto)
    # NOTE: c (x, y), where x and y are both factors gives junk, so explicit
    # conversion required here: TODO: Find a better way?
    if (is.factor (graph [[gr_cols$from]])) {
        graph [[gr_cols$from]] <- paste0 (graph [[gr_cols$from]])
    }
    if (is.factor (graph [[gr_cols$to]])) {
        graph [[gr_cols$to]] <- paste0 (graph [[gr_cols$to]])
    }

    if (is_graph_spatial (graph)) {
        verts <- data.frame (
            id = c (
                graph [[gr_cols$from]],
                graph [[gr_cols$to]]
            ),
            x = c (
                graph [[gr_cols$xfr]],
                graph [[gr_cols$xto]]
            ),
            y = c (
                graph [[gr_cols$yfr]],
                graph [[gr_cols$yto]]
            ),
            stringsAsFactors = FALSE
        )
        if (!is.na (gr_cols$component)) {
            verts$component <- graph [[gr_cols$component]]
        }
    } else {
        verts <- data.frame (
            id = c (
                graph [[gr_cols$from]],
                graph [[gr_cols$to]]
            ),
            stringsAsFactors = FALSE
        )
        if (!is.na (gr_cols$component)) {
            verts$component <- graph [[gr_cols$component]]
        }
    }

    # The next line is the time-killer here, which is why this is cached
    indx <- which (!duplicated (verts$id))
    verts <- verts [indx, , drop = FALSE] # nolint
    verts$n <- seq_len (nrow (verts)) - 1

    return (verts)
}


#' Identify connected components of graph.
#'
#' Identify connected components of graph and add corresponding `component`
#' column to `data.frame`.
#'
#' @param graph A `data.frame` of edges
#' @return Equivalent graph with additional `component` column,
#' sequentially numbered from 1 = largest component.
#' @family modification
#' @export
#' @examples
#' graph <- weight_streetnet (hampi)
#' graph <- dodgr_components (graph)
dodgr_components <- function (graph) {

    graph <- tbl_to_df (graph)

    if ("component" %in% names (graph)) {
        message ("graph already has a component column")
    } else {
        gr_cols <- dodgr_graph_cols (graph)
        graph2 <- convert_graph (graph, gr_cols)
        if (is.na (gr_cols$edge_id)) {
            graph2$edge_id <- seq_len (nrow (graph2))
        }
        cns <- rcpp_get_component_vector (graph2)

        indx <- match (graph2$edge_id, cns$edge_id)
        component <- cns$edge_component [indx]
        # Then re-number in order to decreasing component size:
        graph$component <- match (
            component,
            order (table (component), decreasing = TRUE)
        )
    }

    return (graph)
}

#' Sample a random but connected sub-component of a graph
#'
#' @param graph A flat table of graph edges. Must contain columns labelled
#' `from` and `to`, or `start` and `stop`. May also contain
#' similarly labelled columns of spatial coordinates (for example
#' `from_x`) or `stop_lon`).
#' @param nverts Number of vertices to sample
#'
#' @return A connected sub-component of `graph`
#'
#' @note Graphs may occasionally have `nverts + 1` vertices, rather than
#' the requested `nverts`.
#'
#' @family misc
#' @export
#' @examples
#' graph <- weight_streetnet (hampi)
#' nrow (graph) # 5,742
#' graph <- dodgr_sample (graph, nverts = 200)
#' nrow (graph) # generally around 400 edges
#' nrow (dodgr_vertices (graph)) # 200
dodgr_sample <- function (graph, nverts = 1000) {

    graph <- tbl_to_df (graph)

    fr <- find_fr_id_col (graph)
    to <- find_to_id_col (graph)
    verts <- unique (c (graph [, fr], graph [, to]))
    gr_cols <- dodgr_graph_cols (graph)
    if (is.na (gr_cols$edge_id)) {
        graph$edge_id <- seq_len (nrow (graph))
        gr_cols <- dodgr_graph_cols (graph)
    }

    if (length (verts) > nverts) {
        graph2 <- convert_graph (graph, gr_cols)
        if ("component" %in% names (graph)) {
            graph2$component <- graph$component
        }
        indx <- match (rcpp_sample_graph (graph2, nverts), graph2$edge_id)
        graph <- graph [sort (indx), ]
    }

    attr (graph, "hash") <- get_hash (graph, contracted = FALSE, force = TRUE)

    return (graph)
}

#' Insert a new node or vertex into a network
#'
#' @param v1 Vertex defining start of graph edge along which new vertex is to be
#' inserted
#' @param v2 Vertex defining end of graph edge along which new vertex is to be
#' inserted (order of `v1` and `v2` is not important).
#' @param x The `x`-coordinate of new vertex. If not specified, vertex is
#' created half-way between `v1` and `v2`.
#' @param y The `y`-coordinate of new vertex. If not specified, vertex is
#' created half-way between `v1` and `v2`.
#' @return A modified graph with specified edge between defined start and end
#' vertices split into two edges either side of new vertex.
#' @inheritParams dodgr_vertices
#'
#' @family misc
#' @export
#' @examples
#' graph <- weight_streetnet (hampi)
#' e1 <- sample (nrow (graph), 1)
#' v1 <- graph$from_id [e1]
#' v2 <- graph$to_id [e1]
#' # insert new vertex in the middle of that randomly-selected edge:
#' graph2 <- dodgr_insert_vertex (graph, v1, v2)
#' nrow (graph)
#' nrow (graph2) # new edges added to graph2
dodgr_insert_vertex <- function (graph, v1, v2, x = NULL, y = NULL) {

    graph_t <- tbl_to_df (graph)
    gr_cols <- dodgr_graph_cols (graph_t)
    index12 <- which (graph [[gr_cols$from]] == v1 & graph [[gr_cols$to]] == v2)
    index21 <- which (graph [[gr_cols$from]] == v2 & graph [[gr_cols$to]] == v1)
    if (length (index12) == 0 && length (index21) == 0) {
        stop ("Nominated vertices do not define any edges in graph")
    }
    if ((!is.null (x) && is.null (y)) || (is.null (x) && !is.null (y))) {
        stop ("Either both x and y must be NULL, or both must be specified")
    }


    charvec <- c (letters, LETTERS, 0:9)
    randid <- function (charvec, len = 10) {
        paste0 (sample (charvec, len, replace = TRUE), collapse = "")
    }

    if (length (index12) == 1) {
        graph <- insert_one_edge (graph, index12, x, y, gr_cols)
        graph [index12, gr_cols$to] <-
            graph [index12 + 1, gr_cols$from] <- randid (charvec, 10)
        index21 <- which (graph [[gr_cols$from]] == v2 &
            graph [[gr_cols$to]] == v1)
    }
    if (length (index21) == 1) {
        graph <- insert_one_edge (graph, index21, x, y, gr_cols)
        graph [index21, gr_cols$to] <-
            graph [index21 + 1, gr_cols$from] <- randid (charvec, 10)
    }

    attr (graph, "hash") <- get_hash (graph, contracted = FALSE, force = TRUE)

    return (graph)
}

insert_one_edge <- function (graph, index, x, y, gr_cols) {

    if (is.null (x) && is.null (y)) {
        x <- (graph [[gr_cols$xfr]] [index] +
            graph [[gr_cols$xto]] [index]) / 2
        y <- (graph [[gr_cols$yfr]] [index] +
            graph [[gr_cols$yto]] [index]) / 2
    }
    expand_index <- c (1:index, index, (index + 1):nrow (graph))
    graph <- graph [expand_index, ]
    graph [index, gr_cols$xto] <- x
    graph [index, gr_cols$yto] <- y
    graph [index + 1, gr_cols$xfr] <- x
    graph [index + 1, gr_cols$yfr] <- y

    xy1 <- c (
        x = graph [[gr_cols$xfr]] [index],
        y = graph [[gr_cols$yfr]] [index]
    )
    xy2 <- c (
        x = graph [[gr_cols$xto]] [index + 1],
        y = graph [[gr_cols$yto]] [index + 1]
    )
    if (is_graph_spatial (graph)) {
        requireNamespace ("geodist")
        d1 <- geodist::geodist (xy1, c (x = x, y = y), measure = "geodesic")
        d2 <- geodist::geodist (xy2, c (x = x, y = y), measure = "geodesic")
    } else {
        d1 <- sqrt ((xy1 [1] - x)^2 + (xy1 [2] - y)^2)
        d2 <- sqrt ((x - xy2 [1])^2 + (y - xy2 [2])^2)
    }
    wt <- graph [index, gr_cols$d_weighted] /
        graph [index, gr_cols$d]
    if (!is.na (gr_cols$time)) {
        time_scale <- graph [index, gr_cols$time] /
            graph [index, gr_cols$d]
        time_wt <- graph [index, gr_cols$time_weighted] /
            graph [index, gr_cols$time]
    }
    graph [index, gr_cols$d] <- d1
    graph [index, gr_cols$d_weighted] <- d1 * wt
    graph [index + 1, gr_cols$d] <- d2
    graph [index + 1, gr_cols$d_weighted] <- d2 * wt

    if (!is.na (gr_cols$time)) {
        graph [index, gr_cols$time] <- graph [index, gr_cols$d] *
            time_scale
        graph [index, gr_cols$time_weighted] <-
            graph [index, gr_cols$time] * time_wt
    }

    graph$edge_id [index] <- paste0 (graph$edge_id [index], "_a")
    graph$edge_id [index + 1] <- paste0 (graph$edge_id [index + 1], "_b")

    return (graph)
}
