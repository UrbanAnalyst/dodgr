null_to_na <- function (x)
{
    if (length (x) == 0)
        x <- NA
    return (x)
}

#' get_graph_cols
#'
#' Identify the essential columns of the graph table (data.frame, tibble,
#' whatever) to be analysed in the C++ routines.
#'
#' @param graph A `data.frame` containing the edges of the graph
#' @return A list of column numbers of `edge_id`, `from`,
#' `to`, `d`, `w`, `time`, `xfr`, `yfr`, `xto`, `yto`, and `component`, some of
#' which may be NA.
#'
#' @noRd
dodgr_graph_cols <- function (graph)
{
    nms <- names (graph)
    component <- grep ("comp", nms) %>% null_to_na ()
    if (is (graph, "dodgr_streetnet") & 
        !is (graph, "dodgr_streetnet_sc") & ncol (graph) >= 11)
    {
        # columns are always identically structured
        edge_id <- which (nms == "edge_id")
        fr_col <- which (nms == "from_id") %>% null_to_na ()
        to_col <- which (nms == "to_id") %>% null_to_na ()
        d_col <- which (nms == "d")
        w_col <- which (nms == "d_weighted")

        xfr <- which (nms == "from_lon")
        if (length (xfr) == 0) xfr <- NA
        yfr <- which (nms == "from_lat")
        if (length (yfr) == 0) yfr <- NA
        xto <- which (nms == "to_lon")
        if (length (xto) == 0) xto <- NA
        yto <- which (nms == "to_lat")
        if (length (yto) == 0) yto <- NA
    } else
    {
        edge_id <- grep ("edge_id|edge_$", nms) %>% null_to_na ()

        d_col <- find_d_col (graph)
        w_col <- find_w_col (graph)
        if (length (w_col) == 0)
            w_col <- d_col

        fr_col <- find_fr_id_col (graph)
        to_col <- find_to_id_col (graph)

        xfr <- yfr <- xto <- yto <- NA
        # TODO: Modify for other complex but non-spatial types of graph
        if (is_graph_spatial (graph))
        {
            spcols <- find_spatial_cols (graph)
            graph <- tbl_to_df (graph)

            if (!(all (apply (graph [, spcols$fr_col], 2, is.numeric)) |
                  all (apply (graph [, spcols$to_tol], 2, is.numeric))))
                stop (paste0 ("graph appears to have non-numeric ",
                              "longitudes and latitudes"))

            xfr <- spcols$fr_col [1]
            yfr <- spcols$fr_col [2]
            xto <- spcols$to_col [1]
            yto <- spcols$to_col [2]
        } else
        {
            if (length (fr_col) != 1 & length (to_col) != 1)
                stop ("Unable to determine from and to columns in graph")
        }
    }

    time_col <- grep ("time", nms)
    if (length (time_col) != 1)
    {
        time_col <- grep ("time$", nms)
        if (length (time_col) != 1)
            time_col <- NA
    }
    timew_col <- grep ("time_w|timew|tw", nms)
    if (length (timew_col) != 1)
    {
        timew_col <- grep ("time_w|timew|^tw", nms)
        if (length (timew_col) != 1)
            timew_col <- NA
    }

    ret <- c (edge_id, fr_col, to_col, d_col, w_col, time_col, timew_col,
              xfr, yfr, xto, yto, component)
    names (ret) <- c ("edge_id", "from", "to", "d", "w", "time", "time_weighted",
                      "xfr", "yfr", "xto", "yto", "component")
    class (ret) <- c (class (ret), "graph_columns")

    # This is passed to many C++ routines, in which case it needs to be
    # converted to a vector (`do.call (c, gr_cols)`), and the R-style 1-indexeso
    # need to be converted to equivalent 0-indexed forms
    return (as.list (ret))
}

#' convert_graph
#'
#' Convert graph to a standard form suitable for submission to C++ routines
#' @noRd
convert_graph <- function (graph, gr_cols)
{
    keep_cols <- c ("edge_id", "from", "to", "d", "w", "time", "time_weighted")
    index <- do.call (c, gr_cols [keep_cols])
    index <- index [!is.na (index)]
    graph <- graph [, index]
    names (graph) <- names (index)

    if ("edge_id" %in% names (graph))
        graph$edge_id <- convert_to_char (graph$edge_id)
    graph$from <- convert_to_char (graph$from)
    graph$to <- convert_to_char (graph$to)

    if (!"time_weighted" %in% names (graph))
        graph$time_weighted <- graph$time

    return (graph)
}

convert_to_char <- function (x)
{
    if (!is.character (x)) x <- paste0 (x)
    return (x)
}

tbl_to_df <- function (graph)
{
    if (methods::is (graph, "tbl"))
    {
        classes <- class (graph) [!grepl ("tbl", class (graph))]
        graph <- as.data.frame (graph)
        class (graph) <- classes
    }
    return (graph)
}



#' dodgr_vertices
#'
#' Extract vertices of graph, including spatial coordinates if included
#'
#' @param graph A flat table of graph edges. Must contain columns labelled
#' `from` and `to`, or `start` and `stop`. May also contain
#' similarly labelled columns of spatial coordinates (for example
#' `from_x`) or `stop_lon`).
#' @return A `data.frame` of vertices with unique numbers (`n`).
#'
#' @note Values of `n` are 0-indexed
#'
#' @export
#' @examples
#' graph <- weight_streetnet (hampi)
#' v <- dodgr_vertices (graph)
dodgr_vertices <- function (graph)
{
    graph <- tbl_to_df (graph)

    gr_cols <- dodgr_graph_cols (graph)
    nms <- names (gr_cols)
    # cols are (edge_id, from, to, d, w, component, xfr, yfr, xto, yto)
    # NOTE: c (x, y), where x and y are both factors gives junk, so explicit
    # conversion required here: TODO: Find a better way?
    if (is.factor (graph [[gr_cols$from]]))
        graph [[gr_cols$from]] <- paste0 (graph [[gr_cols$from]])
    if (is.factor (graph [[gr_cols$to]]))
        graph [[gr_cols$to]] <- paste0 (graph [[gr_cols$to]])
    if (is_graph_spatial (graph))
    {
        verts <- data.frame (id = c (graph [[gr_cols$from]],
                                     graph [[gr_cols$to]]),
                             x = c (graph [[gr_cols$xfr]],
                                    graph [[gr_cols$xto]]),
                             y = c (graph [[gr_cols$yfr]],
                                    graph [[gr_cols$yto]]),
                             stringsAsFactors = FALSE)
        if (!is.na (gr_cols$component))
            verts$component <- graph [[gr_cols$component]]
    } else
    {
        verts <- data.frame (id = c (graph [[gr_cols$from]],
                                     graph [[gr_cols$to]]),
                             stringsAsFactors = FALSE)
        if (!is.na (gr_cols$component))
            verts$component <- graph [[gr_cols$component]]
    }

    indx <- which (!duplicated (verts))
    verts <- verts [indx, , drop = FALSE] #nolint
    verts$n <- seq (nrow (verts)) - 1

    return (verts)
}


#' dodgr_components
#'
#' Identify connected components of graph and add corresponding `component`
#' column to `data.frame`.
#'
#' @param graph A `data.frame` of edges
#' @return Equivalent graph with additional `component` column,
#' sequentially numbered from 1 = largest component.
#' @export
#' @examples
#' graph <- weight_streetnet (hampi)
#' graph <- dodgr_components (graph)
dodgr_components <- function (graph)
{
    graph <- tbl_to_df (graph)

    if ("component" %in% names (graph))
        message ("graph already has a component column")
    else
    {
        gr_cols <- dodgr_graph_cols (graph)
        graph2 <- convert_graph (graph, gr_cols)
        if (is.na (gr_cols$edge_id))
            graph2$edge_id <- seq (nrow (graph2))
        cns <- rcpp_get_component_vector (graph2)

        indx <- match (graph2$edge_id, cns$edge_id)
        component <- cns$edge_component [indx]
        # Then re-number in order to decreasing component size:
        graph$component <- match (component,
                                  order (table (component), decreasing = TRUE))
    }

    return (graph)
}

#' dodgr_sample
#'
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
#' @note Graphs may occassionally have `nverts + 1` vertices, rather than
#' the requested `nverts`.
#'
#' @export
#' @examples
#' graph <- weight_streetnet (hampi)
#' nrow (graph) # 5,742
#' graph <- dodgr_sample (graph, nverts = 200)
#' nrow (graph) # generally around 400 edges
#' nrow (dodgr_vertices (graph)) # 200
dodgr_sample <- function (graph, nverts = 1000)
{
    classes <- class (graph)
    graph <- tbl_to_df (graph)

    fr <- find_fr_id_col (graph)
    to <- find_to_id_col (graph)
    verts <- unique (c (graph [, fr], graph [, to]))
    if (length (verts) > nverts)
    {
        gr_cols <- dodgr_graph_cols (graph)
        graph2 <- convert_graph (graph, gr_cols)
        indx <- match (rcpp_sample_graph (graph2, nverts), graph2$edge_id)
        graph <- graph [sort (indx), ]
    }

    class (graph) <- classes
    return (graph)
}
