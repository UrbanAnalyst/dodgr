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
#' @param components If FALSE, components are not calculated (will generally
#' result in faster processing).
#' @return A named vector of column numbers of `edge_id`, `from`,
#' `to`, `d`, `w`, `xfr`, `yfr`, `xto`,
#' `yto`, and `component`, some of which may be NA.
#'
#' @noRd
dodgr_graph_cols <- function (graph)
{
    nms <- names (graph)
    component <- which (grepl ("comp", nms)) %>% null_to_na ()
    if (is (graph, "dodgr_streetnet") & ncol (graph) >= 11)
    {
        # columns are always identically structured
        edge_id <- which (nms == "edge_id")
        fr_col <- which (nms == "from_id") %>% null_to_na ()
        to_col <- which (nms == "to_id") %>% null_to_na ()
        d_col <- which (nms == "d")
        w_col <- which (nms == "d_weighted")
        xfr <- which (nms == "from_lon")
        yfr <- which (nms == "from_lat")
        xto <- which (nms == "to_lon")
        yto <- which (nms == "to_lat")
    } else
    {
        edge_id <- which (grepl ("edge_id", nms)) %>% null_to_na ()

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

    # This is NOT a list because it's much easier to pass as vector to C++
    ret <- c (edge_id, fr_col, to_col, d_col, w_col, xfr, yfr, xto, yto)
    names (ret) <- c ("edge_id", "from", "to", "d", "w",
                      "xfr", "yfr", "xto", "yto")
    if (!is.na (component))
    {
        ret <- c (ret, component)
        names (ret) [length (ret)] <- "component"
    }
    class (ret) <- c (class (ret), "graph_columns")

    # Note that these are R-style 1-indexed, so need to be converted in C++ to
    # equivalent 0-indexed forms
    return (ret)
}

#' convert_graph
#'
#' Convert graph to a standard form suitable for submission to C++ routines
#' @noRd
convert_graph <- function (graph, gr_cols)
{
    indx <- which (!is.na (gr_cols [1:5]))
    graph <- graph [, gr_cols [1:5] [indx]]
    names (graph) <- c ("edge_id", "from", "to", "d", "w") [indx]
    if ("edge_id" %in% names (graph))
        if (!is.character (graph$edge_id))
            graph$edge_id <- paste0 (graph$edge_id)
    if (!is.character (graph$from))
        graph$from <- paste0 (graph$from)
    if (!is.character (graph$to))
        graph$to <- paste0 (graph$to)
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
    cols <- dodgr_graph_cols (graph)
    nms <- names (cols)
    # cols are (edge_id, from, to, d, w, component, xfr, yfr, xto, yto)
    # NOTE: c (x, y), where x and y are both factors gives junk, so explicit
    # conversion required here: TODO: Find a better way?
    from_id <- graph [[cols [which (nms == "from")] ]]
    to_id <- graph [[cols [which (nms == "to")] ]]
    if (is.factor (from_id))
        from_id <- paste0 (from_id)
    if (is.factor (to_id))
        to_id <- paste0 (to_id)
    if (is_graph_spatial (graph))
    {
        verts <- data.frame (id = c (from_id, to_id),
                             x = c (graph [[cols [which (nms == "xfr")] ]],
                                    graph [[cols [which (nms == "xto")] ]]),
                             y = c (graph [[cols [which (nms == "yfr")] ]],
                                    graph [[cols [which (nms == "yto")] ]]),
                             stringsAsFactors = FALSE)
        if ("component" %in% nms)
            verts$component <- graph [[cols [which (nms == "component")] ]]
    } else
    {
        verts <- data.frame (id = c (from_id, to_id),
                             stringsAsFactors = FALSE)
        if ("component" %in% nms)
            verts <- cbind (verts, rep (graph$component, 2))
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
    if ("component" %in% names (graph))
        message ("graph already has a component column")
    else
    {
        gr_cols <- dodgr_graph_cols (graph)
        graph2 <- convert_graph (graph, gr_cols)
        if (is.na (gr_cols [which (names (gr_cols) == "edge_id")]))
            graph2$edge_id <- seq (nrow (graph2))
        cns <- rcpp_get_component_vector (graph2)

        indx <- match (graph2$edge_id, cns$edge_id)
        component <- cns$edge_component [indx]
        # Then re-number in order to decreasing component size:
        graph$component <- match (component, order (table (component),
                                                    decreasing = TRUE))
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
    fr <- find_fr_id_col (graph)
    to <- find_to_id_col (graph)
    verts <- unique (c (graph [, fr], graph [, to]))
    if (length (verts) > nverts)
    {
        gr_cols <- dodgr_graph_cols (graph)
        graph2 <- convert_graph (graph, gr_cols)
        indx <- match (rcpp_sample_graph (graph2, nverts), graph$edge_id)
        graph <- graph [sort (indx), ]
    }

    return (graph)
}
