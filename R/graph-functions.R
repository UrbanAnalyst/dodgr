#' dodgr_convert_graph
#'
#' Convert a graph represented as an arbitarily-structured \code{data.frame} to
#' standard 4 or 5-column format for submission to C++ routines
#'
#' @param graph A \code{data.frame} containing the edges of the graph
#' @param components If FALSE, components are not calculated (will generally
#' result in faster processing).
#' @return A list of two components: (i) \code{graph}, a \code{data.frame} with
#' the same number of rows but with columns of \code{edge_id}, \code{from},
#' \code{to}, \code{d}, \code{w}, and \code{component}, not all of which may be
#' included; and (ii) \code{xy}, a matrix of coordinates of all vertices in
#' \code{graph}.
#'
#' @export
#' @examples
#' graph <- weight_streetnet (hampi)
#' names (graph)
#' names (dodgr_convert_graph (graph)$graph)
dodgr_convert_graph <- function (graph, components = TRUE)
{
    if (is (graph, "graph_converted"))
        return (list (graph = graph))

    if (any (grepl ("edge", names (graph))))
        edge_id <- graph [, which (grepl ("edge", names (graph)))]
    else
        edge_id <- seq (nrow (graph))

    component <- NULL
    if (any (grepl ("comp", names (graph))))
    {
        component <- graph [, which (grepl ("comp", names (graph)))]
        components <- TRUE
    }

    d_col <- find_d_col (graph)
    w_col <- find_w_col (graph)
    if (length (w_col) == 0)
        w_col <- d_col

    fr_col <- find_fr_col (graph)
    to_col <- find_to_col (graph)
    if (length (fr_col) != length (to_col))
        stop (paste0 ("from and to columns in graph appear ",
                      "to have different strutures"))

    xy <- NULL
    # TODO: Modify for other complex but non-spatial types of graph
    if (is_graph_spatial (graph))
    {
        spcols <- find_spatial_cols (graph)

        xy_fr <- graph [, spcols$fr_col]
        xy_to <- graph [, spcols$to_col]
        if (!(all (apply (xy_fr, 2, is.numeric)) |
              all (apply (xy_to, 2, is.numeric))))
            stop (paste0 ("graph appears to have non-numeric ",
                          "longitudes and latitudes"))

        # This same indx is created in vert_map to ensure it follows
        # same order as xy
        indx <- which (!duplicated (c (spcols$xy_id$xy_fr_id,
                                       spcols$xy_id$xy_to_id)))
        xy <- data.frame ("x" = c (graph [, spcols$fr_col [1]],
                                   graph [, spcols$to_col [1]]),
                          "y" = c (graph [, spcols$fr_col [2]],
                                   graph [, spcols$to_col [2]])) [indx, ]

        # then replace 4 xy from/to cols with 2 from/to cols
        graph <- data.frame ("edge_id" = edge_id,
                             "from" = spcols$xy_id$xy_fr_id,
                             "to" = spcols$xy_id$xy_to_id,
                             "d" = graph [, d_col],
                             "w" = graph [, w_col],
                             stringsAsFactors = FALSE)
    } else
    {
        if (length (fr_col) != 1 & length (to_col) != 1)
            stop ("Unable to determine from and to columns in graph")

        graph <- data.frame ("edge_id" = edge_id,
                             "from" = graph [, fr_col],
                             "to" = graph [, to_col],
                             "d" = graph [, d_col],
                             "w" = graph [, w_col],
                             stringsAsFactors = FALSE)
    }

    if (!is.character (graph$from))
        graph$from <- paste0 (graph$from)
    if (!is.character (graph$to))
        graph$to <- paste0 (graph$to)
    if (!is.character (graph$edge_id))
        graph$edge_id <- paste0 (graph$edge_id)

    if (components)
    {
        if (is.null (component))
        {
            cns <- rcpp_get_component_vector (graph)
            component <- cns$edge_component [match (paste0 (graph$edge_id),
                                                    cns$edge_id)]
            # Then re-number in order to decreasing component size:
            component <- match (component, order (table (component),
                                                  decreasing = TRUE))
        }
        graph$component <- component
    }

    class (graph) <- c (class (graph), "graph_converted")

    return (list (graph = graph, xy = xy))
}


#' dodgr_vertices
#'
#' Extract vertices of graph, including spatial coordinates if included
#'
#' @param graph A flat table of graph edges. Must contain columns labelled
#' \code{from} and \code{to}, or \code{start} and \code{stop}. May also contain
#' similarly labelled columns of spatial coordinates (for example
#' \code{from_x}) or \code{stop_lon}).
#' @return A \code{data.frame} of vertices with unique numbers (\code{n}).
#'
#' @note Values of \code{n} are 0-indexed
#'
#' @export
#' @examples
#' graph <- weight_streetnet (hampi)
#' v <- dodgr_vertices (graph)
dodgr_vertices <- function (graph)
{
    fr_col <- find_fr_col (graph)
    to_col <- find_to_col (graph)

    if (is_graph_spatial (graph))
    {
        spcols <- find_spatial_cols (graph)

        xy_fr <- graph [, spcols$fr_col]
        xy_to <- graph [, spcols$to_col]
        if (!(all (apply (xy_fr, 2, is.numeric)) |
              all (apply (xy_to, 2, is.numeric))))
            stop (paste0 ("graph appears to have non-numeric ",
                          "longitudes and latitudes"))

        verts <- data.frame (id = c (spcols$xy_id$xy_fr_id,
                                     spcols$xy_id$xy_to_id),
                             x = c (xy_fr [, 1], xy_to [, 1]),
                             y = c (xy_fr [, 2], xy_to [, 2]),
                             stringsAsFactors = FALSE)
    } else # non-spatial graph
    {
        if (!(length (fr_col) == 1 & length (to_col) == 1))
            stop ("Graph appears to be non-spatial, yet unable to ",
                  "determine vertex columns")

        verts <- data.frame (id = c (graph [, fr_col], graph [, to_col]),
                             stringsAsFactors = FALSE)
    }
    indx <- which (!duplicated (verts))
    verts <- verts [indx, , drop = FALSE] #nolint
    verts$n <- seq (nrow (verts)) - 1

    return (verts)
}


#' dodgr_components
#'
#' Identify connected components of graph and add corresponding \code{component}
#' column to \code{data.frame}.
#'
#' @param graph A \code{data.frame} of edges
#' @return Equivalent graph with additional \code{component} column,
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
        component <- dodgr_convert_graph (graph, components = TRUE)
        graph$component <- component$graph$component
    }

    return (graph)
}

#' dodgr_sample
#'
#' Sample a random but connected sub-component of a graph
#'
#' @param graph A flat table of graph edges. Must contain columns labelled
#' \code{from} and \code{to}, or \code{start} and \code{stop}. May also contain
#' similarly labelled columns of spatial coordinates (for example
#' \code{from_x}) or \code{stop_lon}).
#' @param nverts Number of vertices to sample
#'
#' @return A connected sub-component of \code{graph}
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
        grc <- dodgr_convert_graph (graph)$graph
        indx <- match (rcpp_sample_graph (grc, nverts), grc$edge_id)
        graph <- graph [sort (indx), ]
    }

    return (graph)
}
