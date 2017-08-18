#' is_graph_spatial
#'
#' Is the graph spatial or not?
#' @param graph A \code{data.frame} of edges
#' @return \code{TRUE} is \code{graph} is spatial, otherwise \code{FALSE}
#' @noRd
is_graph_spatial <- function (graph)
{
    ncol (graph) > 4 &
        (any (grepl ("x", names (graph), ignore.case = TRUE)) |
         any (grepl ("y", names (graph), ignore.case = TRUE)) |
         any (grepl ("lon", names (graph), ignore.case = TRUE)) |
         any (grepl ("lat", names (graph), ignore.case = TRUE)))
}


#' convert_graph
#'
#' Convert graph to standard 4-column format for submission to C++ routines
#' @noRd
convert_graph <- function (graph)
{
    if (any (grepl ("edge", names (graph))))
        edge_id <- graph [, which (grepl ("edge", names (graph)))]
    else
        edge_id <- seq (nrow (graph))

    graph$edge_id <- seq (nrow (graph))

    d_col <- which (tolower (substring (names (graph), 1, 1)) == "d" &
                    tolower (substring (names (graph), 2, 2)) != "w" &
                    tolower (substring (names (graph), 2, 2)) != "_")
    w_col <- which (tolower (substring (names (graph), 1, 1)) == "w" |
                    tolower (substring (names (graph), 1, 2)) == "dw" |
                    tolower (substring (names (graph), 1, 3)) == "d_w")
    if (length (d_col) > 1 | length (w_col) > 1)
        stop ("Unable to determine distance and/or weight columns in graph")
    else if (length (d_col) != 1)
        stop ("Unable to determine distance column in graph")

    if (length (w_col) == 0)
        w_col <- d_col

    fr_col <- which (grepl ("fr", names (graph), ignore.case = TRUE))
    to_col <- which (grepl ("to", names (graph), ignore.case = TRUE))

    xy <- NULL
    if (ncol (graph) > 4)
    {
        if (is_graph_spatial (graph))
        {
            if (length (fr_col) != length (to_col))
                stop (paste0 ("from and to columns in graph appear ",
                              "to have different strutures"))
            else if (length (fr_col) >= 2 & length (to_col) >= 2)
            {
                if (length (fr_col) == 3)
                {
                    frx_col <- find_xy_col (graph, fr_col, x = TRUE)
                    fry_col <- find_xy_col (graph, fr_col, x = FALSE)
                    frid_col <- fr_col [which (!fr_col %in%
                                               c (frx_col, fry_col))]
                    fr_col <- c (frx_col, fry_col)
                    xy_fr_id <- graph [, frid_col]
                    if (!is.character (xy_fr_id))
                        xy_fr_id <- paste0 (xy_fr_id)

                    tox_col <- find_xy_col (graph, to_col, x = TRUE)
                    toy_col <- find_xy_col (graph, to_col, x = FALSE)
                    toid_col <- to_col [which (!to_col %in%
                                               c (tox_col, toy_col))]
                    to_col <- c (tox_col, toy_col)
                    xy_to_id <- graph [, toid_col]
                    if (!is.character (xy_to_id))
                        xy_to_id <- paste0 (xy_to_id)
                } else # len == 2, so must be only x-y
                {
                    xy_fr_id <- paste0 (graph [, fr_col [1]],
                                        graph [, fr_col [2]])
                    xy_to_id <- paste0 (graph [, to_col [1]],
                                        graph [, to_col [2]])
                }

                xy_fr <- graph [, fr_col]
                xy_to <- graph [, to_col]
                if (!(all (apply (xy_fr, 2, is.numeric)) |
                      all (apply (xy_to, 2, is.numeric))))
                    stop (paste0 ("graph appears to have non-numeric ",
                                  "longitudes and latitudes"))

                # This same indx is created in vert_map to ensure it follows
                # same order as xy
                indx <- which (!duplicated (c (xy_fr_id, xy_to_id)))
                xy <- data.frame ("x" = c (graph [, fr_col [1]],
                                           graph [, to_col [1]]),
                                  "y" = c (graph [, fr_col [2]],
                                           graph [, to_col [2]])) [indx, ]

                # then replace 4 xy from/to cols with 2 from/to cols
                graph <- data.frame ("edge_id" = edge_id,
                                     "from" = xy_fr_id,
                                     "to" = xy_to_id,
                                     "d" = graph [, d_col],
                                     "w" = graph [, w_col],
                                     stringsAsFactors = FALSE)
            }
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
    } else if (ncol (graph == 4))
    {
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

    return (list (graph = graph, xy = xy))
}

#' find_xy_col
#' @noRd
find_xy_col <- function (graph, indx, x = TRUE)
{
    if (x)
        coli <- which (grepl ("x", names (graph) [indx], ignore.case = TRUE) |
                       grepl ("lon", names (graph) [indx], ignore.case = TRUE))
    else
        coli <- which (grepl ("y", names (graph) [indx], ignore.case = TRUE) |
                       grepl ("lat", names (graph) [indx], ignore.case = TRUE))

    indx [coli]
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
dodgr_vertices <- function (graph)
{
    fr_col <- which (grepl ("fr", names (graph), ignore.case = TRUE) |
                     grepl ("sta", names (graph), ignore.case = TRUE))
    to_col <- which (grepl ("to", names (graph), ignore.case = TRUE) |
                     grepl ("sto", names (graph), ignore.case = TRUE))

    if (is_graph_spatial (graph))
    {
        # graph is spatial
        if (length (fr_col) != length (to_col))
            stop (paste0 ("from and to columns in graph appear ",
                          "to have different strutures"))
        if (length (fr_col) < 2 | length (to_col) < 2)
            stop (paste0 ("Graph appears to be spatial yet unable to ",
                          "extract coordinates."))

        if (length (fr_col) == 3)
        {
            frx_col <- find_xy_col (graph, fr_col, x = TRUE)
            fry_col <- find_xy_col (graph, fr_col, x = FALSE)
            frid_col <- fr_col [which (!fr_col %in%
                                       c (frx_col, fry_col))]
            fr_col <- c (frx_col, fry_col)
            xy_fr_id <- graph [, frid_col]
            if (!is.character (xy_fr_id))
                xy_fr_id <- paste0 (xy_fr_id)

            tox_col <- find_xy_col (graph, to_col, x = TRUE)
            toy_col <- find_xy_col (graph, to_col, x = FALSE)
            toid_col <- to_col [which (!to_col %in%
                                       c (tox_col, toy_col))]
            to_col <- c (tox_col, toy_col)
            xy_to_id <- graph [, toid_col]
            if (!is.character (xy_to_id))
                xy_to_id <- paste0 (xy_to_id)
        } else # len == 2, so must be only x-y
        {
            xy_fr_id <- paste0 (graph [, fr_col [1]],
                                graph [, fr_col [2]])
            xy_to_id <- paste0 (graph [, to_col [1]],
                                graph [, to_col [2]])
        }

        xy_fr <- graph [, fr_col]
        xy_to <- graph [, to_col]
        if (!(all (apply (xy_fr, 2, is.numeric)) |
              all (apply (xy_to, 2, is.numeric))))
            stop (paste0 ("graph appears to have non-numeric ",
                          "longitudes and latitudes"))

        verts <- data.frame (id = c (xy_fr_id, xy_to_id),
                             x = c (graph [, fr_col [1]],
                                    graph [, to_col [1]]),
                             y = c (graph [, fr_col [2]],
                                    graph [, to_col [2]]),
                             stringsAsFactors = FALSE)
    } else # non-spatial graph
    {
        if (!(length (fr_col) == 1 & length (to_col) == 1))
            stop ("Graph appears to be non-spatial, yet unable to ",
                  "determine vertex columns")

        verts <- data.frame (from = graph [, fr_col],
                             to = graph [, to_col])
    }
    indx <- which (!duplicated (verts))
    verts <- verts [indx, ]
    verts$n <- seq (nrow (verts)) - 1

    return (verts)
}

#' match_pts_to_graph
#'
#' Match spatial points to a spatial graph which contains vertex coordindates
#'
#' @param verts A \code{data.frame} of vertices obtained from
#' \code{dodgr_vertices(graph)}.
#' @param xy coordinates of points to be matched to the vertices
#'
#' @return A vector index into verts
#' @noRd
match_pts_to_graph <- function (verts, xy)
{
    if (!(is.matrix (xy) | is.data.frame (xy)))
        stop ("xy must be a matrix or data.frame")
    if (ncol (xy) != 2)
        stop ("xy must have only two columns")

    nms <- names (verts)
    if (is.null (nms))
        nms <- colnames (verts)
    if (!is.null (nms))
    {
        ix <- which (grepl ("x", nms, ignore.case = TRUE) |
                     grepl ("lon", nms, ignore.case = TRUE))
        iy <- which (grepl ("y", nms, ignore.case = TRUE) |
                     grepl ("lat", nms, ignore.case = TRUE))
    } else
    {
        message ("xy has no named columns; assuming order is x then y")
        ix <- 1
        iy <- 2
    }

    verts <- data.frame (x = verts [, ix], y = verts [, iy])

    rcpp_points_index (verts, xy)
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
dodgr_components <- function (graph)
{
    if ("component" %in% names (graph))
        message ("graph already has a component column")
    else
    {
        graphc <- convert_graph (graph)$graph
        if (!(any (grepl ("wt", names (graph))) |
              any (grepl ("weight", names (graph)))))
        cn <- rcpp_get_component_vector (graphc)
        # Then re-number in order to decreasing component size:
        cn <- match (cn, order (table (cn), decreasing = TRUE))
        graph$component <- cn
    }

    return (graph)
}

#' dodgr_contract_graph
#'
#' Removes redundant (straight-line) vertices from graph, leaving only junction
#' vertices.
#'
#' @param graph A flat table of graph edges. Must contain columns labelled
#' \code{from} and \code{to}, or \code{start} and \code{stop}. May also contain
#' similarly labelled columns of spatial coordinates (for example
#' \code{from_x}) or \code{stop_lon}).
#' @param verts Option \code{data.frame} of vertices obtained from
#' \code{dodgr_vertices} (submitting this will simply speed up conversion to
#' compact graph).
#' @param quiet If \code{FALSE}, display progress on screen
#'
#' @return A complex object with both the original graph and its compact verion
#' (\code{$original} and \code{$compact}, respectively), along with several
#' indexes used to map vertices and edges between the two.
#' @export
dodgr_contract_graph <- function (graph, verts = NULL, quiet = TRUE)
{
    if (is.null (verts))
        verts <- dodgr_vertices (graph)

    rcpp_contract_graph (graph, quiet = quiet)
}

#' dodgr_sample_graph
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
dodgr_sample <- function (graph, nverts = 1000)
{
    verts <- unique (c (graph$from_id, graph$to_id))
    if (length (verts) > nverts)
    {
        indx <- rcpp_sample_graph (convert_graph (graph)$graph, nverts)
        graph <- graph [sort (indx), ]
    }

    return (graph)
}
