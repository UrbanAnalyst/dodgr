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

#' Get graph columns containing the from vertex
#' @noRd
find_fr_col <- function (graph)
{
    which (grepl ("fr", names (graph), ignore.case = TRUE) |
           grepl ("sta", names (graph), ignore.case = TRUE))
}

#' Get graph columns containing the to vertex
#' @noRd
find_to_col <- function (graph)
{
    which (grepl ("to", names (graph), ignore.case = TRUE) |
           grepl ("sto", names (graph), ignore.case = TRUE))
}

#' Get single graph column containing the ID of the from vertex
#' @noRd
find_fr_id_col <- function (graph)
{
    fr_col <- find_fr_col (graph)
    if (is_graph_spatial (graph))
    {
        frx_col <- find_xy_col (graph, fr_col, x = TRUE)
        fry_col <- find_xy_col (graph, fr_col, x = FALSE)
        fr_col <- fr_col [which (!fr_col %in%
                                 c (frx_col, fry_col))]
    }
    if (length (fr_col) != 1)
        stop ("Unable to determine column with ID of from vertices")
    return (fr_col)
}

#' Get single graph column containing the ID of the to vertex
#' @noRd
find_to_id_col <- function (graph)
{
    to_col <- find_to_col (graph)
    if (is_graph_spatial (graph))
    {
        tox_col <- find_xy_col (graph, to_col, x = TRUE)
        toy_col <- find_xy_col (graph, to_col, x = FALSE)
        to_col <- to_col [which (!to_col %in%
                                 c (tox_col, toy_col))]
    }
    if (length (to_col) != 1)
        stop ("Unable to determine column with ID of from vertices")
    return (to_col)
}

#' find_xy_col
#'
#' Find columns in graph containing lon and lat coordinates
#' @param indx columns of graph containing either to or from values, so xy
#' columns can be returned separately for each case
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

#' find_spatial_cols
#' @noRd
find_spatial_cols <- function (graph)
{

    fr_col <- find_fr_col (graph)
    to_col <- find_to_col (graph)

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

    list (fr_col = fr_col,
          to_col = to_col,
          xy_id = data.frame (xy_fr_id = xy_fr_id,
                              xy_to_id = xy_to_id,
                              stringsAsFactors = FALSE))
}

find_d_col <- function (graph)
{
    d_col <- which (tolower (substring (names (graph), 1, 1)) == "d" &
                    tolower (substring (names (graph), 2, 2)) != "w" &
                    tolower (substring (names (graph), 2, 2)) != "_")
    if (length (d_col) != 1)
        stop ("Unable to determine distance column in graph")
    return (d_col)
}

find_w_col <- function (graph)
{
    w_col <- which (tolower (substring (names (graph), 1, 1)) == "w" |
                    tolower (substring (names (graph), 1, 2)) == "dw" |
                    tolower (substring (names (graph), 1, 3)) == "d_w")
    if (length (w_col) > 1)
        stop ("Unable to determine weight column in graph")
    return (w_col)
}

#' convert_graph
#'
#' Convert a graph represented as an arbitarily-structured \code{data.frame} to
#' standard 4 or 5-column format for submission to C++ routines
#'
#' @param graph A \code{data.frame} containing the edges of the graph
#' @param components If FALSE, components are not calculated (will generally
#' result in faster processing).
#' @return A list of two components: \code{graph}, a \code{data.frame} with the
#' same number of rows but with columns of \code{edge_id}, \code{from},
#' \code{to}, \code{d}, \code{w}, and \code{component}, not all of which may be
#' included; and \code{xy}, A matrix of coordinates of all vertices in
#' \code{graph}.
#'
#' @export
#' @examples
#' graph <- weight_streetnet (hampi)
#' names (graph)
#' names (convert_graph (graph)$graph)
convert_graph <- function (graph, components = TRUE)
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
#' @examples
#' graph <- weight_streetnet (hampi)
#' graph <- dodgr_components (graph)
dodgr_components <- function (graph)
{
    if ("component" %in% names (graph))
        message ("graph already has a component column")
    else
        graph$component <- convert_graph (graph)$graph$component

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
#' @return A list with the contracted graph (\code{graph}), and a map
#' (\code{map}) between the IDs of edges in the contracted graph and those in
#' the original, full graph.
#' @export
#' @examples
#' graph <- weight_streetnet (hampi)
#' nrow (graph) # 5,742
#' graph <- dodgr_contract_graph (graph)
#' nrow (graph$graph) # size of contracted graph = 2,878
#' nrow (graph$map) # 3,692 = number of old edges replaced by new ones
dodgr_contract_graph <- function (graph, verts = NULL, quiet = TRUE)
{
    graph_converted <- convert_graph (graph)$graph
    graph_contracted <- rcpp_contract_graph (graph_converted, quiet = quiet)

    # re-insert spatial coordinates:
    if (is_graph_spatial (graph))
    {
        spcols <- find_spatial_cols (graph)
        fr_col <- find_fr_id_col (graph)
        indx <- match (graph_contracted$graph$from, graph [, fr_col])
        fr_xy <- graph [indx, spcols$fr_col]
        to_col <- find_to_id_col (graph)
        indx <- match (graph_contracted$graph$to, graph [, to_col])
        to_xy <- graph [indx, spcols$to_col]
        graph_contracted$graph <- cbind (graph_contracted$graph,
                                         fr_xy, to_xy)
    }

    return (graph_contracted)
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
        grc <- convert_graph (graph)$graph
        indx <- match (rcpp_sample_graph (grc, nverts), grc$edge_id)
        graph <- graph [sort (indx), ]
    }

    return (graph)
}
