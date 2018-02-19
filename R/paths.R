#' dodgr_paths
#'
#' Calculate lists of pair-wise shortest paths between points.
#'
#' @param graph \code{data.frame} or equivalent object representing the network
#' graph (see Details)
#' @param from Vector or matrix of points **from** which route distances are to
#' be calculated (see Details)
#' @param to Vector or matrix of points **to** which route distances are to be
#' calculated (see Details)
#' @param vertices If \code{TRUE}, return lists of lists of vertices for each
#' path, otherwise return corresponding lists of edge numbers from \code{graph}.
#' @param wt_profile Name of weighting profile for street networks (one of foot,
#' horse, wheelchair, bicycle, moped, motorcycle, motorcar, goods, hgv, psv).
#' @param heap Type of heap to use in priority queue. Options include
#' Fibonacci Heap (default; \code{FHeap}), Binary Heap (\code{BHeap}),
#' \code{Radix}, Trinomial Heap (\code{TriHeap}), Extended Trinomial Heap
#' (\code{TriHeapExt}, and 2-3 Heap (\code{Heap23}).
#' @param quiet If \code{FALSE}, display progress messages on screen.
#' @return List of list of paths tracing all connections between nodes such that
#' if \code{x <- dodgr_paths (graph, from, to)}, then the path between
#' \code{from[i]} and \code{to[j]} is \code{x [[i]] [[j]]}.
#'
#' @note \code{graph} must minimally contain four columns of \code{from},
#' \code{to}, \code{dist}. If an additional column named \code{weight} or
#' \code{wt} is present, shortest paths are calculated according to values
#' specified in that column; otherwise according to \code{dist} values. Either
#' way, final distances between \code{from} and \code{to} points are calculated
#' according to values of \code{dist}. That is, paths between any pair of points
#' will be calculated according to the minimal total sum of \code{weight}
#' values (if present), while reported distances will be total sums of
#' \code{dist} values.
#'
#' The \code{from} and \code{to} columns of \code{graph} may be either single
#' columns of numeric or character values specifying the numbers or names of
#' graph vertices, or combinations to two columns specifying geographical
#' (longitude and latitude) coordinates. In the latter case, almost any sensible
#' combination of names will be accepted (for example, \code{fromx, fromy},
#' \code{from_x, from_y}, or \code{fr_lat, fr_lon}.)
#'
#' \code{from} and \code{to} values can be either two-column matrices of
#' equivalent of longitude and latitude coordinates, or else single columns
#' precisely matching node numbers or names given in \code{graph$from} or
#' \code{graph$to}. If \code{to} is missing, pairwise distances are calculated
#' between all points specified in \code{from}. If neither \code{from} nor
#' \code{to} are specified, pairwise distances are calculated between all nodes
#' in \code{graph}.
#'
#' @export
#' @examples
#' graph <- weight_streetnet (hampi)
#' from <- sample (graph$from_id, size = 100)
#' to <- sample (graph$to_id, size = 50)
#' dp <- dodgr_paths (graph, from = from, to = to)
#' # dp is a list with 100 items, and each of those 100 items has 30 items, each
#' # of which is a single path listing all vertiex IDs as taken from \code{graph}.
dodgr_paths <- function (graph, from, to, vertices = TRUE,
                         wt_profile = "bicycle", heap = 'BHeap', quiet = TRUE)
{
    if (missing (graph) & (!missing (from) | !missing (to)))
        graph <- graph_from_pts (from, to, expand = 0.1,
                                 wt_profile = wt_profile, quiet = quiet)

    hps <- get_heap (heap, graph)
    heap <- hps$heap
    graph <- hps$graph

    gr_cols <- dodgr_graph_cols (graph)
    # cols are (edge_id, from, to, d, w, component, xfr, yfr, xto, yto)
    vert_map <- make_vert_map (graph, gr_cols)

    index_id <- get_index_id_cols (graph, gr_cols, vert_map, from)
    from_index <- index_id$index - 1 # 0-based
    if (!is.null (index_id$id))
        from_id <- index_id$id # can be null
    else
        from_id <- vert_map$vert

    index_id <- get_index_id_cols (graph, gr_cols, vert_map, to)
    to_index <- index_id$index - 1 # 0-based
    if (!is.null (index_id$id))
        to_id <- index_id$id # can be null
    else
        to_id <- vert_map$vert

    graph <- convert_graph (graph, gr_cols)

    if (!quiet)
        message ("Calculating shortest paths ... ", appendLF = FALSE)
    paths <- rcpp_get_paths (graph, vert_map, from_index, to_index, heap)

    # convert 1-based indices back into vertex IDs:
    paths <- lapply (paths, function (i)
                     lapply (i, function (j)
                             vert_map$vert [j] ))
    # name path lists
    if (!is.null (from_id) & !is.null (to_id))
    {
        for (i in seq (from_id))
            names (paths [[i]]) <- paste0 (from_id [i], "-", to_id)
        names (paths) <- from_id
    }

    if (!vertices)
    {
        # convert vertex IDs to corresponding sequences of edge numbers
        graph_verts <- paste0 ("f", graph$from, "t", graph$to)

        paths <- lapply (paths, function (i)
                         lapply (i, function (j)
                                 if (length (j) > 0)
                                 {
                                     indx <- 2:length (j)
                                     pij <- paste0 ("f", j [indx - 1],
                                                    "t", j [indx])
                                     match (pij, graph_verts)
                                 } ))
    }

    return (paths)
}
