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

    if (!quiet)
        message ("Converting network to dodgr graph ... ",
                 appendLF = FALSE)
    graph <- dodgr_convert_graph (graph, components = FALSE)
    xy <- graph$xy
    graph <- graph$graph
    vert_map <- make_vert_map (graph)
    # vert_map$vert is char vertex ID; vert_map$id is 0-indexed integer

    from_index <- get_tofrom_index (vert_map, xy, from) # 0-indexed
    to_index <- get_tofrom_index (vert_map, xy, to)

    if (!quiet)
        message ("done\nCalculating shortest paths ... ", appendLF = FALSE)
    paths <- rcpp_get_paths (graph, vert_map, from_index, to_index, heap)

    # convert 1-based indices back into vertex IDs:
    paths <- lapply (paths, function (i)
                     lapply (i, function (j)
                             vert_map$vert [j] ))

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



#' dodgr_flows
#'
#' Aggregate flows throughout a network based on an input matrix of flows
#' between all pairs of \code{from} and \code{to} points.
#'
#' @param graph \code{data.frame} or equivalent object representing the network
#' graph (see Details)
#' @param from Vector or matrix of points **from** which route distances are to
#' be calculated (see Details)
#' @param to Vector or matrix of points **to** which route distances are to be
#' calculated (see Details)
#' @param flows Matrix of flows with \code{nrow(flows)==length(from)} and
#' \code{ncol(flows)==length(to)}.
#' @param wt_profile Name of weighting profile for street networks (one of foot,
#' horse, wheelchair, bicycle, moped, motorcycle, motorcar, goods, hgv, psv).
#' @param heap Type of heap to use in priority queue. Options include
#' Fibonacci Heap (default; \code{FHeap}), Binary Heap (\code{BHeap}),
#' \code{Radix}, Trinomial Heap (\code{TriHeap}), Extended Trinomial Heap
#' (\code{TriHeapExt}, and 2-3 Heap (\code{Heap23}).
#' @param quiet If \code{FALSE}, display progress messages on screen.
#' @return Modified version of graph with additonal \code{flow} column added.
#'
#' @export
#' @examples
#' graph <- weight_streetnet (hampi)
#' from <- sample (graph$from_id, size = 10)
#' to <- sample (graph$to_id, size = 5)
#' to <- to [!to %in% from]
#' flows <- matrix (10 * runif (length (from) * length (to)),
#'                  nrow = length (from))
#' graph <- dodgr_flows (graph, from = from, to = to, flows = flows)
#' # graph then has an additonal 'flows` column of aggregate flows along all
#' # edges
dodgr_flows <- function (graph, from, to, flows, directed = TRUE,
                         wt_profile = "bicycle", heap = 'BHeap', quiet = TRUE)
{
    if (missing (graph) & (!missing (from) | !missing (to)))
        graph <- graph_from_pts (from, to, expand = 0.1,
                                 wt_profile = wt_profile, quiet = quiet)

    if ("flow" %in% names (graph))
        warning ("graph already has a 'flow' column; ",
                  "this will be overwritten")

    hps <- get_heap (heap, graph)
    heap <- hps$heap
    graph <- hps$graph

    if (!quiet)
        message ("Converting network to dodgr graph ... ",
                 appendLF = FALSE)
    graph <- dodgr_convert_graph (graph, components = FALSE)
    xy <- graph$xy
    graph <- graph$graph
    vert_map <- make_vert_map (graph)
    # vert_map$vert is char vertex ID; vert_map$id is 0-indexed integer

    if (!is.matrix (flows))
        flow <- as.matrix (flows)
    if (!(length (from) == 1 | nrow (flows) == length (from)))
        stop ("flows must have number of rows equal to length of from")
    if (!(length (to) == 1 | ncol (flows) == length (to)))
        stop ("flows must have number of columns equal to length of to")

    from_index <- get_tofrom_index (vert_map, xy, from) # 0-indexed
    to_index <- get_tofrom_index (vert_map, xy, to)

    if (!quiet)
        message ("done\nAggregating flows ... ", appendLF = FALSE)

    graph$flow <- rcpp_aggregate_flows (graph, vert_map, from_index, to_index,
                                        flows, heap)
    return (graph)
}

#' merge_directed_flows
#'
#' The \code{dodgr_flows} function returns a column of aggregated flows directed
#' along each edge of a graph, so the aggregated flow from vertex A to vertex B
#' will not necessarily equal that from B to A, and the total flow in both
#' directions will be the sum of flow from A to B plus that from B to A. This
#' function converts a directed graph to undirected form through reducing all
#' pairs of directed edges to a single edge, and aggregating flows from both
#' directions.
#'
#' @param graph A graph containing a \code{flow} column as returned from
#' \code{dodgr_flows}
#' @return An equivalent graph in which all directed edges have been reduced to
#' single, undirected edges, and all directed flows aggregated to undirected
#' flows.
#' @export
merge_directed_flows <- function (graph)
{
    flows <- rcpp_merge_flows (graph)
    indx <- which (flows > 0)
    graph <- graph [indx, , drop = FALSE]
    graph$flow <- flows [indx]
    return (graph)
}
