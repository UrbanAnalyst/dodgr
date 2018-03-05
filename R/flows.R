nodes_arg_to_pts <- function (nodes, graph)
{
    if (!is.matrix (nodes))
        nodes <- as.matrix (nodes)
    if (ncol (nodes) == 2)
    {
        verts <- dodgr_vertices (graph)
        nodes <- verts$id [match_pts_to_graph (verts, nodes)]
    }
    return (nodes)
}


# keep from and to routing points in contracted graph
contract_graph_with_pts <- function (graph, from, to)
{
    pts <- NULL
    if (!missing (from))
        pts <- c (pts, from)
    if (!missing (to))
        pts <- c (pts, to)
    graph_full <- graph
    graph <- dodgr_contract_graph (graph, unique (pts))
    graph$graph_full <- graph_full
    return (graph)
}

# map contracted flows back onto full graph
uncontract_graph <- function (graph, edge_map, graph_full)
{
    indx_to_full <- match (edge_map$edge_old, graph_full$edge_id)
    indx_to_contr <- match (edge_map$edge_new, graph$edge_id)
    # edge_map only has the contracted edges; flows from the original
    # non-contracted edges also need to be inserted
    edges <- graph$edge_id [which (!graph$edge_id %in% edge_map$edge_new)]
    indx_to_full <- c (indx_to_full, match (edges, graph_full$edge_id))
    indx_to_contr <- c (indx_to_contr, match (edges, graph$edge_id))
    graph_full$flow <- 0
    graph_full$flow [indx_to_full] <- graph$flow [indx_to_contr]

    return (graph_full)
}

#' dodgr_flows_aggregate
#'
#' Aggregate flows throughout a network based on an input matrix of flows
#' between all pairs of \code{from} and \code{to} points.
#'
#' @param graph \code{data.frame} or equivalent object representing the network
#' graph (see Details)
#' @param from Vector or matrix of points **from** which aggregate flows are to
#' be calculated (see Details)
#' @param to Vector or matrix of points **to** which aggregate flows are to be
#' calculated (see Details)
#' @param flows Matrix of flows with \code{nrow(flows)==length(from)} and
#' \code{ncol(flows)==length(to)}.
#' @param wt_profile Name of weighting profile for street networks (one of foot,
#' horse, wheelchair, bicycle, moped, motorcycle, motorcar, goods, hgv, psv;
#' only used if \code{graph} is not provided, in which case a street network is
#' downloaded and correspondingly weighted).
#' @param contract If \code{TRUE}, calculate flows on contracted graph before
#' mapping them back on to the original full graph (recommended as this will
#' generally be much faster).
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
#' graph <- dodgr_flows_aggregate (graph, from = from, to = to, flows = flows)
#' # graph then has an additonal 'flows` column of aggregate flows along all
#' # edges. These flows are directed, and can be aggregated to equivalent
#' # undirected flows on an equivalent undirected graph with:
#' graph_undir <- merge_directed_flows (graph)
#' # This graph will only include those edges having non-zero flows, and so:
#' nrow (graph); nrow (graph_undir) # the latter is much smaller
dodgr_flows_aggregate <- function (graph, from, to, flows, wt_profile =
                                   "bicycle", contract = FALSE, heap = 'BHeap',
                                   quiet = TRUE)
{
    if (missing (graph) & (!missing (from) | !missing (to)))
        graph <- graph_from_pts (from, to, expand = 0.1,
                                 wt_profile = wt_profile, quiet = quiet)

    if ("flow" %in% names (graph))
        warning ("graph already has a 'flow' column; ",
                  "this will be overwritten")

    if (any (is.na (flows))) {
        flows [is.na (flows)] <- 0
    }
    hps <- get_heap (heap, graph)
    heap <- hps$heap
    graph <- hps$graph

    # change from and to just to check conformity
    if (!missing (from))
        from <- nodes_arg_to_pts (from, graph)
    if (!missing (to))
        to <- nodes_arg_to_pts (to, graph)

    if (contract)
    {
        graph <- contract_graph_with_pts (graph, from, to)
        graph_full <- graph$graph_full
        edge_map <- graph$edge_map
        graph <- graph$graph
    }

    gr_cols <- dodgr_graph_cols (graph)
    vert_map <- make_vert_map (graph, gr_cols)

    index_id <- get_index_id_cols (graph, gr_cols, vert_map, from)
    from_index <- index_id$index - 1 # 0-based
    index_id <- get_index_id_cols (graph, gr_cols, vert_map, to)
    to_index <- index_id$index - 1 # 0-based

    if (!is.matrix (flows))
        flows <- t (as.matrix (flows))

    graph2 <- convert_graph (graph, gr_cols)

    if (!quiet)
        message ("\nAggregating flows ... ", appendLF = FALSE)

    # parallel results are dumped in tempdir, which is read here with an extra
    # char to get the terminal dir separator character:
    dirtxt <- file.path (tempdir (), "a")
    dirtxt <- substr (dirtxt, 1, nchar (dirtxt) - 1)
    rcpp_flows_aggregate_par (graph2, vert_map, from_index, to_index,
                              flows, dirtxt, heap)
    files <- list.files (tempdir (), pattern = "flow_", full.names = TRUE)
    graph$flow <- rcpp_aggregate_files (files, nrow (graph))
    junk <- file.remove (files) # nolint

    if (contract) # map contracted flows back onto full graph
        graph <- uncontract_graph (graph, edge_map, graph_full)

    return (graph)
}

#' dodgr_flows_disperse
#'
#' Disperse flows throughout a network based on a input vectors of origin points
#' and associated densities
#'
#' @param graph \code{data.frame} or equivalent object representing the network
#' graph (see Details)
#' @param from Vector or matrix of points **from** which aggregate dispersed
#' flows are to be calculated (see Details)
#' @param dens Vectors of densities correponsing to the \code{from} points
#' @param wt_profile Name of weighting profile for street networks (one of foot,
#' horse, wheelchair, bicycle, moped, motorcycle, motorcar, goods, hgv, psv).
#' @param contract If \code{TRUE}, calculate flows on contracted graph before
#' mapping them back on to the original full graph (recommended as this will
#' generally be much faster).
#' @param k Width coefficient of exponential diffusion function defined as
#' \code{exp(-d/k)}.  If value of \code{k<0} is given, a standard logistic
#' polynomial will be used.
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
#' dens <- rep (1, length (from)) # Uniform densities
#' graph <- dodgr_flows_disperse (graph, from = from, dens = dens)
#' # graph then has an additonal 'flows` column of aggregate flows along all
#' # edges. These flows are directed, and can be aggregated to equivalent
#' # undirected flows on an equivalent undirected graph with:
#' graph_undir <- merge_directed_flows (graph)
dodgr_flows_disperse <- function (graph, from, dens, wt_profile = "bicycle",
                         contract = FALSE, k = 2, heap = 'BHeap', quiet = TRUE)
{
    if (missing (graph) & !missing (from))
        graph <- graph_from_pts (from, from, expand = 0.1,
                                 wt_profile = wt_profile, quiet = quiet)

    if ("flow" %in% names (graph))
        warning ("graph already has a 'flow' column; ",
                  "this will be overwritten")

    if (any (is.na (dens))) {
        dens [is.na (dens)] <- 0
    }
    hps <- get_heap (heap, graph)
    heap <- hps$heap
    graph <- hps$graph

    if (!missing (from))
        from <- nodes_arg_to_pts (from, graph)

    if (contract)
    {
        graph <- contract_graph_with_pts (graph, from, from)
        graph_full <- graph$graph_full
        edge_map <- graph$edge_map
        graph <- graph$graph
    }

    gr_cols <- dodgr_graph_cols (graph)
    vert_map <- make_vert_map (graph, gr_cols)

    index_id <- get_index_id_cols (graph, gr_cols, vert_map, from)
    from_index <- index_id$index - 1 # 0-based

    if (!is.matrix (dens))
        dens <- as.matrix (dens)

    graph2 <- convert_graph (graph, gr_cols)

    if (!quiet)
        message ("\nAggregating flows ... ", appendLF = FALSE)

    graph$flow <- rcpp_flows_disperse (graph2, vert_map, from_index,
                                       k, dens, heap)

    if (contract) # map contracted flows back onto full graph
        graph <- uncontract_graph (graph, edge_map, graph_full)

    return (graph)
}

#' merge_directed_flows
#'
#' The \link{dodgr_flows_aggregate} and \link{dodgr_flows_disperse} functions
#' return a column of aggregated flows directed along each edge of a graph, so
#' the aggregated flow from vertex A to vertex B will not necessarily equal that
#' from B to A, and the total flow in both directions will be the sum of flow
#' from A to B plus that from B to A. This function converts a directed graph to
#' undirected form through reducing all pairs of directed edges to a single
#' edge, and aggregating flows from both directions.
#'
#' @param graph A graph containing a \code{flow} column as returned from
#' \link{dodgr_flows_aggregate} or \link{dodgr_flows_disperse}
#' @return An equivalent graph in which all directed edges have been reduced to
#' single, undirected edges, and all directed flows aggregated to undirected
#' flows.
#' @export
#' @examples
#' graph <- weight_streetnet (hampi)
#' from <- sample (graph$from_id, size = 10)
#' to <- sample (graph$to_id, size = 5)
#' to <- to [!to %in% from]
#' flows <- matrix (10 * runif (length (from) * length (to)),
#'                  nrow = length (from))
#' graph <- dodgr_flows_aggregate (graph, from = from, to = to, flows = flows)
#' # graph then has an additonal 'flows` column of aggregate flows along all
#' # edges. These flows are directed, and can be aggregated to equivalent
#' # undirected flows on an equivalent undirected graph with:
#' graph_undir <- merge_directed_flows (graph)
#' # This graph will only include those edges having non-zero flows, and so:
#' nrow (graph); nrow (graph_undir) # the latter is much smaller
merge_directed_flows <- function (graph)
{
    if (!"flow" %in% names (graph))
        stop ("graph does not have any flows to merge")

    gr_cols <- dodgr_graph_cols (graph)
    graph2 <- convert_graph (graph, gr_cols)
    graph2$flow <- graph$flow

    flows <- rcpp_merge_flows (graph2)

    indx <- which (flows > 0)
    graph <- graph [indx, , drop = FALSE] #nolint
    graph$flow <- flows [indx]
    class (graph) <- c (class (graph), "dodgr_merged_flows")
    return (graph)
}
