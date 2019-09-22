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
    dodgr_contract_graph (graph, unique (pts))
}

#' dodgr_flows_aggregate
#'
#' Aggregate flows throughout a network based on an input matrix of flows
#' between all pairs of `from` and `to` points.
#'
#' @param graph `data.frame` or equivalent object representing the network
#' graph (see Details)
#' @param from Vector or matrix of points **from** which aggregate flows are to
#' be calculated (see Details)
#' @param to Vector or matrix of points **to** which aggregate flows are to be
#' calculated (see Details)
#' @param flows Matrix of flows with `nrow(flows)==length(from)` and
#' `ncol(flows)==length(to)`.
#' @param contract If `TRUE`, calculate flows on contracted graph before
#' mapping them back on to the original full graph (recommended as this will
#' generally be much faster).
#' @param heap Type of heap to use in priority queue. Options include
#' Fibonacci Heap (default; `FHeap`), Binary Heap (`BHeap`),
#' `Radix`, Trinomial Heap (`TriHeap`), Extended Trinomial Heap
#' (`TriHeapExt`, and 2-3 Heap (`Heap23`).
#' @param tol Relative tolerance below which flows towards `to` vertices are not
#' considered. This will generally have no effect, but can provide speed gains
#' when flow matrices represent spatial interaction models, in which case this
#' parameter effectively reduces the radius from each `from` point over which
#' flows are aggregated. To remove any such effect, set `tol = 0`.
#' @param quiet If `FALSE`, display progress messages on screen.
#' @return Modified version of graph with additonal `flow` column added.
#'
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
#'
#' # The following code can be used to convert the resultant graph to an `sf`
#' # object suitable for plotting
#' \dontrun{
#' geoms <- dodgr_to_sfc (graph_undir)
#' gc <- dodgr_contract_graph (graph_undir)
#' gsf <- sf::st_sf (geoms)
#' gsf$flow <- gc$flow
#'
#' # example of plotting with the 'mapview' package
#' library (mapview)
#' flow <- gsf$flow / max (gsf$flow)
#' ncols <- 30
#' cols <- colorRampPalette (c ("lawngreen", "red")) (ncols) [ceiling (ncols * flow)]
#' mapview (gsf, color = cols, lwd = 10 * flow)
#' }
#'
#' # An example of flow aggregation across a generic (non-OSM) highway,
#' # represented as the `routes_fast` object of the \pkg{stplanr} package,
#' # which is a SpatialLinesDataFrame containing commuter densities along
#' # components of a street network.
#' \dontrun{
#' library (stplanr)
#' # merge all of the 'routes_fast' lines into a single network
#' r <- overline (routes_fast, attrib = "length", buff_dist = 1)
#' r <- sf::st_as_sf (r)
#' # then extract the start and end points of each of the original 'routes_fast'
#' # lines and use these for routing with `dodgr`
#' l <- lapply (routes_fast@lines, function (i)
#'              c (sp::coordinates (i) [[1]] [1, ],
#'                 tail (sp::coordinates (i) [[1]], 1)))
#' l <- do.call (rbind, l)
#' xy_start <- l [, 1:2]
#' xy_end <- l [, 3:4]
#' # Then just specify a generic OD matrix with uniform values of 1:
#' flows <- matrix (1, nrow = nrow (l), ncol = nrow (l))
#' # We need to specify both a `type` and `id` column for the
#' # \link{weight_streetnet} function.
#' r$type <- 1
#' r$id <- seq (nrow (r))
#' graph <- weight_streetnet (r, type_col = "type", id_col = "id",
#'                            wt_profile = 1)
#' f <- dodgr_flows_aggregate (graph, from = xy_start, to = xy_end, flows = flows)
#' # Then merge directed flows and convert to \pkg{sf} for plotting as before:
#' f <- merge_directed_flows (f)
#' geoms <- dodgr_to_sfc (f)
#' gc <- dodgr_contract_graph (f)
#' gsf <- sf::st_sf (geoms)
#' gsf$flow <- gc$flow
#' # sf plot:
#' plot (gsf ["flow"])
#' }
#' @export
dodgr_flows_aggregate <- function (graph, from, to, flows, contract = FALSE,
                                   heap = "BHeap", tol = 1e-12, quiet = TRUE)
{
    if ("flow" %in% names (graph))
        warning ("graph already has a 'flow' column; ",
                  "this will be overwritten")

    if (any (is.na (flows))) {
        flows [is.na (flows)] <- 0
    }
    hps <- get_heap (heap, graph)
    heap <- hps$heap
    graph <- hps$graph

    gr_cols <- dodgr_graph_cols (graph)

    # change from and to just to check conformity
    tp <- attr (graph, "turn_penalty")
    tp <- ifelse (is.null (tp), 0, tp)
    if (!missing (from))
    {
        # remove any routing points not in edge start nodes:
        from <- nodes_arg_to_pts (from, graph)
        if (is (graph, "dodgr_streetnet_sc") & tp > 0)
            from <- remap_verts_with_turn_penalty (graph, from, from = TRUE)
        from <- from [from %in% graph [[gr_cols$from]] ]
    }
    if (!missing (to))
    {
        # remove any routing points not in edge end nodes:
        to <- nodes_arg_to_pts (to, graph)
        if (is (graph, "dodgr_streetnet_sc") & tp > 0)
            to <- remap_verts_with_turn_penalty (graph, to, from = FALSE)
        to <- to [to %in% graph [[gr_cols$to]] ]
    }

    if (contract)
    {
        graph_full <- graph
        graph <- contract_graph_with_pts (graph, from, to)
        hashc <- get_hash (graph, hash = FALSE)
        fname_c <- file.path (tempdir (),
                              paste0 ("dodgr_edge_map_", hashc, ".Rds"))
        if (!file.exists (fname_c))
            stop ("something went wrong extracting the edge_map ... ") # nocov
        edge_map <- readRDS (fname_c)
    }

    vert_map <- make_vert_map (graph, gr_cols)

    index_id <- get_index_id_cols (graph, gr_cols, vert_map, from)
    from_index <- index_id$index - 1 # 0-based
    index_id <- get_index_id_cols (graph, gr_cols, vert_map, to)
    to_index <- index_id$index - 1 # 0-based

    if (!is.matrix (flows))
        flows <- matrix (flows, nrow = length (from_index))

    graph2 <- convert_graph (graph, gr_cols)

    if (!quiet)
        message ("\nAggregating flows ... ", appendLF = FALSE)

    dirtxt <- get_random_prefix ()
    rcpp_flows_aggregate_par (graph2, vert_map, from_index, to_index,
                              flows, tol, dirtxt, heap)
    f <- list.files (tempdir (), full.names = TRUE)
    files <- f [grep (dirtxt, f)]
    graph$flow <- rcpp_aggregate_files (files, nrow (graph))
    junk <- file.remove (files) # nolint

    if (contract) # map contracted flows back onto full graph
        graph <- uncontract_graph (graph, edge_map, graph_full)

    return (graph)
}

get_random_prefix <- function (n = 5)
{
    charvec <- c (letters, LETTERS, 0:9)
    prefix <- paste0 (sample (charvec, n, replace = TRUE), collapse = "")
    file.path (tempdir (), paste0 ("flow_", prefix))
}

#' dodgr_flows_disperse
#'
#' Disperse flows throughout a network based on a input vectors of origin points
#' and associated densities
#'
#' @param graph `data.frame` or equivalent object representing the network
#' graph (see Details)
#' @param from Vector or matrix of points **from** which aggregate dispersed
#' flows are to be calculated (see Details)
#' @param dens Vectors of densities correponsing to the `from` points
#' @param k Width coefficient of exponential diffusion function defined as
#' `exp(-d/k)`, in units of distance column of `graph` (metres by default). Can
#' also be a vector with same length as `from`, giving dispersal coefficients
#' from each point. If value of `k<0` is given, a standard logistic polynomial
#' will be used.
#' @param contract If `TRUE`, calculate flows on contracted graph before
#' mapping them back on to the original full graph (recommended as this will
#' generally be much faster).
#' @param heap Type of heap to use in priority queue. Options include
#' Fibonacci Heap (default; `FHeap`), Binary Heap (`BHeap`),
#' `Radix`, Trinomial Heap (`TriHeap`), Extended Trinomial Heap
#' (`TriHeapExt`, and 2-3 Heap (`Heap23`).
#' @param tol Relative tolerance below which dispersal is considered to have
#' finished. This parameter can generally be ignored; if in doubt, its effect
#' can be removed by setting `tol = 0`.
#' @param quiet If `FALSE`, display progress messages on screen.
#' @return Modified version of graph with additonal `flow` column added.
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
dodgr_flows_disperse <- function (graph, from, dens, k = 500, contract = FALSE, 
                                  heap = 'BHeap', tol = 1e-12, quiet = TRUE)
{
    if ("flow" %in% names (graph))
        warning ("graph already has a 'flow' column; ",
                  "this will be overwritten")

    if (!(length (k) == 1 | length (k) == length (from)))
        stop ("'k' must be either single value or vector ",
              "of same length as 'from'")
    if (length (k) == 1)
        k <- rep (k, length (from))

    if (any (is.na (dens))) {
        dens [is.na (dens)] <- 0
    }
    hps <- get_heap (heap, graph)
    heap <- hps$heap
    graph <- hps$graph

    gr_cols <- dodgr_graph_cols (graph)

    tp <- attr (graph, "turn_penalty")
    tp <- ifelse (is.null (tp), 0, tp)
    if (!missing (from))
    {
        # remove any routing points not in edge start nodes:
        from <- nodes_arg_to_pts (from, graph)
        if (is (graph, "dodgr_streetnet_sc") & tp > 0)
            from <- remap_verts_with_turn_penalty (graph, from, from = TRUE)
        from <- from [from %in% graph [[gr_cols$from]] ]
    }

    if (contract)
    {
        graph_full <- graph
        graph <- contract_graph_with_pts (graph, from)
        hashc <- get_hash (graph, hash = FALSE)
        fname_c <- file.path (tempdir (),
                              paste0 ("dodgr_edge_map_", hashc, ".Rds"))
        if (!file.exists (fname_c))
            stop ("something went wrong extracting the edge_map ... ") # nocov
        edge_map <- readRDS (fname_c)
    }

    vert_map <- make_vert_map (graph, gr_cols)

    index_id <- get_index_id_cols (graph, gr_cols, vert_map, from)
    from_index <- index_id$index - 1 # 0-based

    if (!is.matrix (dens))
        dens <- as.matrix (dens)

    graph2 <- convert_graph (graph, gr_cols)

    if (!quiet)
        message ("\nAggregating flows ... ", appendLF = FALSE)

    # parallel results are dumped in tempdir
    dirtxt <- get_random_prefix ()
    rcpp_flows_disperse_par (graph2, vert_map, from_index,
                             k, dens, tol, dirtxt, heap)
    f <- list.files (tempdir (), full.names = TRUE)
    files <- f [grep (dirtxt, f)]
    graph$flow <- rcpp_aggregate_files (files, nrow (graph))
    junk <- file.remove (files) # nolint

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
#' @param graph A graph containing a `flow` column as returned from
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

    attr (graph, "hash") <- digest::digest (graph [[gr_cols$edge_id]])

    return (graph)
}
