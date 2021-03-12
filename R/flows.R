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
#' Trinomial Heap (`TriHeap`), Extended Trinomial Heap
#' (`TriHeapExt`, and 2-3 Heap (`Heap23`).
#' @param tol Relative tolerance below which flows towards `to` vertices are not
#' considered. This will generally have no effect, but can provide speed gains
#' when flow matrices represent spatial interaction models, in which case this
#' parameter effectively reduces the radius from each `from` point over which
#' flows are aggregated. To remove any such effect, set `tol = 0`.
#' @param quiet If `FALSE`, display progress messages on screen.
#' @inheritParams dodgr_flows_si
#' @return Modified version of graph with additional `flow` column added.
#'
#' @note Spatial Interaction models are often fitted through trialling a range
#' of values of 'k'. The specification above allows fitting multiple values of
#' 'k' to be done with a single call, in a way that is far more efficient than
#' making multiple calls. A matrix of 'k' values may be entered, with each
#' column holding a different vector of values, one for each 'from' point. For a
#' matrix of 'k' values having 'n' columns, the return object will be a modified
#' version in the input 'graph', with an additional 'n' columns, named 'flow1',
#' 'flow2', ... up to 'n'. These columns must be subsequently matched by the
#' user back on to the corresponding columns of the matrix of 'k' values.
#'
#' @note The `norm_sums` parameter should be used whenever densities at origins
#' and destinations are absolute values, and ensures that the sum of resultant
#' flow values throughout the entire network equals the sum of densities at all
#' origins. For example, with `norm_sums = TRUE` (the default), a flow from a
#' single origin with density one to a single destination along two edges will
#' allocate flows of one half to each of those edges, such that the sum of flows
#' across the network will equal one, or the sum of densities from all origins.
#' The `norm_sums = TRUE` option is appropriate where densities are relative
#' values, and ensures that each edge maintains relative proportions. In the
#' above example, flows along each of two edges would equal one, for a network
#' sum of two, or greater than the sum of densities.
#'
#' Flows are calculated by default using parallel computation with the maximal
#' number of available cores or threads. This number can be reduced by
#' specifying a value via
#' `RcppParallel::setThreadOptions (numThreads = <desired_number>)`.
#'
#' @examples
#' graph <- weight_streetnet (hampi)
#' from <- sample (graph$from_id, size = 10)
#' to <- sample (graph$to_id, size = 5)
#' to <- to [!to %in% from]
#' flows <- matrix (10 * runif (length (from) * length (to)),
#'                  nrow = length (from))
#' graph <- dodgr_flows_aggregate (graph, from = from, to = to, flows = flows)
#' # graph then has an additonal 'flows' column of aggregate flows along all
#' # edges. These flows are directed, and can be aggregated to equivalent
#' # undirected flows on an equivalent undirected graph with:
#' graph_undir <- merge_directed_graph (graph)
#' # This graph will only include those edges having non-zero flows, and so:
#' nrow (graph); nrow (graph_undir) # the latter is much smaller
#'
#' # The following code can be used to convert the resultant graph to an `sf`
#' # object suitable for plotting
#' \dontrun{
#' gsf <- dodgr_to_sf (graph_undir)
#'
#' # example of plotting with the 'mapview' package
#' library (mapview)
#' flow <- gsf$flow / max (gsf$flow)
#' ncols <- 30
#' cols <- c ("lawngreen", "red")
#' colranmp <- colorRampPalette (cols) (ncols) [ceiling (ncols * flow)]
#' mapview (gsf, color = colranmp, lwd = 10 * flow)
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
#' f <- dodgr_flows_aggregate (graph, from = xy_start, to = xy_end,
#'                             flows = flows)
#' # Then merge directed flows and convert to \pkg{sf} for plotting as before:
#' f <- merge_directed_graph (f)
#' geoms <- dodgr_to_sfc (f)
#' gc <- dodgr_contract_graph (f)
#' gsf <- sf::st_sf (geoms)
#' gsf$flow <- gc$flow
#' # sf plot:
#' plot (gsf ["flow"])
#' }
#' @export
dodgr_flows_aggregate <- function (graph,
                                   from,
                                   to,
                                   flows,
                                   contract = FALSE,
                                   heap = "BHeap",
                                   tol = 1e-12,
                                   norm_sums = TRUE,
                                   quiet = TRUE) {

    if (any (is.na (flows))) {
        flows [is.na (flows)] <- 0
    }
    hps <- get_heap (heap, graph)
    heap <- hps$heap
    graph <- hps$graph

    if (contract) {
        graph_full <- graph
        graph <- contract_graph_with_pts (graph, from, to)
        hashc <- get_hash (graph, hash = FALSE)
        fname_c <- file.path (tempdir (),
                              paste0 ("dodgr_edge_map_", hashc, ".Rds"))
        if (!file.exists (fname_c))
            stop ("something went wrong extracting the edge_map ... ") # nocov
        edge_map <- readRDS (fname_c)
    }

    g <- prepare_graph (graph, from, to)
    if (!is.matrix (flows))
        flows <- matrix (flows, nrow = length (g$from_index))

    if (!quiet)
        message ("\nAggregating flows ... ", appendLF = FALSE)

    graph$flow <- rcpp_flows_aggregate_par (g$graph, g$vert_map,
                                            g$from_index, g$to_index,
                                            flows, norm_sums, tol, heap)

    if (contract) # map contracted flows back onto full graph
        graph <- uncontract_graph (graph, edge_map, graph_full)

    return (graph)
}

#' dodgr_flows_disperse
#'
#' Disperse flows throughout a network based on a input vectors of origin points
#' and associated densities
#'
#' @inheritParams dodgr_flows_aggregate
#' @param graph `data.frame` or equivalent object representing the network
#' graph (see Details)
#' @param from Vector or matrix of points **from** which aggregate dispersed
#' flows are to be calculated (see Details)
#' @param dens Vectors of densities corresponding to the `from` points
#' @param k Width coefficient of exponential diffusion function defined as
#' `exp(-d/k)`, in units of distance column of `graph` (metres by default). Can
#' also be a vector with same length as `from`, giving dispersal coefficients
#' from each point. If value of `k<0` is given, a standard logistic polynomial
#' will be used.
#' @param tol Relative tolerance below which dispersal is considered to have
#' finished. This parameter can generally be ignored; if in doubt, its effect
#' can be removed by setting `tol = 0`.
#' @return Modified version of graph with additional `flow` column added.
#'
#' @note Spatial Interaction models are often fitted through trialling a range
#' of values of 'k'. The specification above allows fitting multiple values of
#' 'k' to be done with a single call, in a way that is far more efficient than
#' making multiple calls. A matrix of 'k' values may be entered, with each
#' column holding a different vector of values, one for each 'from' point. For a
#' matrix of 'k' values having 'n' columns, the return object will be a modified
#' version in the input 'graph', with an additional 'n' columns, named 'flow1',
#' 'flow2', ... up to 'n'. These columns must be subsequently matched by the
#' user back on to the corresponding columns of the matrix of 'k' values.
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
#' graph_undir <- merge_directed_graph (graph)
dodgr_flows_disperse <- function (graph,
                                  from,
                                  dens,
                                  k = 500,
                                  contract = FALSE,
                                  heap = "BHeap",
                                  tol = 1e-12,
                                  quiet = TRUE) {

    res <- check_k (k, from)
    k <- res$k
    nk <- res$nk

    if (any (is.na (dens))) {
        dens [is.na (dens)] <- 0
    }

    hps <- get_heap (heap, graph)
    heap <- hps$heap
    graph <- hps$graph

    if (contract) {
        graph_full <- graph
        graph <- contract_graph_with_pts (graph, from)
        hashc <- get_hash (graph, hash = FALSE)
        fname_c <- file.path (tempdir (),
                              paste0 ("dodgr_edge_map_", hashc, ".Rds"))
        if (!file.exists (fname_c))
            stop ("something went wrong extracting the edge_map ... ") # nocov
        edge_map <- readRDS (fname_c)
    }

    g <- prepare_graph (graph, from)

    if (!is.matrix (dens))
        dens <- as.matrix (dens)

    if (!quiet)
        message ("\nAggregating flows ... ", appendLF = FALSE)

    f <- rcpp_flows_disperse_par (g$graph, g$vert_map, g$from_index,
                                  k, dens, tol, heap)
    if (nk == 1)
        graph$flow <- f
    else {
        flowmat <- data.frame (matrix (f, ncol = nk))
        names (flowmat) <- paste0 ("flow", seq (nk))
        graph <- cbind (graph, flowmat)
    }

    if (contract) # map contracted flows back onto full graph
        graph <- uncontract_graph (graph, edge_map, graph_full)

    return (graph)
}

#' dodgr_flows_si
#'
#' Aggregate flows throughout a network based using an exponential Spatial
#' Interaction (SI) model between a specified set of origin and destination
#' points, and associated vectors of densities.
#'
#' @inheritParams dodgr_flows_aggregate
#' @param k Width of exponential spatial interaction function (exp (-d / k)),
#' in units of 'd', specified in one of 3 forms: (i) a single value; (ii) a
#' vector of independent values for each origin point (with same length as
#' 'from' points); or (iii) an equivalent matrix with each column holding values
#' for each 'from' point, so 'nrow(k)==length(from)'. See Note.
#' @param dens_from Vector of densities at origin ('from') points
#' @param dens_to Vector of densities at destination ('to') points
#' @param norm_sums Standardise sums from all origin points, so sum of flows
#' throughout entire network equals sum of densities from all origins (see
#' Note).
#' @return Modified version of graph with additional `flow` column added.
#'
#' @note Spatial Interaction models are often fitted through trialling a range
#' of values of 'k'. The specification above allows fitting multiple values of
#' 'k' to be done with a single call, in a way that is far more efficient than
#' making multiple calls. A matrix of 'k' values may be entered, with each
#' column holding a different vector of values, one for each 'from' point. For a
#' matrix of 'k' values having 'n' columns, the return object will be a modified
#' version in the input 'graph', with an additional 'n' columns, named 'flow1',
#' 'flow2', ... up to 'n'. These columns must be subsequently matched by the
#' user back on to the corresponding columns of the matrix of 'k' values.
#'
#' @note The `norm_sums` parameter should be used whenever densities at origins
#' and destinations are absolute values, and ensures that the sum of resultant
#' flow values throughout the entire network equals the sum of densities at all
#' origins. For example, with `norm_sums = TRUE` (the default), a flow from a
#' single origin with density one to a single destination along two edges will
#' allocate flows of one half to each of those edges, such that the sum of flows
#' across the network will equal one, or the sum of densities from all origins.
#' The `norm_sums = TRUE` option is appropriate where densities are relative
#' values, and ensures that each edge maintains relative proportions. In the
#' above example, flows along each of two edges would equal one, for a network
#' sum of two, or greater than the sum of densities.
#'
#' With `norm_sums = TRUE`, the sum of network flows (`sum(output$flow)`) should
#' equal the sum of origin densities (`sum(dens_from)`). This may nevertheless
#' not always be the case, because origin points may simply be too far from any
#' destination (`to`) points for an exponential model to yield non-zero values
#' anywhere in a network within machine tolerance. Such cases may result in sums
#' of output flows being less than sums of input densities.
#'
#' @examples
#' graph <- weight_streetnet (hampi)
#' from <- sample (graph$from_id, size = 10)
#' to <- sample (graph$to_id, size = 5)
#' to <- to [!to %in% from]
#' flows <- matrix (10 * runif (length (from) * length (to)),
#'                  nrow = length (from))
#' graph <- dodgr_flows_aggregate (graph, from = from, to = to, flows = flows)
#' # graph then has an additonal 'flows' column of aggregate flows along all
#' # edges. These flows are directed, and can be aggregated to equivalent
#' # undirected flows on an equivalent undirected graph with:
#' graph_undir <- merge_directed_graph (graph)
#' # This graph will only include those edges having non-zero flows, and so:
#' nrow (graph); nrow (graph_undir) # the latter is much smaller
#' @export
dodgr_flows_si <- function (graph,
                            from,
                            to,
                            k = 500,
                            dens_from = NULL,
                            dens_to = NULL,
                            contract = FALSE,
                            norm_sums = TRUE,
                            heap = "BHeap",
                            tol = 1e-12,
                            quiet = TRUE) {

    if (missing (from))
        stop ("'from' must be provided for spatial interaction models.")

    hps <- get_heap (heap, graph)
    heap <- hps$heap
    graph <- hps$graph

    res <- check_k (k, from)
    k <- res$k
    nk <- res$nk

    if (contract) {
        graph_full <- graph
        graph <- contract_graph_with_pts (graph, from, to)
        hashc <- get_hash (graph, hash = FALSE)
        fname_c <- file.path (tempdir (),
                              paste0 ("dodgr_edge_map_", hashc, ".Rds"))
        if (!file.exists (fname_c))
            stop ("something went wrong extracting the edge_map ... ") # nocov
        edge_map <- readRDS (fname_c)
    }

    g <- prepare_graph (graph, from, to)

    if (!quiet)
        message ("\nAggregating flows ... ", appendLF = FALSE)

    f <- rcpp_flows_si (g$graph, g$vert_map, g$from_index, g$to_index,
                        k, dens_from, dens_to, norm_sums, tol, heap)

    if (nk == 1)
        graph$flow <- f
    else {
        flowmat <- data.frame (matrix (f, ncol = nk))
        names (flowmat) <- paste0 ("flow", seq (nk))
        graph <- cbind (graph, flowmat)
    }


    if (contract) # map contracted flows back onto full graph
        graph <- uncontract_graph (graph, edge_map, graph_full)

    return (graph)
}

check_k <- function (k,
                     from) {

    nk <- 1

    if (is.data.frame (k))
        k <- as.matrix (k)

    if (is.matrix (k)) {
        if (nrow (k) != length (from))
            stop ("nrow(k) must equal length of 'from' points")
        nk <- ncol (k)
    } else if (is.numeric (k)) {
        if (length (k) == 1)
            k <- rep (k, length (from))
        else if (length (k) != length (from)) {
            # convert to matrix
            nk <- length (k)
            k <- array (rep (k, each = length (from)),
                        dim = c (length (from), nk))
        }
    } else
        stop ("'k' must be either a single value, a vector, or a matrix")

    list (k = k, nk = nk)
}

# transform input graph and (from, to) arguments to standard forms for passing
# to C++ routines
prepare_graph <- function (graph,
                           from,
                           to) {

    if ("flow" %in% names (graph))
        warning ("graph already has a 'flow' column; ",
                  "this will be overwritten")

    gr_cols <- dodgr_graph_cols (graph)
    vert_map <- make_vert_map (graph, gr_cols)

    # change from and to just to check conformity
    tp <- attr (graph, "turn_penalty")
    tp <- ifelse (is.null (tp), 0, tp)

    # remove any routing points not in edge start nodes:
    from <- nodes_arg_to_pts (from, graph)
    if (methods::is (graph, "dodgr_streetnet_sc") & tp > 0)
        from <- remap_verts_with_turn_penalty (graph, from, from = TRUE)
    from <- from [from %in% graph [[gr_cols$from]] ] # nolint
    index_id <- get_index_id_cols (graph, gr_cols, vert_map, from)
    from_index <- index_id$index - 1 # 0-based

    to_index <- NULL
    if (!missing (to)) {
        # remove any routing points not in edge end nodes:
        to <- nodes_arg_to_pts (to, graph)
        if (methods::is (graph, "dodgr_streetnet_sc") & tp > 0)
            to <- remap_verts_with_turn_penalty (graph, to, from = FALSE)
        to <- to [to %in% graph [[gr_cols$to]] ] # nolint
        index_id <- get_index_id_cols (graph, gr_cols, vert_map, to)
        to_index <- index_id$index - 1 # 0-based
    }

    graph2 <- convert_graph (graph, gr_cols)

    list (graph = graph2, vert_map = vert_map,
          from_index = from_index, to_index = to_index)
}

get_random_prefix <- function (prefix = "flow",
                               n = 5) {

    charvec <- c (letters, LETTERS, 0:9)
    rand_prefix <- paste0 (sample (charvec, n, replace = TRUE), collapse = "")
    file.path (tempdir (), paste0 (prefix, "_", rand_prefix))
}

nodes_arg_to_pts <- function (nodes,
                              graph) {

    if (!is.matrix (nodes))
        nodes <- as.matrix (nodes)
    if (ncol (nodes) == 2) {
        verts <- dodgr_vertices (graph)
        nodes <- verts$id [match_pts_to_graph (verts, nodes)]
    }
    return (nodes)
}


# keep from and to routing points in contracted graph
contract_graph_with_pts <- function (graph,
                                     from,
                                     to) {

    pts <- NULL
    if (!missing (from))
        pts <- c (pts, from)
    if (!missing (to))
        pts <- c (pts, to)
    dodgr_contract_graph (graph, unique (pts))
}
