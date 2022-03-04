#' dodgr_centrality
#'
#' Calculate betweenness centrality for a 'dodgr' network, in either vertex- or
#' edge-based form.
#'
#' @param graph 'data.frame' or equivalent object representing the network
#' graph (see Details)
#' @param contract If 'TRUE', centrality is calculated on contracted graph
#' before mapping back on to the original full graph. Note that for street
#' networks, in particular those obtained from the \pkg{osmdata} package, vertex
#' placement is effectively arbitrary except at junctions; centrality for such
#' graphs should only be calculated between the latter points, and thus
#' 'contract' should always be 'TRUE'.
#' @param edges If 'TRUE', centrality is calculated for graph edges, returning
#' the input 'graph' with an additional 'centrality' column; otherwise
#' centrality is calculated for vertices, returning the equivalent of
#' 'dodgr_vertices(graph)', with an additional vertex-based 'centrality' column.
#' @param column Column of graph defining the edge properties used to calculate
#' centrality (see Note).
#' @param vert_wts Optional vector of length equal to number of vertices
#' (`nrow(dodgr_vertices(graph))`), to enable centrality to be calculated in
#' weighted form, such that centrality measured from each vertex will be
#' weighted by the specified amount.
#' @param dist_threshold If not 'NULL', only calculate centrality for each point
#' out to specified threshold. Setting values for this will result in
#' approximate estimates for centrality, yet with considerable gains in
#' computational efficiency. For sufficiently large values, approximations will
#' be accurate to within some constant multiplier. Appropriate values can be
#' established via the \link{estimate_centrality_threshold} function.
#' @param heap Type of heap to use in priority queue. Options include
#' Fibonacci Heap (default; 'FHeap'), Binary Heap ('BHeap'),
#' Trinomial Heap ('TriHeap'), Extended Trinomial Heap
#' ('TriHeapExt', and 2-3 Heap ('Heap23').
#' @return Modified version of graph with additional 'centrality' column added.
#'
#' @note The `column` parameter is by default `d_weighted`, meaning centrality
#' is calculated by routing according to weighted distances. Other possible
#' values for this parameter are
#' \itemize{
#' \item `d` for unweighted distances
#' \item `time` for unweighted time-based routing
#' \item `time_weighted` for weighted time-based routing
#' }
#'
#' @note Centrality is calculated by default using parallel computation with the
#' maximal number of available cores or threads. This number can be reduced by
#' specifying a value via
#' `RcppParallel::setThreadOptions (numThreads = <desired_number>)`.
#'
#' @family centrality
#' @export
#' @examples
#' graph_full <- weight_streetnet (hampi)
#' graph <- dodgr_contract_graph (graph_full)
#' graph <- dodgr_centrality (graph)
#' # 'graph' is then the contracted graph with an additional 'centrality' column
#' # Same calculation via 'igraph':
#' igr <- dodgr_to_igraph (graph)
#' library (igraph)
#' cent <- edge_betweenness (igr)
#' identical (cent, graph$centrality) # TRUE
#' # Values of centrality between all junctions in the contracted graph can then
#' # be mapped back onto the original full network by "uncontracting":
#' graph_full <- dodgr_uncontract_graph (graph)
#' # For visualisation, it is generally necessary to merge the directed edges to
#' # form an equivalent undirected graph. Conversion to 'sf' format via
#' # 'dodgr_to_sf()' is also useful for many visualisation routines.
#' graph_sf <- merge_directed_graph (graph_full) %>%
#'     dodgr_to_sf ()
#'
#' \dontrun{
#' library (mapview)
#' centrality <- graph_sf$centrality / max (graph_sf$centrality)
#' ncols <- 30
#' cols <- c ("lawngreen", "red")
#' cols <- colorRampPalette (cols) (ncols) [ceiling (ncols * centrality)]
#' mapview (graph_sf, color = cols, lwd = 10 * centrality)
#' }
#'
#' # An example of flow aggregation across a generic (non-OSM) highway,
#' # represented as the 'routes_fast' object of the \pkg{stplanr} package,
#' # which is a SpatialLinesDataFrame containing commuter densities along
#' # components of a street network.
#' \dontrun{
#' library (stplanr)
#' # merge all of the 'routes_fast' lines into a single network
#' r <- overline (routes_fast, attrib = "length", buff_dist = 1)
#' r <- sf::st_as_sf (r)
#' # Convert to a 'dodgr' network, for which we need to specify both a 'type'
#' # and 'id' column.
#' r$type <- 1
#' r$id <- seq (nrow (r))
#' graph_full <- weight_streetnet (r, type_col = "type", id_col = "id",
#'                                 wt_profile = 1)
#' # convert to contracted form, retaining junction vertices only, and append
#' # 'centrality' column
#' graph <- dodgr_contract_graph (graph_full) %>%
#'     dodgr_centrality ()
#' #' expand back to full graph; merge directed flows; and convert result to
#' # 'sf'-format for plotting
#' graph_sf <- dodgr_uncontract_graph (graph) %>%
#'     merge_directed_graph () %>%
#'     dodgr_to_sf ()
#' plot (graph_sf ["centrality"])
#' }
#'
dodgr_centrality <- function (graph,
                              contract = TRUE,
                              edges = TRUE,
                              column = "d_weighted",
                              vert_wts = NULL,
                              dist_threshold = NULL,
                              heap = "BHeap") {

    column <- match.arg (column, c ("d_weighted",
                                    "d",
                                    "time",
                                    "time_weighted"))
    if ("centrality" %in% names (graph))
        warning ("graph already has a 'centrality' column; ",
                  "this will be overwritten")

    if (!is.null (vert_wts) &
        (!is.vector (vert_wts) |
        !is.numeric (vert_wts) |
        length (vert_wts) != nrow (dodgr_vertices (graph)))) {

        stop ("vert_wts must be a vector of same length as ",
              "nrow (dodgr_vertices (graph))")
    }

    if (is.null (dist_threshold))
        dist_threshold <- .Machine$double.xmax

    hps <- get_heap (heap, graph)
    heap <- hps$heap
    graph <- hps$graph

    gr_cols <- dodgr_graph_cols (graph)

    if (contract & methods::is (graph, "dodgr_contracted"))
        contract <- FALSE
    graph_full <- edge_map <- NULL
    if (contract & !methods::is (graph, "dodgr_contracted")) {
        graph_full <- graph
        graph <- dodgr_contract_graph (graph)
        hashc <- get_hash (graph, hash = FALSE)
        fname_c <- fs::path (fs::path_temp (),
                             paste0 ("dodgr_edge_map_", hashc, ".Rds"))
        if (!fs::file_exists (fname_c))
            stop ("something went wrong extracting the edge_map ... ") # nocov
        edge_map <- readRDS (fname_c)
    }

    vert_map <- make_vert_map (graph, gr_cols)
    if (!is.null (vert_wts))
        vert_map$vert_wts <- vert_wts

    graph2 <- convert_graph (graph, gr_cols)
    if (column != "d_weighted") {
        gr_cols2 <- dodgr_graph_cols (graph2)
        column <- gr_cols2 [[column]]
        graph2$d_weighted <- graph2 [[column]]
    }

    # final '0' is for sampling calculation to estimate speed - non-zero values
    # used only in 'estimate_centrality_time'
    centrality <- rcpp_centrality (graph2,
                                   vert_map,
                                   heap,
                                   dist_threshold,
                                   edges,
                                   0)

    # attach result to edge or vertex objects:
    if (edges) {
        graph$centrality <- centrality
        if (contract)
            graph <- uncontract_graph (graph, edge_map, graph_full)
        res <- graph
    } else {
        res <- dodgr_vertices (graph)
        res$centrality <- centrality
    }

    # re-cache graph:
    if (is_dodgr_cache_on ())
        attr (graph, "px") <- cache_graph (graph, gr_cols$edge_id)

    return (res)
}

#' estimate_centrality_threshold
#'
#' Estimate a value for the 'dist_threshold' parameter of the
#' \link{dodgr_centrality} function. Providing distance thresholds to this
#' function generally provides considerably speed gains, and results in
#' approximations of centrality. This function enables the determination of
#' values of 'dist_threshold' corresponding to specific degrees of accuracy.
#'
#' @inheritParams dodgr_centrality
#' @param tolerance Desired maximal degree of inaccuracy in centrality estimates
#' - values will be accurate to within this amount, subject to a constant
#' scaling factor. Note that threshold values increase non-linearly with
#' decreasing values of 'tolerance'
#' @return A single value for 'dist_threshold' giving the required tolerance.
#'
#' @note This function may take some time to execute. While running, it displays
#' ongoing information on screen of estimated values of 'dist_threshold' and
#' associated errors. Thresholds are progressively increased until the error is
#' reduced below the specified tolerance.
#'
#' @family centrality
#' @export
estimate_centrality_threshold <- function (graph, tolerance = 0.001) {
    # nocov start - can't be tested on any sample data

    # estimate absolute maximum distance in graph
    v <- dodgr_vertices (graph)
    vs <- sample (v$id, 1000)
    dabsmax <- max (dodgr_dists (graph, from = vs, to = vs), na.rm = TRUE)

    nverts <- min (10000, nrow (v))
    if (nverts == 10000)
        graphs <- dodgr_sample (graph, nverts = nverts)
    else
        graphs <- graph

    # estimate max dist
    vs <- dodgr_vertices (graphs)
    vsample <- sample (vs$id, 1000)
    dmax <- max (dodgr_dists (graphs, from = vsample, to = vsample),
                 na.rm = TRUE)
    if (dmax > (0.75 * dabsmax)) {
        message ("dist_threshold approaches size of graph;\n",
                 "Recommended value of 'dist_threshold' remains 'NULL',\n",
                 "with centrality calculated across entire graph")
        return (NULL)
    }

    d <- 100
    mult <- 1.1
    ss <- 9999
    g_old <- dodgr_centrality (graphs, dist_threshold = d)
    while (ss > tolerance) {
        d <- signif (d * mult, 2)
        g <- dodgr_centrality (graphs, dist_threshold = d)
        mod <- stats::lm (g$centrality ~ g_old$centrality)
        ss <- sum (mod$residuals ^ 2) / sum (g_old$centrality ^ 2)
        message ("d = ", round (d), "; error = ",
                 formatC (ss, format = "f", digits = 4))
        g_old <- g

        if (d > dmax) {
            message ("re-sampling graph to ", nverts, " vertices")
            nverts <- signif (nverts * mult, 2)
            graphs <- dodgr_sample (graph, nverts)
            if (nverts > nrow (v))
                break
        }
    }
    if (nverts > nrow (v)) {
        message ("Failed to converge within tolerance;\n",
                 "Recommended value of 'dist_threshold' remains 'NULL',\n",
                 "with centrality calculated across entire graph")
        d <- NULL
    } else
        message ("converged on distance threshold of ", d)

    return (d)
    # nocov end
}

#' estimate_centrality_time
#'
#' The 'dodgr' centrality functions are designed to be applied to potentially
#' very large graphs, and may take considerable time to execute. This helper
#' function estimates how long a centrality function may take for a given graph
#' and given value of 'dist_threshold' estimated via the
#' \link{estimate_centrality_threshold} function.
#'
#' @inheritParams dodgr_centrality
#' @return An estimated calculation time for calculating centrality for the
#' given value of 'dist_threshold'
#'
#' @note This function may take some time to execute. While running, it displays
#' ongoing information on screen of estimated values of 'dist_threshold' and
#' associated errors. Thresholds are progressively increased until the error is
#' reduced below the specified tolerance.
#'
#' @family centrality
#' @export
estimate_centrality_time <- function (graph,
                                      contract = TRUE,
                                      edges = TRUE,
                                      dist_threshold = NULL,
                                      heap = "BHeap") {

    # copies all code from dodgr_centrality, but uses the otherwise non-exposed
    # 'sample' parameter passed through to C++ routines
    if (is.null (dist_threshold))
        dist_threshold <- .Machine$double.xmax

    hps <- get_heap (heap, graph)
    heap <- hps$heap
    graph <- hps$graph

    gr_cols <- dodgr_graph_cols (graph)

    if (contract & methods::is (graph, "dodgr_contracted"))
        contract <- FALSE
    if (contract & !methods::is (graph, "dodgr_contracted")) {
        graph <- dodgr_contract_graph (graph)
        hashc <- get_hash (graph, hash = FALSE)
        fname_c <- fs::path (fs::path_temp (),
                             paste0 ("dodgr_edge_map_", hashc, ".Rds"))
        if (!fs::file_exists (fname_c))
            stop ("something went wrong extracting the edge_map ... ") # nocov
    }

    vert_map <- make_vert_map (graph, gr_cols)

    graph2 <- convert_graph (graph, gr_cols)

    # final '0' is for sampling calculation to estimate speed - non-zero values
    # used only in 'estimate_centrality_time'
    st <- system.time (
                       x <- rcpp_centrality (graph2, vert_map, heap,
                                             dist_threshold, edges, 100)
                       ) [3]

    # convert to estimated time for full graph
    st <- st * nrow (graph) / 100

    hh <- floor (st / 3600)
    st <- st - 3600 * hh
    mm <- floor (st / 60)
    ss <- round (st - 60 * mm)

    hh <- sprintf ("%02d", hh)
    mm <- sprintf ("%02d", mm)
    ss <- sprintf ("%02d", ss)
    res <- paste0 (hh, ":", mm, ":", ss)
    message ("Estimated time to calculate centrality for full graph is ",
             res)
    invisible (res)
}
