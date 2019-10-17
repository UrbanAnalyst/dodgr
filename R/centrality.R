#' dodgr_centrality
#'
#' Calculate betweenness centrality for a `dodgr` network, in either vertex- or
#' edge-based form.
#'
#' @param graph `data.frame` or equivalent object representing the network
#' graph (see Details)
#' @param contract If `TRUE`, centrality is calculated on contracted graph
#' before mapping back on to the original full graph. Note that for street
#' networks, in particular those obtained from the \pkg{osmdata} package, vertex
#' placement is effectively arbitrary except at junctions; centrality for such
#' graphs should only be calculated between the latter points, and thus
#' `contract` should always be `TRUE`.
#' @param edges If `TRUE`, centrality is calculated for graph edges, returning
#' the input `graph` with an additional `centrality` column; otherwise
#' centrality is calculated for vertices, returning the equivalent of
#' `dodgr_vertices(graph)`, with an additional vertex-based `centrality` column.
#' @param heap Type of heap to use in priority queue. Options include
#' Fibonacci Heap (default; `FHeap`), Binary Heap (`BHeap`),
#' `Radix`, Trinomial Heap (`TriHeap`), Extended Trinomial Heap
#' (`TriHeapExt`, and 2-3 Heap (`Heap23`).
#' @return Modified version of graph with additonal `centrality` column added.
#'
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
#' cols <- colorRampPalette (c ("lawngreen", "red")) (ncols) [ceiling (ncols * centrality)]
#' mapview (graph_sf, color = cols, lwd = 10 * centrality)
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
#' # Convert to a 'dodgr' network, for which we need to specify both a `type` and
#' # `id` column.
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
#' @export
dodgr_centrality <- function (graph, contract = TRUE, edges = TRUE,
                              heap = "BHeap")
{
    if ("centrality" %in% names (graph))
        warning ("graph already has a 'centrality' column; ",
                  "this will be overwritten")

    hps <- get_heap (heap, graph)
    heap <- hps$heap
    graph <- hps$graph

    gr_cols <- dodgr_graph_cols (graph)

    if (contract & methods::is (graph, "dodgr_contracted"))
        contract <- FALSE
    if (contract & !methods::is (graph, "dodgr_contracted"))
    {
        graph_full <- graph
        graph <- dodgr_contract_graph (graph)
        hashc <- get_hash (graph, hash = FALSE)
        fname_c <- file.path (tempdir (),
                              paste0 ("dodgr_edge_map_", hashc, ".Rds"))
        if (!file.exists (fname_c))
            stop ("something went wrong extracting the edge_map ... ") # nocov
        edge_map <- readRDS (fname_c)
    }

    vert_map <- make_vert_map (graph, gr_cols)

    graph2 <- convert_graph (graph, gr_cols)

    # centrality calculation, done in parallel with each thread dumping results
    # to files in tempdir()
    if (edges)
    {
        dirtxt <- get_random_prefix ("centrality_edge")
        rcpp_centrality_edge (graph2, vert_map, heap, dirtxt)
    } else
    {
        dirtxt <- get_random_prefix ("centrality_vert")
        rcpp_centrality_vertex (graph2, vert_map, heap, dirtxt)
    }

    # aggregate results from the threads:
    f <- list.files (tempdir (), full.names = TRUE)
    files <- f [grep (dirtxt, f)]
    if (edges)
        centrality <- rcpp_aggregate_files (files, nrow (graph))
    else
    {
        v <- dodgr_vertices (graph)
        centrality <- rcpp_aggregate_files (files, nrow (v))
    }
    junk <- file.remove (files) # nolint

    # attach result to edge or vertex objects:
    if (edges)
    {
        graph$centrality <- centrality
        if (contract)
            graph <- uncontract_graph (graph, edge_map, graph_full)
        res <- graph
    } else
    {
        res <- dodgr_vertices (graph)
        res$centrality <- centrality
    }

    return (res)
}
