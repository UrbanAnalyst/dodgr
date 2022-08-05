#' merge_directed_graph
#'
#' Merge directed edges into equivalent undirected values by aggregating across
#' directions. This function is primarily intended to aid visualisation of
#' directed graphs, particularly visualising the results of the
#' \link{dodgr_flows_aggregate} and \link{dodgr_flows_disperse} functions, which
#' return columns of aggregated flows directed along each edge of a graph.
#'
#' @param graph A undirected graph in which directed edges of the input graph
#' have been merged through aggregation to yield a single, undirected edge
#' between each pair of vertices.
#' @param col_names Names of columns to be merged through aggregation. Values
#' for these columns in resultant undirected graph will be aggregated from
#' directed values.
#' @return An equivalent graph in which all directed edges have been reduced to
#' single, undirected edges, and all values of the specified column(s) have been
#' aggregated across directions to undirected values.
#' @export
#' @family misc
#' @examples
#' graph <- weight_streetnet (hampi)
#' from <- sample (graph$from_id, size = 10)
#' to <- sample (graph$to_id, size = 5)
#' to <- to [!to %in% from]
#' flows <- matrix (10 * runif (length (from) * length (to)),
#'     nrow = length (from)
#' )
#' graph <- dodgr_flows_aggregate (graph, from = from, to = to, flows = flows)
#' # graph then has an additonal 'flows` column of aggregate flows along all
#' # edges. These flows are directed, and can be aggregated to equivalent
#' # undirected flows on an equivalent undirected graph with:
#' graph_undir <- merge_directed_graph (graph)
#' # This graph will only include those edges having non-zero flows, and so:
#' nrow (graph)
#' nrow (graph_undir) # the latter is much smaller
merge_directed_graph <- function (graph, col_names = c ("flow")) {

    # auto-detect either flow or centrality as col_names:
    if (length (col_names) == 1) {
        if (col_names == "flow" && !"flow" %in% names (graph) &&
            "centrality" %in% names (graph)) {
            col_names <- "centrality"
        }
    }
    if (!all (col_names %in% names (graph))) {
        stop (paste0 (
            "col_names [",
            paste (col_names, collapse = ", "),
            "] do not match columns in graph"
        ))
    }

    gr_cols <- dodgr_graph_cols (graph)
    graph2 <- convert_graph (graph, gr_cols) # nolint
    res <- lapply (col_names, function (i) {
        graph2$merge <- graph [[i]]
        rcpp_merge_cols (graph2)
    })
    res <- do.call (cbind, res)
    index <- which (rowSums (res) > 0)
    graph <- graph [index, , drop = FALSE] # nolint
    for (i in seq (col_names)) {
        graph [[col_names [i]]] <- res [index, i]
    } # nolint
    class (graph) <- c (class (graph), "dodgr_merged")

    attr (graph, "hash") <- digest::digest (graph [[gr_cols$edge_id]])

    return (graph)
}
