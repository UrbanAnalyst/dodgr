#' Create a map of `dodgr` flows.
#'
#' Create a map of the output of \link{dodgr_flows_aggregate} or
#' \link{dodgr_flows_disperse}
#'
#' @param net A street network with a `flow` column obtained from
#' \link{dodgr_flows_aggregate} or \link{dodgr_flows_disperse}
#' @param bbox If given, scale the map to this bbox, otherwise use entire extend
#' of `net`
#' @param linescale Maximal thickness of plotted lines
#'
#' @note `net` should be first passed through `merge_directed_graph`
#' prior to plotting, otherwise lines for different directions will be overlaid.
#' @family misc
#' @export
#' @examples
#' graph <- weight_streetnet (hampi)
#' from <- sample (graph$from_id, size = 10)
#' to <- sample (graph$to_id, size = 5)
#' to <- to [!to %in% from]
#' flows <- matrix (
#'     10 * runif (length (from) * length (to)),
#'     nrow = length (from)
#' )
#' graph <- dodgr_flows_aggregate (graph, from = from, to = to, flows = flows)
#' # graph then has an additonal 'flows` column of aggregate flows along all
#' # edges. These flows are directed, and can be aggregated to equivalent
#' # undirected flows on an equivalent undirected graph with:
#' graph_undir <- merge_directed_graph (graph)
#' \dontrun{
#' dodgr_flowmap (graph_undir)
#' }
dodgr_flowmap <- function (net, bbox = NULL, linescale = 1) {

    if (!"flow" %in% names (net)) {
        net$flow <- 1
    }
    gr_cols <- dodgr_graph_cols (net)
    names (net) [gr_cols$xfr] <- "from_lon"
    names (net) [gr_cols$yfr] <- "from_lat"
    names (net) [gr_cols$xto] <- "to_lon"
    names (net) [gr_cols$yto] <- "to_lat"

    if (is.null (bbox)) {
        bbox <- c (
            min (net$from_lon), min (net$from_lat),
            max (net$from_lon), max (net$from_lat)
        )
    }

    xlims <- c (bbox [1], bbox [3])
    ylims <- c (bbox [2], bbox [4])
    cols <- colorRampPalette (c ("lawngreen", "red")) (30)
    plot (NULL, xlim = xlims, ylim = ylims, xlab = "lon", ylab = "lat")
    net <- net [which (net$flow > 0), ]
    net$flow <- net$flow / max (net$flow)
    ncols <- 30
    cols <- colorRampPalette (c ("lawngreen", "red")) (ncols)
    cols <- cols [ceiling (net$flow * ncols)]

    # Suppress 'no visible binding' lint messages.
    from_lon <- from_lat <- to_lon <- to_lat <- NULL
    with (net, segments (from_lon, from_lat, to_lon, to_lat,
        col = cols, lwd = linescale * net$flow
    ))
}
