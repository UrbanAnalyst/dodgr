#' dodgr_flowmap
#'
#' Map the output of \link{dodgr_flows_aggregate} or \link{dodgr_flows_disperse}
#'
#' @param net A street network with a \code{flow} column obtained from
#' \link{dodgr_flows_aggregate} or \link{dodgr_flows_disperse}
#' @param bbox If given, scale the map to this bbox, otherwise use entire extend
#' of \code{net}
#' @param linescale Maximal thickness of plotted lines
#'
#' @note \code{net} should be first passed through \code{merge_directed_flows}
#' prior to plotting, otherwise lines for different directions will be overlaid.
#' @export
dodgr_flowmap <- function (net, bbox = NULL, linescale = 1)
{
    if (!"flow" %in% names (net))
        net$flow <- 1
    fmax <- max (net$flow)
    if (is.null (bbox))
        bbox <- c (min (net$from_lon), min (net$from_lat),
                 max (net$from_lon), max (net$from_lat))

    xlims <- c (bbox [1], bbox [3])
    ylims <- c (bbox [2], bbox [4])
    cols <- colorRampPalette (c ("lawngreen", "red")) (30)
    plot.new ()
    plot (NULL, xlim = xlims, ylim = ylims, xlab = "lon", ylab = "lat")
    net <- net [which (net$flow > 0), ]
    net$flow <- net$flow / max (net$flow)
    ncols <- 30
    cols <- colorRampPalette (c ("lawngreen", "red")) (ncols)
    cols <- cols [ceiling (net$flow * ncols)]

    with (net, segments (from_lon, from_lat, to_lon, to_lat,
                            col = cols, lwd = linescale * net$flow))
}
