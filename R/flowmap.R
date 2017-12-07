#' dodgr_flowmap
#'
#' Map the output of \link{dodgr_flows}
#'
#' @param net A street network with a \code{flow} column obtained from
#' \code{dodgr_flows}
#' @param filename Name of \code{.png} file to write to
#' @param linescale Maximal thickness of plotted lines
#'
#' @note \code{net} should be first passed through \code{merge_directed_flows}
#' prior to plotting, otherwise lines for different directions will be overlaid.
#' @export
dodgr_flowmap <- function (net, filename, linescale = 5)
{
    fmax <- max (net$flow)
    bb <- c (min (net$from_lon), min (net$from_lat),
             max (net$from_lon), max (net$from_lat)) %>%
            osmplotr::get_bbox ()

    cols <- colorRampPalette (c ("lawngreen", "red")) (30)
    map <- osmplotr::osm_basemap (bb, bg = "gray10") %>%
        osmplotr::add_colourbar (colour = "gray95", zlims = c (0, fmax),
                                 col = cols) %>%
        osmplotr::add_axes (colour = "gray95")

    net$flow <- linescale * net$flow / fmax
    from_lon <- from_lat <- to_lon <- to_lat <- flow <- NULL # suppress warn
    map <- map +
        ggplot2::geom_segment (ggplot2::aes (x = from_lon, y = from_lat,
                           xend = to_lon, yend = to_lat,
                           colour = flow, size = flow),
                      size = net$flow, data = net) +
        ggplot2::scale_colour_gradient (low = "lawngreen", high = "red",
                                        guide = "none",
                                        limits = c (0, max (net$flow)))
    osmplotr::print_osm_map (map, file = paste0 (filename, ".png"), dpi = 72)
}
