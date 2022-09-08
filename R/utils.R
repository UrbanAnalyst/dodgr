null_to_na <- function (x) {

    if (length (x) == 0) {
        x <- NA
    }
    return (x)
}

#' Issue a warning if graph has duplicated edges
#'
#' Currently only used in centrality; see #186.
#'
#' @return Logical flag indicating whether or not graph has duplicated edges;
#' `TRUE` denotes a graph with duplicated edges; `FALSE` denotes a graph which
#' passes this check.
#' @noRd
duplicated_edge_check <- function (graph) {

    gr_cols <- dodgr_graph_cols (graph)
    fr_name <- names (graph) [gr_cols$from]
    to_name <- names (graph) [gr_cols$to]
    fr <- graph [[fr_name]]
    to <- graph [[to_name]]

    fr_to <- paste0 (fr, "-", to)

    has_duplicates <- any (duplicated (fr_to))

    if (has_duplicates) {

        if (interactive ()) {
            message (paste0 (
                "Graph has duplicated edges. Only the first will be used here, ",
                "but it is better to manually remove them first."
            ))
            x <- readline ("Do you want to proceed (y/n)? ")
            if (tolower (substring (x, 1, 1) != "y")) { # nocov
                stop ("Okay, we'll stop there", call. = FALSE)
            }
        } else {
            warning ("graph has duplicated edges; only the first will be used here.")
        }
    }

    invisible (has_duplicates)
}
