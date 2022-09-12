#' Deduplicate edges in a graph
#'
#' Graph may have duplicated edges, particularly when extracted as
#' \link{dodgr_streetnet} objects. This function de-duplicates any repeated
#' edges, reducing weighted distances and times to the minimal values from all
#' duplicates.
#' @param graph Any 'dodgr' graph or network.
#' @return A potentially modified version of graph, with any formerly duplicated
#' edges reduces to single rows containing minimal weighted distances and times.
#' @family conversion
#' @export
dodgr_deduplicate_graph <- function (graph) {

    gr_cols <- dodgr_graph_cols (graph)
    fr_col <- names (graph) [gr_cols$from]
    to_col <- names (graph) [gr_cols$to]
    d_col <- names (graph) [gr_cols$d_weighted]
    t_col <- names (graph) [gr_cols$time_weighted]
    t_col <- ifelse (is.na (t_col), "", t_col)

    fr_to <- paste0 (graph [[fr_col]], "-", graph [[to_col]])
    index <- which (duplicated (fr_to))

    if (length (index) == 0L) {
        return (graph)
    }

    has_times <- nzchar (t_col)

    res <- rcpp_deduplicate (graph, fr_col, to_col, d_col, t_col)

    if (has_times) {

        n <- as.integer (nrow (res) / 2)
        res_t <- res [seq (n) + n, ]
        res <- res [seq (n), ]
        ft <- paste0 (res$from, "-", res$to)
        ft_t <- paste0 (res_t$from, "-", res_t$to)
        index_t <- match (ft, ft_t)
        res$t <- NA # necessary to define as numeric
        res$t [index_t] <- res_t$d
    }

    graph <- graph [-index, ]
    fr_to <- fr_to [-index]

    fr_to_res <- paste0 (res$from, "-", res$to)

    index_to_gr <- match (fr_to_res, fr_to)
    graph [[d_col]] [index_to_gr] <- res$d
    if (has_times) {
        graph [[t_col]] [index_to_gr] <- res$t
    }

    return (graph)
}

#' Issue a warning if graph has duplicated edges
#'
#' Currently only used in centrality; see #186.
#'
#' @return Logical flag indicating whether or not graph has duplicated edges;
#' `TRUE` denotes a graph with duplicated edges; `FALSE` denotes a graph which
#' passes this check.
#' @noRd
duplicated_edge_check <- function (graph, proceed = FALSE) {

    gr_cols <- dodgr_graph_cols (graph)
    fr_name <- names (graph) [gr_cols$from]
    to_name <- names (graph) [gr_cols$to]
    fr <- graph [[fr_name]]
    to <- graph [[to_name]]

    fr_to <- paste0 (fr, "-", to)

    has_duplicates <- any (duplicated (fr_to))

    if (has_duplicates) {

        msg <- paste0 (
            "Graph has duplicated edges. Only the first will be used here,\n",
            "but it is better to remove them first with the ",
            "'dodgr_deduplicate_graph() function."
        )
        if (interactive ()) {
            message (msg)
            if (!proceed) {
                x <- readline ("Do you want to proceed (y/n)? ")
                if (tolower (substring (x, 1, 1) != "y")) { # nocov
                    stop ("Okay, we'll stop there", call. = FALSE)
                }
            }
        } else {
            warning (msg)
        }
    }

    invisible (has_duplicates)
}
