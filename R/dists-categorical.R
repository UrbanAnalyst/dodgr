#' Cumulative distances along different edge categories
#'
#' @inheritParams dodgr_dists
#' @param graph `data.frame` or equivalent object representing the network
#' graph which must have a column named "edge_type" which labels categories of
#' edge types along which categorical distances are to be aggregated (see
#' Note).
#' @param proportions_only If `FALSE`, return distance matrices for full
#' distances and for each edge category; if `TRUE`, return single vector of
#' proportional distances, like the `summary` function applied to full
#' results. See Note.
#' @param dlimit If no value to `to` is given, distances are aggregated from
#' each `from` point out to the specified distance limit (in the same units as
#' the edge distances of the input graph). `dlimit` only has any effect if `to`
#' is not specified, in which case the `proportions_only` argument has no
#' effect.
#' @return If `to` is specified, a list of distance matrices of equal dimensions
#' (length(from), length(to)), the first of which ("distance") holds the final
#' distances, while the rest are one matrix for each unique value of
#' "edge_type", holding the distances traversed along those types of edges only.
#' Otherwise, a single matrix of total distances along all ways from each point
#' out to the specified value of `dlimit`, along with distances along each of
#' the different kinds of ways specified in the "edge_type" column of the input
#' graph.
#'
#' @note The "edge_type" column in the graph can contain any kind of discrete or
#' categorical values, although integer values of 0 are not permissible. `NA`
#' values are ignored. The function requires one full distance
#' matrix to be stored for each category of "edge_type" (unless
#' `proportions_only = TRUE`). It is wise to keep numbers of discrete types as
#' low as possible, especially for large distance matrices.
#'
#' @note Setting the `proportions_only` flag to `TRUE` may be advantageous for
#' large jobs, because this avoids construction of the full matrices. This may
#' speed up calculations, but perhaps more importantly it may make possible
#' calculations which would otherwise require distance matrices too large to be
#' directly stored.
#'
#' @note Calculations are not able to be interrupted (for example, by `Ctrl-C`),
#' and can only be stopped by killing the R process.
#' @examples
#' # Prepare a graph for categorical routing by including an "edge_type" column
#' graph <- weight_streetnet (hampi, wt_profile = "foot")
#' graph <- graph [graph$component == 1, ]
#' graph$edge_type <- graph$highway
#' # Define start and end points for categorical distances; using all vertices
#' # here.
#' length (unique (graph$edge_type)) # Number of categories
#' v <- dodgr_vertices (graph)
#' from <- to <- v$id [1:100]
#' d <- dodgr_dists_categorical (graph, from, to)
#' class (d)
#' length (d)
#' sapply (d, dim)
#' # 9 distance matrices, all of same dimensions, first of which is standard
#' # distance matrix
#' # s <- summary (d) # print summary as proportions along each "edge_type"
#' # or directly calculate proportions only
#' dodgr_dists_categorical (graph, from, to,
#'     proportions_only = TRUE
#' )
#'
#' # The 'dlimit' parameter can be used to calculate total distances along each
#' # category of edges from a set of points:
#' dlimit <- 2000 # in metres
#' d <- dodgr_dists_categorical (graph, from, dlimit = dlimit)
#' dim (d) # length(from), length(unique(edge_type)) + 1
#' @family distances
#' @export
dodgr_dists_categorical <- function (graph,
                                     from = NULL,
                                     to = NULL,
                                     proportions_only = FALSE,
                                     dlimit = NULL,
                                     heap = "BHeap",
                                     quiet = TRUE) {

    if (!"edge_type" %in% names (graph)) {
        stop ("graph must have a column named 'edge_type'")
    }
    if (is.integer (graph$edge_type) && any (graph$edge_type == 0L)) {
        stop ("graphs with integer edge_type columns may not contain 0s")
    }
    if (is.null (to)) {
        if (is.null (dlimit)) {
            stop ("'dlimit' must be specified if no 'to' points are given.")
        }
        if (!(is.numeric (dlimit) && length (dlimit) == 1L)) {
            stop ("'dlimit' must be a single number.")
        }
    }

    graph <- tbl_to_df (graph)

    edge_type <- graph$edge_type

    hps <- get_heap (heap, graph)
    heap <- hps$heap
    graph <- hps$graph

    gr_cols <- dodgr_graph_cols (graph)
    if (is.na (gr_cols$from) || is.na (gr_cols$to)) {
        scols <- find_spatial_cols (graph)
        graph$from_id <- scols$xy_id$xy_fr_id
        graph$to_id <- scols$xy_id$xy_to_id
        gr_cols <- dodgr_graph_cols (graph)
    }
    is_spatial <- is_graph_spatial (graph)
    if (!is_spatial) {
        stop ("Categorical distances only implemented for spatial graphs")
    }

    vert_map <- make_vert_map (graph, gr_cols, is_spatial)

    from_index <- fill_to_from_index (graph, vert_map, gr_cols, from, from = TRUE)
    to_index <- fill_to_from_index (graph, vert_map, gr_cols, to, from = FALSE)

    if (get_turn_penalty (graph) > 0.0) {
        if (methods::is (graph, "dodgr_contracted")) {
            warning (
                "graphs with turn penalties should be submitted in full, ",
                "not contracted form;\nsubmitting contracted graphs may ",
                "produce unexpected behaviour."
            )
        }
        graph <- create_compound_junctions (graph)$graph
        edge_type <- graph$edge_type

        # remap any 'from' and 'to' vertices to compound junction versions:
        vert_map <- make_vert_map (graph, gr_cols, is_spatial)

        from_index <- remap_tf_index_for_tp (from_index, vert_map, from = TRUE)
        to_index <- remap_tf_index_for_tp (to_index, vert_map, from = FALSE)
    }

    graph <- convert_graph (graph, gr_cols)
    edge_type_table <- table (edge_type)
    graph$edge_type <- match (edge_type, names (edge_type_table))
    graph$edge_type [is.na (graph$edge_type)] <- 0L

    if (!quiet) {
        message ("Calculating shortest paths ... ", appendLF = FALSE)
    }

    if (!is.null (to)) {

        d <- rcpp_get_sp_dists_categorical (
            graph,
            vert_map,
            from_index$index,
            to_index$index,
            heap,
            proportions_only
        )

        n <- length (to)

        if (!proportions_only) {

            res <- process_categorical_dmat (
                d,
                from_index,
                to_index,
                vert_map,
                edge_type_table
            )

        } else {

            res <- apply (d, 2, sum)
            res [2:length (res)] <- res [2:length (res)] / res [1]
            res <- res [-1]
            names (res) <- names (edge_type_table)
        }
    } else {

        d <- rcpp_get_sp_dists_cat_threshold (
            graph,
            vert_map,
            from_index$index,
            dlimit,
            heap
        )

        res <- process_threshold_dmat (d, from_index, vert_map, edge_type_table)
    }


    return (res)
}

process_categorical_dmat <- function (d, from_index, to_index, vert_map,
                                      edge_type_table) {

    n <- length (to_index$index)

    if (is.null (from_index$id)) {
        rnames <- vert_map$vert
    } else {
        rnames <- from_index$id
    }
    if (is.null (to_index$id)) {
        cnames <- vert_map$vert
    } else {
        cnames <- to_index$id
    }

    # compound turn-penalty junctions:
    rnames <- gsub ("\\_start$", "", rnames)
    cnames <- gsub ("\\_end$", "", cnames)

    d0 <- d [, seq (n)]
    rownames (d0) <- rnames
    colnames (d0) <- cnames

    d0 <- list ("distances" = d0)
    d <- lapply (seq_along (edge_type_table), function (i) {
        index <- i * n + seq (n) - 1
        res <- d [, index]
        rownames (res) <- rnames
        colnames (res) <- cnames
        return (res)
    })
    names (d) <- names (edge_type_table)

    res <- c (d0, d)
    class (res) <- append (class (res), "dodgr_dists_categorical")

    return (res)
}

process_threshold_dmat <- function (d, from_index, vert_map, edge_type_table) {

    if (is.null (from_index$id)) {
        rownames (d) <- vert_map$vert
    } else {
        rownames (d) <- from_index$id
    }
    rownames (d) <- gsub ("\\_start$", "", rownames (d))

    res <- data.frame (d)
    names (res) <- c ("distance", names (edge_type_table))

    return (res)
}

#' Transform a result from 'dodgr_dists_categorical' to summary statistics
#'
#' @param object A 'dodgr_dists_categorical' object
#' @param ... Extra parameters currently not used
#' @return The summary statistics (invisibly)
#' @family misc
#' @export
summary.dodgr_dists_categorical <- function (object, ...) {

    d0 <- object$distances # first list item
    sum_d0 <- sum (d0, na.rm = TRUE)
    object <- object [-1]

    dprop <- vapply (
        object, function (i) {
            sum (i, na.rm = TRUE) / sum_d0
        },
        numeric (1)
    )

    message ("Proportional distances along each kind of edge:")
    for (i in seq_along (dprop)) {
        message (
            "  ", names (dprop) [i],
            ": ", round (dprop [i], digits = 4)
        )
    }

    invisible (dprop)
}
