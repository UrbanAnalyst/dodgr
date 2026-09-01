#' Calculate lists of pair-wise shortest paths between points.
#'
#' @param graph `data.frame` or equivalent object representing the network
#' graph (see Details)
#' @param from Vector or matrix of points **from** which route paths are to
#' be calculated (see Details)
#' @param to Vector or matrix of points **to** which route paths are to be
#' calculated (see Details)
#' @param vertices If `TRUE`, return lists of lists of vertices for each
#' path, otherwise return corresponding lists of edge numbers from `graph`.
#' @param pairwise If `TRUE`, calculate paths only between the ordered
#' pairs of `from` and `to`. In this case, each of these must be the
#' same length, and the output will contain paths the i-th members of each, and
#' thus also be of that length.
#' @param heap Type of heap to use in priority queue. Options include
#' Fibonacci Heap (default; `FHeap`), Binary Heap (`BHeap`),
#' `Radix`, Trinomial Heap (`TriHeap`), Extended Trinomial Heap
#' (`TriHeapExt`, and 2-3 Heap (`Heap23`).
#' @param quiet If `FALSE`, display progress messages on screen.
#' @return List of list of paths tracing all connections between nodes such that
#' if `x <- dodgr_paths (graph, from, to)`, then the path between
#' `from[i]` and `to[j]` is `x [[i]] [[j]]`. Each individual path is then a
#' vector of integers indexing into the rows of `graph` if `vertices = FALSE`,
#' or into the rows of `dodgr_vertices (graph)` if `vertices = TRUE`. The
#' returned list also has a `vertices` attribute recording the value of the
#' `vertices` argument used to generate it, so that other functions (such as
#' \link{dodgr_paths_expand}) can subsequently interpret individual paths
#' correctly without that argument having to be separately re-specified.
#'
#' @note `graph` must minimally contain four columns of `from`,
#' `to`, `dist`. If an additional column named `weight` or
#' `wt` is present, shortest paths are calculated according to values
#' specified in that column; otherwise according to `dist` values. Either
#' way, final distances between `from` and `to` points are calculated
#' according to values of `dist`. That is, paths between any pair of points
#' will be calculated according to the minimal total sum of `weight`
#' values (if present), while reported distances will be total sums of
#' `dist` values.
#'
#' The `from` and `to` columns of `graph` may be either single
#' columns of numeric or character values specifying the numbers or names of
#' graph vertices, or combinations to two columns specifying geographical
#' (longitude and latitude) coordinates. In the latter case, almost any sensible
#' combination of names will be accepted (for example, `fromx, fromy`,
#' `from_x, from_y`, or `fr_lat, fr_lon`.)
#'
#' `from` and `to` values can be either two-column matrices of
#' equivalent of longitude and latitude coordinates, or else single columns
#' precisely matching node numbers or names given in `graph$from` or
#' `graph$to`. If `to` is missing, pairwise distances are calculated
#' between all points specified in `from`. If neither `from` nor
#' `to` are specified, pairwise distances are calculated between all nodes
#' in `graph`.
#'
#' @family distances
#' @export
#' @examples
#' graph <- weight_streetnet (hampi)
#' from <- sample (graph$from_id, size = 100)
#' to <- sample (graph$to_id, size = 50)
#' dp <- dodgr_paths (graph, from = from, to = to)
#' # dp is a list with 100 items, and each of those 100 items has 30 items, each
#' # of which is a single path listing all vertiex IDs as taken from `graph`.
#'
#' # it is also possible to calculate paths between pairwise start and end
#' # points
#' from <- sample (graph$from_id, size = 5)
#' to <- sample (graph$to_id, size = 5)
#' dp <- dodgr_paths (graph, from = from, to = to, pairwise = TRUE)
#' # dp is a list of 5 items, each of which just has a single path between each
#' # pairwise from and to point.
dodgr_paths <- function (graph,
                         from,
                         to,
                         vertices = TRUE,
                         pairwise = FALSE,
                         heap = "BHeap",
                         quiet = TRUE) {

    is_contracted <- methods::is (graph, "dodgr_contracted")

    hps <- get_heap (heap, graph)
    heap <- hps$heap
    graph <- hps$graph

    gr_cols <- dodgr_graph_cols (graph)
    # cols are (edge_id, from, to, d, w, component, xfr, yfr, xto, yto)
    is_spatial <- is_graph_spatial (graph)
    vert_map <- make_vert_map (graph, gr_cols, is_spatial)

    if (missing (from)) {
        from <- vert_map$vert
    }
    from_index <- get_path_indices (graph, gr_cols, vert_map, from)

    if (missing (to)) {
        to <- vert_map$vert
    }
    to_index <- get_path_indices (graph, gr_cols, vert_map, to)

    graph <- convert_graph (graph, gr_cols)

    if (!quiet) {
        message ("Calculating shortest paths ... ", appendLF = FALSE)
    }
    if (pairwise) {
        if (length (from_index$index) != length (to_index$index)) {
            stop ("pairwise paths require from and to to have same length")
        }

        do_bidirectional <- !is_contracted &&
            nrow (vert_map) > 1000

        paths <- rcpp_get_paths_pairwise (
            graph,
            vert_map,
            from_index$index,
            to_index$index,
            heap,
            is_spatial,
            do_bidirectional = do_bidirectional
        )
    } else {
        paths <- rcpp_get_paths (
            graph,
            vert_map,
            from_index$index,
            to_index$index,
            heap
        )
    }

    # convert 1-based indices back into vertex IDs. Note both paths that can not
    # be traced and single-step paths are returned from the above as NULL. The
    # former are retained as NULL, while the following converts the latter to
    # appropriate start-end vertices.
    paths <- lapply (paths, function (i) {
        lapply (i, function (j) {
            if (is.null (j)) {
                return (j)
            } # nocov
            vert_map$vert [j]
        })
    }) # nolint


    # name path lists
    if (!is.null (from_index$id) && !is.null (to_index$id)) {
        if (!pairwise) {
            for (i in seq_along (from_index$id)) {
                names (paths [[i]]) <- paste0 (
                    from_index$id [i],
                    "-",
                    to_index$id
                )
            }
        }
        names (paths) <- from_index$id
    }

    if (!vertices) {
        graph_verts <- paste0 ("f", graph$from, "t", graph$to)

        # convert vertex IDs to corresponding sequences of edge numbers
        paths <- lapply (paths, function (i) {
            lapply (i, function (j) {
                if (length (j) > 1) {
                    indx <- 2:length (j)
                    pij <- paste0 (
                        "f", j [indx - 1],
                        "t", j [indx]
                    )
                    res <- match (pij, graph_verts)
                    res <- res [which (!is.na (res))]
                    return (if (length (res) == 0) {
                        NULL
                    } else {
                        res
                    })
                }
            })
        }) # nolint
    }

    attr (paths, "vertices") <- vertices

    return (paths)
}

get_path_indices <- function (graph, gr_cols, vert_map, to_from) {

    index_id <- get_index_id_cols (graph, gr_cols, vert_map, to_from)

    index <- index_id$index - 1 # 0-based
    if (!is.null (index_id$id)) {
        id <- index_id$id
    } else {
        id <- vert_map$vert # nocov
    }

    return (list (index = index, id = id))
}

sort_transitions <- function (graph) {

    cols <- dodgr_graph_cols (graph)

    to_neighbors <- match (graph [[cols$to]], graph [[cols$from]])
    pos <- which (!(graph [[cols$from]] %in% graph [[cols$to]]))

    if (length (pos) != 1L) {
        stop (
            "graph does not represent a single, contiguous chain of edges.",
            call. = FALSE
        )
    }

    perm <- integer (nrow (graph))
    for (i in seq_along (perm)) {
        perm [[i]] <- pos
        pos <- to_neighbors [[pos]]
    }

    perm
}

#' Expand contracted paths back on to full, uncontracted graphs.
#'
#' Paths calculated with \link{dodgr_paths} on graphs contracted with
#' \link{dodgr_contract_graph} trace and return only the intervening junction
#' vertices. This function expands contracted paths back onto full sequences
#' of edges of original, uncontracted graphs.
#'
#' @param paths A nested list of paths between pairs of (origin, destination)
#' points, generated by \link{dodgr_paths} on a graph previously contracted
#' with \link{dodgr_contract_graph}, so that `paths [[i]] [[j]]` is the path
#' between the i-th origin and j-th destination points. Each individual path
#' may be given either as a character vector of vertex IDs (as returned by
#' \link{dodgr_paths} with `vertices = TRUE`), or as an integer vector of row
#' indices into `graph_c` (as returned with `vertices = FALSE`). `paths`
#' should generally carry the `vertices` attribute attached by
#' \link{dodgr_paths}, used here to determine the type of the return value
#' (see 'Value' below); if that attribute is absent, the type of the first
#' non-`NULL` path is used instead. Individual paths which are `NULL` or
#' which trace no transitions are ignored.
#' @param graph The full, uncontracted graph from which `graph_c` was
#' generated with \link{dodgr_contract_graph}.
#' @param graph_c The contracted graph on which `paths` was calculated.
#' @param edge_map Optional edge map generated internally by
#' \link{dodgr_contract_graph}, and used to map contracted edges back on to
#' their corresponding sequences of edges in `graph`. If not given, this is
#' extracted directly from the dodgr cache.
#'
#' @return A nested list of the same structure as `paths`, in which
#' `result [[i]] [[j]]` is the expansion of `paths [[i]] [[j]]` on to the
#' full, uncontracted sequence of edges of `graph`, ordered from start to end
#' of that path. Each expanded path is returned in a form structurally
#' identical to the corresponding input: a character vector of vertex IDs
#' from `graph` if `paths [[i]] [[j]]` was itself a character vector (that
#' is, if `paths` was generated with `vertices = TRUE`), or an integer vector
#' of row indices into `graph` if `paths [[i]] [[j]]` was an integer vector
#' (`vertices = FALSE`). Entries for which the corresponding input path is
#' `NULL` or has no transitions are returned as `NULL`. The result also
#' carries the same `vertices` attribute as `paths`.
#' @family distances
#' @export
#' @examples
#' graph <- weight_streetnet (hampi)
#' graph_c <- dodgr_contract_graph (graph)
#' verts <- dodgr_vertices (graph_c)
#' from <- verts$id [1]
#' to <- verts$id [nrow (verts)]
#' paths_c <- dodgr_paths (graph_c, from = from, to = to)
#' paths <- dodgr_paths_expand (paths_c, graph, graph_c)
dodgr_paths_expand <- function (paths, graph, graph_c, edge_map = NULL) {

    vertices <- attr (paths, "vertices")

    if (is.null (edge_map)) {
        edge_map <- get_edge_map (graph)
    }

    result <- lapply (paths, function (paths_i) {
        lapply (paths_i, function (path_c) {

            vertices_i <- vertices
            if (is.null (vertices_i)) {
                vertices_i <- is.character (path_c)
            }
            min_len <- if (vertices_i) 2L else 1L

            if (is.null (path_c) || length (path_c) < min_len) {
                return (NULL)
            }
            dodgr_one_path_expand (path_c, graph, graph_c, edge_map = edge_map)
        })
    })

    attr (result, "vertices") <- vertices

    result
}

# Expand a single contracted path back on to the full, uncontracted graph.
# `path_c` is a single path along a contracted graph, as one element of the
# nested list produced by \link{dodgr_paths}, given either as a character
# vector of vertex IDs, or an integer vector of row indices into `graph_c`.
# The return value is structurally identical to `path_c`: a character vector
# of the full, uncontracted sequence of vertex IDs from `graph` if `path_c`
# was character, or an integer vector of row indices into `graph` if
# `path_c` was integer.
dodgr_one_path_expand <- function (path_c, graph, graph_c, edge_map = NULL) {

    vertices <- is.character (path_c)

    # Ensure the path contains at least one transition
    min_len <- if (vertices) 2L else 1L
    if (length (path_c) < min_len) {
        stop ("There are no transitions.")
    }

    # dodgr columns
    cols <- dodgr_graph_cols (graph)

    if (vertices) {

        # Reduce size of contracted graph
        keep <- graph_c [[cols$from]] %in% path_c |
            graph_c [[cols$to]] %in% path_c
        graph_c <- graph_c [keep, ]

        # Match path transitions to contracted graph
        path_c <- match (
            paste (path_c [-length (path_c)], path_c [-1]),
            paste (graph_c [[cols$from]], graph_c [[cols$to]])
        )

        # Verify that all edges were mapped
        if (any (is.na (path_c))) {
            stop ("Not all transitions were matched to the contracted graph.")
        }
    }

    if (!is.integer (path_c)) {
        stop ("Path must be provided as integer or character vector.")
    }

    # Path edges
    edges <- graph_c [[cols$edge_id]] [path_c]

    # Retrieve edge map if not provided
    if (is.null (edge_map)) {
        edge_map <- get_edge_map (graph)
    }

    # Convert edge map to list:
    # maps each contracted edge to a sequence of original edges
    edge_map <- edge_map [edge_map$edge_new %in% edges, ]
    edge_map <- split (edge_map$edge_old, edge_map$edge_new)

    # Add edges that map to themselves (no contraction)
    missing_edges <- setdiff (edges, names (edge_map))
    edge_map [missing_edges] <- as.list (missing_edges)

    # Match expanded edge sequence to rows in original graph, ordered from
    # start to end of the path
    indices <- match (unlist (edge_map [edges]), graph [[cols$edge_id]])
    indices <- indices [sort_transitions (graph [indices, ])]

    # Return the expanded path, in the same form as the input `path_c`
    if (vertices) {
        graph_expanded <- graph [indices, ]
        c (
            graph_expanded [[cols$from]],
            graph_expanded [[cols$to]] [nrow (graph_expanded)]
        )
    } else {
        indices
    }
}
