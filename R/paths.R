#' dodgr_paths
#'
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
#' @param wt_profile Name of weighting profile for street networks (one of foot,
#' horse, wheelchair, bicycle, moped, motorcycle, motorcar, goods, hgv, psv;
#' only used if `graph` is not provided, in which case a street network is
#' downloaded and correspondingly weighted).
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
#' `from[i]` and `to[j]` is `x [[i]] [[j]]`.
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
dodgr_paths <- function (graph, from, to, vertices = TRUE,
                         wt_profile = "bicycle", pairwise = FALSE,
                         heap = 'BHeap', quiet = TRUE)
{
    if (missing (graph) & (!missing (from) | !missing (to)))
        graph <- graph_from_pts (from, to, expand = 0.1,
                                 wt_profile = wt_profile, quiet = quiet)

    hps <- get_heap (heap, graph)
    heap <- hps$heap
    graph <- hps$graph

    gr_cols <- dodgr_graph_cols (graph)
    # cols are (edge_id, from, to, d, w, component, xfr, yfr, xto, yto)
    vert_map <- make_vert_map (graph, gr_cols)

    index_id <- get_index_id_cols (graph, gr_cols, vert_map, from)
    from_index <- index_id$index - 1 # 0-based
    if (!is.null (index_id$id))
        from_id <- index_id$id # can be null
    else
        from_id <- vert_map$vert

    index_id <- get_index_id_cols (graph, gr_cols, vert_map, to)
    to_index <- index_id$index - 1 # 0-based
    if (!is.null (index_id$id))
        to_id <- index_id$id # can be null
    else
        to_id <- vert_map$vert

    graph <- convert_graph (graph, gr_cols)

    if (!quiet)
        message ("Calculating shortest paths ... ", appendLF = FALSE)
    if (pairwise)
    {
        if (length (from_index) != length (to_index))
            stop ("pairwise paths require from and to to have same length")
        paths <- list ()
        for (i in seq (from_index))
        {
            paths [[i]] <- rcpp_get_paths (graph, vert_map,
                                           from_index [i], to_index [i],
                                           heap) [[1]]

        }
    } else
    {
        paths <- rcpp_get_paths (graph, vert_map, from_index, to_index, heap)
    }

    # convert 1-based indices back into vertex IDs. Note both paths that can not
    # be traced and single-step paths are returned from the above as NULL. The
    # former are retained as NULL, while the following converts the latter to
    # appropriate start-end vertices.
    paths <- lapply (paths, function (i)
                     lapply (i, function (j)
                             {
                                 if (is.null (j))
                                     return (j)
                                 vert_map$vert [j]
                             }  )   )


    # name path lists
    if (!is.null (from_id) & !is.null (to_id))
    {
        if (!pairwise)
        {
            for (i in seq (from_id))
                names (paths [[i]]) <- paste0 (from_id [i], "-", to_id)
        }
        names (paths) <- from_id
    }

    if (!vertices)
    {
        graph_verts <- paste0 ("f", graph$from, "t", graph$to)

        # convert vertex IDs to corresponding sequences of edge numbers
        paths <- lapply (paths, function (i)
                         lapply (i, function (j)
                                 if (length (j) > 1)
                                 {
                                     indx <- 2:length (j)
                                     pij <- paste0 ("f", j [indx - 1],
                                                    "t", j [indx])
                                     res <- match (pij, graph_verts)
                                     res <- res [which (!is.na (res))]
                                     return (if (length (res) == 0) NULL
                                             else res)
                                 } ))
    }

    return (paths)
}
