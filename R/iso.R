#' dodgr_isodists
#'
#' Isodistances from specified points. Function is fully vectorized to calculate
#' accept vectors of central points and vectors of isodistance thresholds.
#'
#' @param graph `data.frame` or equivalent object representing the network
#' graph (see Notes)
#' @param from Vector or matrix of points **from** which isodistances are to
#' be calculated (see Notes)
#' @param dlim Desired limits of isodistances in metres.
#' @param heap Type of heap to use in priority queue. Options include
#' Fibonacci Heap (default; `FHeap`), Binary Heap (`BHeap`),
#' `Radix`, Trinomial Heap (`TriHeap`), Extended Trinomial Heap
#' (`TriHeapExt`, and 2-3 Heap (`Heap23`).
#' @return A single `data.frame` of isodistances as points sorted anticlockwise
#' around each origin (`from`) point, with columns denoting the `from` points
#' and `dlim` value(s). The isodistances are given as `id` values of the series
#' of points from each from point at the specified isodistance.
#'
#' @export 
#' @examples
#' graph <- weight_streetnet (hampi)
#' from <- sample (graph$from_id, size = 100)
#' dlim <- c (1, 2, 5, 10, 20) * 100
#' d <- dodgr_isodists (graph, from = from, dlim)
dodgr_isodists <- function (graph, from = NULL, dlim = NULL, heap = 'BHeap')
{
    if (is.null (dlim))
        stop ("dlim must be specified")
    if (!is.numeric (dlim))
        stop ("dlim must be numeric")

    v <- dodgr_vertices (graph)
    graph <- tbl_to_df (graph)

    hps <- get_heap (heap, graph)
    heap <- hps$heap
    graph <- hps$graph

    gr_cols <- dodgr_graph_cols (graph)
    if (is.na (gr_cols$from) | is.na (gr_cols$to))
    {
        scols <- find_spatial_cols (graph)
        graph$from_id <- scols$xy_id$xy_fr_id
        graph$to_id <- scols$xy_id$xy_to_id
        gr_cols <- dodgr_graph_cols (graph)
    }
    vert_map <- make_vert_map (graph, gr_cols, FALSE)

    tp <- attr (graph, "turn_penalty")
    tp <- ifelse (is.null (tp), 0, tp)
    if (is (graph, "dodgr_streetnet_sc") & tp > 0)
    {
        if (!is.null (from))
        {
            from <- nodes_arg_to_pts (from, graph)
            from <- remap_verts_with_turn_penalty (graph, from, from = TRUE)
        }
    }

    from_index <- get_to_from_index (graph, vert_map, gr_cols, from)

    graph <- convert_graph (graph, gr_cols)

    d <- rcpp_get_iso (graph, vert_map, from_index$index, dlim, heap)

    if (!is.null (from_index$id))
        rownames (d) <- from_index$id
    else
        rownames (d) <- vert_map$vert
    colnames (d) <- vert_map$vert

    return (dmat_to_pts (d, from, v, dlim))
}

# convert distance matrix with values equal to various isodistances into list of
# lists of points ordered around the central points
dmat_to_pts <- function (d, from, v, dlim)
{
    pt_names <- colnames (d)
    pts <- list ()
    for (i in seq (nrow (d)))
    {
        o <- v [match (from [i], v$id), ]
        pts [[i]] <- lapply (dlim, function (j) {
                                 res <- pt_names [which (d [i, ] == j)]
                                 res <- v [match (res, v$id), ]
                                 if (nrow (res) > 0)
                                 {
                                     res$from <- o$id
                                     res <- order_points (res, o)
                                     res$dlim <- j
                                     res <- res [, c ("from", "dlim",
                                                      "id", "x", "y")]
                                 }
                                 return (res)
            })
        names (pts [[i]]) <- paste (dlim)
    }
    names (pts) <- rownames (d)

    # flatten lists
    pts <- do.call (rbind, lapply (pts, function (i) do.call (rbind, i)))
    rownames (pts) <- NULL

    return (pts)
}

# order points around circle
order_points <- function (pts, origin)
{
    dx <- pts$x - origin$x
    dy <- pts$y - origin$y
    theta <- rep (NA, nrow (pts))

    index <- which (dx > 0 & dy >= 0)
    theta [index] <- atan (dy [index] / dx [index])
    index <- which (dx > 0 & dy < 0)
    theta [index] <- 2 * pi + atan (dy [index] / dx [index])
    index <- which (dx < 0)
    theta [index] <- pi + atan (dy [index] / dx [index])
    index <- which (dx == 0 & dy >= 0)
    theta [index] <- 0
    index <- which (dx == 0 & dy < 0)
    theta [index] <- 3 * pi / 2

    pts <- pts [order (theta), c ("from", "id", "x", "y")]
    rbind (pts, pts [1, ])
}
