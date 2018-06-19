#' dodgr_dists
#'
#' Calculate matrix of pair-wise distances between points.
#'
#' @param graph \code{data.frame} or equivalent object representing the network
#' graph (see Details)
#' @param from Vector or matrix of points **from** which route distances are to
#' be calculated (see Details)
#' @param to Vector or matrix of points **to** which route distances are to be
#' calculated (see Details)
#' @param wt_profile Name of weighting profile for street networks (one of foot,
#' horse, wheelchair, bicycle, moped, motorcycle, motorcar, goods, hgv, psv).
#' @param expand Only when \code{graph} not given, the multiplicative factor by
#' which to expand the street network surrounding the points defined by
#' \code{from} and/or \code{to}.
#' @param heap Type of heap to use in priority queue. Options include
#' Fibonacci Heap (default; \code{FHeap}), Binary Heap (\code{BHeap}),
#' \code{Radix}, Trinomial Heap (\code{TriHeap}), Extended Trinomial Heap
#' (\code{TriHeapExt}, and 2-3 Heap (\code{Heap23}).
#' @param parallel If \code{TRUE}, perform routing calculation in parallel (see
#' details)
#' @param quiet If \code{FALSE}, display progress messages on screen.
#' @return square matrix of distances between nodes
#'
#' @note \code{graph} must minimally contain three columns of \code{from},
#' \code{to}, \code{dist}. If an additional column named \code{weight} or
#' \code{wt} is present, shortest paths are calculated according to values
#' specified in that column; otherwise according to \code{dist} values. Either
#' way, final distances between \code{from} and \code{to} points are calculated
#' according to values of \code{dist}. That is, paths between any pair of points
#' will be calculated according to the minimal total sum of \code{weight}
#' values (if present), while reported distances will be total sums of
#' \code{dist} values.
#'
#' The \code{from} and \code{to} columns of \code{graph} may be either single
#' columns of numeric or character values specifying the numbers or names of
#' graph vertices, or combinations to two columns specifying geographical
#' (longitude and latitude) coordinates. In the latter case, almost any sensible
#' combination of names will be accepted (for example, \code{fromx, fromy},
#' \code{from_x, from_y}, or \code{fr_lat, fr_lon}.)
#'
#' \code{from} and \code{to} values can be either two-column matrices of
#' equivalent of longitude and latitude coordinates, or else single columns
#' precisely matching node numbers or names given in \code{graph$from} or
#' \code{graph$to}. If \code{to} is missing, pairwise distances are calculated
#' between all points specified in \code{from}. If neither \code{from} nor
#' \code{to} are specified, pairwise distances are calculated between all nodes
#' in \code{graph}.
#'
#' Calculations in parallel (\code{parallel = TRUE}) ought very generally be
#' advantageous. For small graphs, Calculating distances in parallel is likely
#' to offer relatively little gain in speed, but increases from parallel
#' computation will generally markedly increase with increasing graph sizes.
#'
#' @export 
#' @examples
#' # A simple graph
#' graph <- data.frame (from = c ("A", "B", "B", "B", "C", "C", "D", "D"),
#'                      to = c ("B", "A", "C", "D", "B", "D", "C", "A"),
#'                      d = c (1, 2, 1, 3, 2, 1, 2, 1))
#' dodgr_dists (graph)
#'
#' # A larger example from the included \code{\link{hampi}} data.
#' graph <- weight_streetnet (hampi)
#' from <- sample (graph$from_id, size = 100)
#' to <- sample (graph$to_id, size = 50)
#' d <- dodgr_dists (graph, from = from, to = to)
#' # d is a 100-by-50 matrix of distances between \code{from} and \code{to}
#'
#' \dontrun{
#' # a more complex street network example, thanks to @chrijo; see
#' # https://github.com/ATFutures/dodgr/issues/47
#'
#' xy <- rbind (c (7.005994, 51.45774), # limbeckerplatz 1 essen germany
#'              c (7.012874, 51.45041)) # hauptbahnhof essen germany
#' xy <- data.frame (lon = xy [, 1], lat = xy [, 2])
#' essen <- dodgr_streetnet (pts = xy, expand = 0.2, quiet = FALSE)
#' graph <- weight_streetnet (essen, wt_profile = "foot")
#' d <- dodgr_dists (graph, from = xy, to = xy)
#' # First reason why this does not work is because the graph has multiple,
#' # disconnected components.
#' table (graph$component)
#' # reduce to largest connected component, which is always number 1
#' graph <- graph [which (graph$component == 1), ]
#' d <- dodgr_dists (graph, from = xy, to = xy)
#' # should work, but even then note that
#' table (essen$level)
#' # There are parts of the network on different building levels (because of
#' # shopping malls and the like). These may or may not be connected, so it may be
#' # necessary to filter out particular levels
#' levs <- paste0 (essen$level) # because sf data are factors
#' index <- which (! (levs == "-1" | levs == "1")) # for example
#' library (sf) # needed for following sub-select operation
#' essen <- essen [index, ]
#' graph <- weight_streetnet (essen, wt_profile = "foot")
#' d <- dodgr_dists (graph, from = xy, to = xy)
#' }
dodgr_dists <- function (graph, from, to, wt_profile = "bicycle", expand = 0,
                         heap = 'BHeap', parallel = TRUE, quiet = TRUE)
{
    if (missing (graph) & (!missing (from) | !missing (to)))
        graph <- graph_from_pts (from, to, expand = expand,
                                 wt_profile = wt_profile, quiet = quiet)

    hps <- get_heap (heap, graph)
    heap <- hps$heap
    graph <- hps$graph

    gr_cols <- dodgr_graph_cols (graph)
    vert_map <- make_vert_map (graph, gr_cols)

    index_id <- get_index_id_cols (graph, gr_cols, vert_map, from)
    from_index <- index_id$index - 1 # 0-based
    from_id <- index_id$id
    index_id <- get_index_id_cols (graph, gr_cols, vert_map, to)
    to_index <- index_id$index - 1 # 0-based
    to_id <- index_id$id

    graph <- convert_graph (graph, gr_cols)

    if (!quiet)
        message ("Calculating shortest paths ... ", appendLF = FALSE)

    if (parallel & heap == "TriHeapExt")
    {
        message ("Extended TriHeaps can not be calculated in parallel; ",
                 "reverting to serial calculation")
        parallel <- FALSE
    }

    if (parallel)
        d <- rcpp_get_sp_dists_par (graph, vert_map, from_index, to_index, heap)
    else
        d <- rcpp_get_sp_dists (graph, vert_map, from_index, to_index, heap)

    if (!is.null (from_id))
        rownames (d) <- from_id
    else
        rownames (d) <- vert_map$vert
    if (!is.null (to_id))
        colnames (d) <- to_id
    else
        colnames (d) <- vert_map$vert

    if (!quiet)
        message ("done.")

    return (d)
}

#' dodgr_distances
#'
#' Alias for \link{dodgr_dists}
#' @inherit dodgr_dists
#' @export
dodgr_distances <- function (graph, from, to, wt_profile = "bicycle", expand = 0,
                         heap = 'BHeap', parallel = TRUE, quiet = TRUE)
{
    dodgr_dists (graph, from, to, wt_profile = wt_profile, expand = expand,
                 heap = heap, parallel = parallel, quiet = quiet)
}

#' get_index_id_cols
#'
#' Get an index of \code{pts} matching \code{vert_map}, as well as the
#' corresonding names of those \code{pts}
#'
#' @return list of \code{index}, which is 0-based for C++, and corresponding
#' \code{id} values.
#' @noRd
get_index_id_cols <- function (graph, gr_cols, vert_map, pts)
{
    index <- -1
    id <- NULL
    if (!missing (pts))
    {
        index <- get_pts_index (graph, gr_cols, vert_map, pts)
        if (is.vector (pts) & length (pts == 2) & is.numeric (pts) &
            ( (any (grepl ("x", names (pts), ignore.case = TRUE)) &
             any (grepl ("y", names (pts), ignore.case = TRUE))) |
             (any (grepl ("lon", names (pts), ignore.case = TRUE) &
                   (any (grepl ("lat", names (pts), ignore.case = TRUE)))))))
            names (pts) <- NULL
        id <- get_id_cols (pts)
        if (is.null (id))
            id <- vert_map$vert [index] # from_index is 1-based
    }
    list (index = index, id = id)
}


#' get_id_cols
#'
#' Get the ID columns from a matrix or data.frame of from or to points
#' @param pts The \code{from} or \code{to} args passed to \code{dodgr_dists}
#' @return Character vector of names of points, if they exist in \code{pts}
#' @noRd
get_id_cols <- function (pts)
{
    ids <- NULL
    if (any (grepl ("id", colnames (pts), ignore.case = TRUE)))
    {
        nmc <- which (grepl ("id", colnames (pts)))
        if (is (pts, "data.frame"))
            ids <- pts [[nmc]]
        else if (is.matrix (pts))
            ids <- pts [, nmc, drop = TRUE]
    } else if (is.vector (pts) & !is.null (names (pts)))
        ids <- names (pts)
    else if (!is.null (rownames (pts)))
        ids <- rownames (pts)
    return (ids)
}

#' make_vert_map
#'
#' Map unique vertex names to sequential numbers in matrix
#' @noRd
make_vert_map <- function (graph, gr_cols)
{
    # gr_cols are (edge_id, from, to, d, w, component, xfr, yfr, xto, yto)
    verts <- c (paste0 (graph [[gr_cols [2] ]]),
                paste0 (graph [[gr_cols [3] ]]))
    indx <- which (!duplicated (verts))
    # Note id has to be 0-indexed:
    data.frame (vert = paste0 (verts [indx]), id = seq (indx) - 1,
                stringsAsFactors = FALSE)
}

#' get_pts_index
#'
#' Convert \code{from} or \code{to} args of \code{dodgr_dists} to indices into
#' \code{vert_map}
#'
#' @param vert_map Two-column \code{data.frame} of unique vertices and
#' corresponding IDs, obtained from \code{make_vert_map}
#' @param xy List of x (longitude) and y (latitude) coordinates of all vertices
#' in \code{vert_map}
#' @param pts Either a vector of names, or a matrix or \code{data.frame} of
#' arbitrary geographical coordinates for which to get index into vertices of
#' graph.
#'
#' @noRd
get_pts_index <- function (graph, gr_cols, vert_map, pts)
{
    if (!(is.matrix (pts) | is.data.frame (pts)))
    {
        if (!is.numeric (pts))
            pts <- matrix (pts, ncol = 1)
        else
        {
            nms <- names (pts)
            if (is.null (nms))
                nms <- c ("x", "y")
            pts <- matrix (pts, nrow = 1) # vector of (x,y) vals
            colnames (pts) <- nms
        }
    }

    if (ncol (pts) == 1)
    {
        pts <- pts [, 1]
        if (!is.numeric (pts))
        {
            indx <- match (pts, vert_map$vert)
            if (any (is.na (indx)))
                stop (paste0 ("from/to are not numeric yet can not be",
                              " matched onto graph vertices"))
            pts <- indx
        }
        if (any (pts < 1 | pts > nrow (vert_map)))
            stop (paste0 ("points exceed numbers of vertices"))
    } else
    {
        nms <- names (pts)
        if (is.null (nms))
            nms <- colnames (pts)
        ix <- which (grepl ("x", nms, ignore.case = TRUE) |
                     grepl ("lon", nms, ignore.case = TRUE))
        iy <- which (grepl ("y", nms, ignore.case = TRUE) |
                     grepl ("lat", nms, ignore.case = TRUE))
        if (length (ix) != 1 | length (iy) != 1)
            stop (paste0 ("Unable to determine geographical ",
                          "coordinates in pts"))

        # gr_cols are (edge_id, from, to, d, w, component, xfr, yfr, xto, yto)
        if (any (is.na (gr_cols [7:10])))
            stop (paste0 ("Cannot determine geographical coordinates ",
                          "against which to match pts"))

        names (pts) [ix] <- "x"
        names (pts) [iy] <- "y"

        # Result of rcpp_points_index is 0-indexed for C++
        pts <- rcpp_points_index (dodgr_vertices (graph), pts) + 1
        # xy has same order as vert_map
    }

    pts
}

#' get_heap
#'
#' Match the heap arg and convert graph is necessary (for Radix)
#' @param heap Name of heap as passed to \code{dodgr_dists}
#' @param graph \code{data.frame} of graph edges
#' @return List of matched heap arg and potentially converted graph
#' @noRd
get_heap <- function (heap, graph)
{
    heaps <- c ("FHeap", "BHeap", "Radix", "TriHeap", "TriHeapExt", "Heap23")
    heap <- match.arg (arg = heap, choices = heaps)
    if (heap == "Radix")
    {
        dfr <- min (abs (c (graph$d %% 1, graph$d %% 1 - 1)))
        if (dfr > 1e-6)
        {
            message (paste0 ("RadixHeap can only be implemented for ",
                             "integer weights;\nall weights will now be ",
                             "rounded"))
            graph$d <- round (graph$d)
            graph$d_weighted <- round (graph$d_weighted)
        }
    }

    list (heap = heap, graph = graph)
}

#' graph_from_pts
#'
#' Download a street network when \code{graph} not passed to \code{dodgr_dists},
#' by using the lists of from and to points.
#' @param from Arg passed to \code{dodgr_dists}
#' @param to Arg passed to \code{dodgr_dists}
#' @param expand Factor by which street network is to be expanded beyond range
#' of \code{from} and \code{to} points.
#' @return Converted graph as \code{data.frame}
#' @noRd
graph_from_pts <- function (from, to, expand = 0.1, wt_profile = "bicycle",
                            quiet = TRUE)
{
    if (!quiet)
        message (paste0 ("No graph submitted to dodgr_dists; ",
                         "downloading street network ... "),
                 appendLF = FALSE)

    pts <- NULL
    if (!missing (from))
        pts  <- from
    if (!missing (to))
        pts <- rbind (pts, to)
    pts <- pts [which (!duplicated (pts)), ]
    graph <- dodgr_streetnet (pts = pts, expand = expand) %>%
        weight_streetnet (wt_profile = wt_profile)

    if (!quiet)
        message ("done")

    return (graph)
}
