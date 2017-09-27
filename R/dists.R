#' dodgr_dists
#'
#' @param graph \code{data.frame} or equivalent object representing the network
#' graph (see Details)
#' @param from Vector or matrix of points **from** which route distances are to
#' be calculated (see Details)
#' @param to Vector or matrix of points **to** which route distances are to be
#' calculated (see Details)
#' @param wt_profile Name of weighting profile for street networks (one of foot,
#' horse, wheelchair, bicycle, moped, motorcycle, motorcar, goods, hgv, psv).
#' @param heap Type of heap to use in priority queue. Options include
#' Fibonacci Heap (default; \code{FHeap}), Binary Heap (\code{BHeap}),
#' \code{Radix}, Trinomial Heap (\code{TriHeap}), Extended Trinomial Heap
#' (\code{TriHeapExt}, and 2-3 Heap (\code{Heap23}).
#' @param quiet If \code{FALSE}, display progress messages on screen.
#' @return square matrix of distances between nodes
#'
#' @note \code{graph} must minimially contain four columns of \code{from},
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
#' \code{to} are specified, pairwise distances are calcualted between all nodes
#' in \code{graph}.
#'
#' @export
#' @examples
#' graph <- weight_streetnet (hampi)
#' from <- sample (graph$from_id, size = 100)
#' to <- sample (graph$to_id, size = 50)
#' d <- dodgr_dists (graph, from = from, to = to)
#' # d is a 100-by-50 matrix of distances between \code{from} and \code{to}
dodgr_dists <- function (graph, from, to, wt_profile = "bicycle",
                         heap = 'BHeap', quiet = TRUE)
{
    if (missing (graph) & !missing (from))
    {
        if (!quiet)
            message (paste0 ("No graph submitted to dodgr_dists; ",
                             "downloading street network ... "),
                     appendLF = FALSE)
        graph <- dodgr_streetnet (pts = from, expand = 0.1) %>%
            weight_streetnet (wt_profile = wt_profile)
    }

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

    if (!quiet)
        message ("done\nConverting network to dodgr graph ... ",
                 appendLF = FALSE)
    graph <- dodgr_convert_graph (graph, components = FALSE)
    xy <- graph$xy
    graph <- graph$graph
    vert_map <- make_vert_map (graph)

    if (missing (from))
        from <- -1
    else
        from <- get_pts_index (vert_map, xy, from) # 0-indexed

    if (missing (to))
        to <- -1
    else
        to <- get_pts_index (vert_map, xy, to)

    if (!quiet)
        message ("done\nCalculating shortest paths ... ", appendLF = FALSE)
    d <- rcpp_get_sp (graph, vert_map, from, to, heap)
    if (any (from < 0))
        from <- seq (nrow (vert_map)) - 1
    if (any (to < 0))
        to <- seq (nrow (vert_map)) - 1
    rownames (d) <- vert_map$vert [from + 1] # coz from is 0-indexed
    colnames (d) <- vert_map$vert [to + 1] # ditto
    if (!quiet)
        message ("done.")

    return (d)
}

#' make_vert_map
#'
#' Map unique vertex names to sequential numbers in matrix
#' @noRd
make_vert_map <- function (graph)
{
    verts <- c (graph$from, graph$to)
    indx <- which (!duplicated (verts))
    # Note id has to be 0-indexed:
    data.frame (vert = verts [indx], id = seq (indx) - 1,
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
#' @param pts Matrix or \code{data.frame} of arbitrary geographical coordinates
#' for which to get index into vertices of graph.
#'
#' @noRd
get_pts_index <- function (vert_map, xy, pts)
{
    if (!(is.matrix (pts) | is.data.frame (pts)))
        pts <- matrix (pts, ncol = 1)

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
        if (is.null (xy))
            stop (paste0 ("xy has no geographical coordinates ",
                          "against which to match pts"))

        pts <- rcpp_points_index (xy, pts)
        # xy has same order as vert_map
    }

    pts - 1 # 0-indexed for C++
}
