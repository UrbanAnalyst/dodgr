# flatten lists of lists to single list
flatten_list <- function (x)
{
    x2 <- list ()
    for (i in x)
        x2 <- c (x2, i)
    return (x2)
}

#' dodgr_fundamental_cycles
#'
#' Calculate fundamental cycles in a graph.
#'
#' @param graph `data.frame` or equivalent object representing the contracted
#' network graph (see Details).
#' @param vertices `data.frame` returned from \link{dodgr_vertices}`(graph)`. Will
#' be calculated if not provided, but it's quicker to pass this if it has
#' already been calculated.
#' @param graph_max_size Maximum size submitted to the internal C++ routines as
#' a single chunk. Warning: Increasing this may lead to computer meltdown!
#' @param expand For large graphs which must be broken into chunks, this factor
#' determines the relative overlap between chunks to ensure all cycles are
#' captured. (This value should only need to be modified in special cases.)
#' @return List of cycle paths, in terms of vertex IDs in `graph` and, for
#' spatial graphs, the corresponding coordinates.
#'
#' @note Calculation of fundamental cycles is VERY computationally demanding,
#' and this function should only be executed on CONTRACTED graphs (that is,
#' graphs returned from \link{dodgr_contract_graph}), and even than may take a
#' long time to execute. Results for full graphs can be obtained with the
#' function \link{dodgr_full_cycles}. The computational complexity can also not
#' be calculated in advance, and so the parameter `graph_max_size` will lead to
#' graphs larger than that (measured in numbers of edges) being cut into smaller
#' parts. (Note that that is only possible for spatial graphs, meaning that it
#' is not at all possible to apply this function to large, non-spatial graphs.)
#' Each of these smaller parts will be expanded by the specified amount
#' (`expand`), and cycles found within. The final result is obtained by
#' aggregating all of these cycles and removing any repeated ones arising due to
#' overlap in the expanded portions. Finally, note that this procedure of
#' cutting graphs into smaller, computationally manageable sub-graphs provides
#' only an approximation and may not yield all fundamental cycles.
#'
#' @examples 
#' net <- weight_streetnet (hampi)
#' graph <- dodgr_contract_graph (net)
#' verts <- dodgr_vertices (graph)
#' cyc <- dodgr_fundamental_cycles (graph, verts)
#' @export 
dodgr_fundamental_cycles <- function (graph, vertices = NULL,
                                      graph_max_size = 10000, expand = 0.05)
{
    if (missing (graph))
        stop ("graph must be provided")
    if (!inherits (graph, "data.frame"))
        stop ("graph must be a data.frame object")

    if (is.null (vertices))
        vertices <- dodgr_vertices (graph)
    if (!"flow" %in% names (graph)) # makes no difference
        graph$flow <- 1
    graph <- merge_directed_flows (graph) # uses fast C++ routines
    graph$flow <- NULL
    bb <- get_graph_bb (graph)
    if (nrow (graph) <= graph_max_size)
    {
        bb_indices <- list (bb)
    } else
    {
        ndivs <- get_ndivs (graph, graph_max_size)
        bb_list <- get_bb_list (bb, ndivs, expand = expand)
        bb_data <- subdivide_bb (graph, bb_list, graph_max_size, expand)
        bb_indices <- bb_data$bb_indices
    }

    graphc <- convert_graph (graph, dodgr_graph_cols (graph))

    if (length (bb_indices) == 1)
    {
        res <- rcpp_fundamental_cycles (graphc, verts = vertices)
    } else
    {
        message ("Now computing fundamental cycles by breaking graph with ",
                 nrow (graphc), " edges into ", length (bb_indices),
                 " components ...")
        pb <- utils::txtProgressBar (style = 3)
        res <- list ()
        for (i in seq (bb_indices))
        {
            graphi <- graphc [bb_indices [[i]], ]
            verti <- dodgr_vertices (graphi)
            res [[i]] <- rcpp_fundamental_cycles (graphi, verts = verti)
            utils::setTxtProgressBar (pb, i / length (bb_indices))
        }
        close (pb)

        # each element of res is a list, so flatten these:
        res <- flatten_list (res)
        # These hash each and remove any duplicated ones:
        dig <- unlist (lapply (res, digest::digest))
        res <- res [which (!duplicated (dig))]
    }

    if (is_graph_spatial (graph))
    {
        if (length (bb_indices) > 1)
            message ("Generating spatial coordinates of polygons ",
                     "(this should be fairly quick ...)")
        res <- lapply (res, function (i) {
                           index <- match (i, vertices$id)
                           data.frame (id = i,
                                       x = vertices$x [index],
                                       y = vertices$y [index],
                                       stringsAsFactors = FALSE)
                                })
    }
    return (res)
}

# ********** FUNCTIONS TO BREAK SPATIAL GRAPHS INTO SUB-COMPONENTS **********

# Initial estimate of how many divisions needed
get_ndivs <- function (graph, graph_max_size)
{
    ndivs <- ceiling (nrow (graph) / graph_max_size)
    ceiling (sqrt (ndivs)) # num of grid rows and cols
}

# get boundingn box of graph
get_graph_bb <- function (graph)
{
    gr_cols <- dodgr_graph_cols (graph)
    from_lon <- graph [, gr_cols$xfr]
    from_lat <- graph [, gr_cols$yfr]
    to_lon <- graph [, gr_cols$xto]
    to_lat <- graph [, gr_cols$yto]
    apply (cbind (c (from_lon, to_lon), c (from_lat, to_lat)), 2, range)
}


# divide a bb into rectangular grid of sub-boxes and return list of
# corresponding bboxes
get_bb_list <- function (bb, ndivs, expand = 0.05)
{
    # divide one column of bb: either lons or lats
    divide_bb_vec <- function (bb, ndivs, colnum = 2, expand)
    {
        bb <- c (bb [1, colnum], vapply (seq (ndivs), function (i)
                                         bb [1, colnum] + i / ndivs *
                                             diff (bb [, colnum]),
                                         numeric (1)))
        bb <- cbind (bb [1:ndivs], bb [2:(ndivs + 1)])
        t (apply (bb, 1, function (i) 
                  mean (i) + c (-0.5 - expand, 0.5 + expand) * diff (i)))
    }
    bb_lons <- divide_bb_vec (bb, ndivs, colnum = 1, expand = expand)
    bb_lats <- divide_bb_vec (bb, ndivs, colnum = 2, expand = expand)
    bb_list <- list ()
    for (i in seq (ndivs))
        for (j in seq (ndivs))
        {
            bb_list [[length (bb_list) + 1]] <-
                cbind (bb_lons [i, ], bb_lats [j, ])
        }
    return (bb_list)
}

# get indices into graph of edges lying within each bit of a bb_list
get_bb_indices <- function (graph, bb_list)
{
    res <- list ()
    for (i in seq (bb_list))
    {
        res [[i]] <- which (graph$from_lon > bb_list [[i]] [1, 1] &
                            graph$from_lon < bb_list [[i]] [2, 1] &
                            graph$from_lat > bb_list [[i]] [1, 2] &
                            graph$from_lat < bb_list [[i]] [2, 2] &
                            graph$to_lon > bb_list [[i]] [1, 1] &
                            graph$to_lon < bb_list [[i]] [2, 1] &
                            graph$to_lat > bb_list [[i]] [1, 2] &
                            graph$to_lat < bb_list [[i]] [2, 2])
    }
    return (res)
}

# Rectangularly subdivide any components of bb_list that are > graph_max_size
# into 4 sub-components.
subdivide_bb <- function (graph, bb_list, graph_max_size, expand)
{
    bb_indices <- get_bb_indices (graph, bb_list)
    lens <- unlist (lapply (bb_indices, length))
    while (any (lens > graph_max_size))
    {
        indx <- which (lens > graph_max_size)
        bbs <- bb_list [indx]
        bb_list [indx] <- NULL
        bb_indices [indx] <- NULL
        for (i in bbs)
        {
            bb_list <- c (bb_list, get_bb_list (i, ndivs = 2, expand = expand))
            bb_indices <- get_bb_indices (graph, bb_list)
            lens <- unlist (lapply (bb_indices, length))
        }
        bb_list [which (lens == 0)] <- NULL
        bb_indices [which (lens == 0)] <- NULL
    }
    list (bb_list = bb_list, bb_indices = bb_indices)
}

#' dodgr_full_cycles
#'
#' Calculate fundamental cycles on a FULL (that is, non-contracted) graph.
#' @inheritParams dodgr_fundamental_cycles
#' @note This function converts the `graph` to its contracted form, calculates
#' the fundamental cycles on that version, and then expands these cycles back
#' onto the original graph. This is far more computationally efficient than
#' calculating fundamental cycles on a full (non-contracted) graph.
#'
#' @examples 
#' net <- weight_streetnet (hampi)
#' graph <- dodgr_contract_graph (net)
#' cyc1 <- dodgr_fundamental_cycles (graph)
#' cyc2 <- dodgr_full_cycles (net)
#' # cyc2 has same number of cycles, but each one is generally longer, through
#' # including all points intermediate to junctions; cyc1 has cycles composed of
#' # junction points only.
#' @export
dodgr_full_cycles <- function (graph, graph_max_size = 10000, expand = 0.05)
{
    graph$flow <- 1
    graph <- merge_directed_flows (graph)
    graph$flow <- NULL
    #graph <- graph [graph$component == 1, ]
    graphc <- dodgr_contract_graph (graph)
    v <- dodgr_vertices (graphc)

    #edge_map <- get_edge_map (graphc) # TODO: Implement this
    hashc <- get_hash (graphc, hash = FALSE)
    fname_c <- file.path (tempdir (), paste0 ("dodgr_edge_map_", hashc, ".Rds"))
    if (!file.exists (fname_c))
        stop ("something unexpected went wrong extracting the edge map") # nocov
    edge_map <- readRDS (fname_c)

    x <- dodgr_fundamental_cycles (graphc,
                                   vertices = v,
                                   graph_max_size = graph_max_size,
                                   expand = expand)

    from_to <- paste0 (graphc$from_id, "-", graphc$to_id)
    to_from <- paste0 (graphc$to_id, "-", graphc$from_id)
    ids <- lapply (x, function (i) {
           idpairs <- paste0 (i$id [-length (i$id)], "-", i$id [-1])
           # Get the edge pairs that match the idpairs, whether as from->to or
           # to->from
           edges_c1 <- graphc$edge_id [match (idpairs, from_to)]
           edges_c2 <- graphc$edge_id [match (idpairs, to_from)]
           edges_c <- apply (rbind (edges_c1, edges_c2), 2, function (i)
                             i [which (!is.na (i))] [1])
           edges_new <- lapply (as.list (edges_c), function (j) {
                    if (j %in% edge_map$edge_new)
                    {
                        index <- which (edge_map$edge_new %in% j)
                        index <- as.numeric (edge_map$edge_old [index])
                        j <- match (index, graph$edge_id)
                        if ((graph$from_id [j [1] ] ==
                             graph$to_id [utils::tail (j, 1)]) ||
                            (graph$from_id [j [1] ] == graph$to_id [j [2] ]))
                            j <- rev (j)
                    } else
                        j <- match (j, graph$edge_id)
                    return (as.numeric (j))
                    }) # end edges_new lapply
           unlist (edges_new)
            }) # end ids lapply
    # ids at that point is a sequences of indices into graph. This is then
    # converted to a sequence of vertex IDs, through just adding the last vertex
    # of the sequence on to close the polygon
    gr_cols <- dodgr_graph_cols (graph) # (from, to) = [, 2:3]
    res <- lapply (ids, function (i)
                   graph [c (i, i [1]), gr_cols$from] )

    if (is_graph_spatial (graph))
    {
        vertices <- dodgr_vertices (graph)
        res <- lapply (res, function (i) {
                           index <- match (i, vertices$id)
                           data.frame (id = i,
                                       x = vertices$x [index],
                                       y = vertices$y [index],
                                       stringsAsFactors = FALSE)
                                })
    }
    return (res)
}



#' dodgr_sflines_to_poly
#'
#' Convert \pkg{sf} `LINESTRING` objects to `POLYGON` objects representing all
#' fundamental cycles within the `LINESTRING` objects.
#'
#' @inheritParams dodgr_fundamental_cycles
#' @param sflines An \pkg{sf} `LINESTRING` object representing a network.
#' @return An `sf::sfc` collection of `POLYGON` objects.
#' @export
dodgr_sflines_to_poly <- function (sflines, graph_max_size = 10000,
                                   expand = 0.05)
{
    if (!(methods::is (sflines, "sf") | methods::is (sflines, "sf")))
        stop ("lines must be an object of class 'sf' or 'sfc'")
    if (!methods::is (sflines [[attr (sflines, "sf_column")]], "sfc_LINESTRING"))
        stop ("lines must be an 'sfc_LINESTRING' object")

    graph <- weight_streetnet (sflines, wt_profile = 1)
    # Different graph components need to be analysed seperately, and an
    # arbitrary decision is made here to only consider components with > 100
    # edges - smaller ones are unlikely to have any cycles.
    comps <- table (graph$component)
    comps <- as.numeric (names (comps) [which (comps > 100)])
    x <- list ()
    for (i in seq (comps))
    {
        graphi <- graph [graph$component == comps [i], ]
        graphi$edge_id <- seq (nrow (graphi))
        attr (graphi, "hash") <- NULL
        x [[i]] <-  dodgr_full_cycles (graphi,
                                       graph_max_size = graph_max_size,
                                       expand = expand)
    }
    x <- flatten_list (x)
    polys_to_sfc (x, sflines)
}

# convert list of polygon coordinates to `sf::sfc` collection
# code adapted from osmdata/tests/testthat/test-sf-construction.R
polys_to_sfc <- function (x, sflines)
{
    g <- sflines [[attr (sflines, "sf_column")]]
    crs <- attr (g, "crs")

    xy <- do.call (rbind, x)
    xvals <- xy [, 2]
    yvals <- xy [, 3]
    bb <- structure (rep (NA_real_, 4),
                     names = c("xmin", "ymin", "xmax", "ymax"))
    bb [1:4] <- c (min (xvals), min (yvals), max (xvals), max (yvals))
    class (bb) <- "bbox"
    attr (bb, "crs") <- crs

    x <- lapply (x, function (i) {
                     res <- as.matrix (i [, 2:3])
                     colnames (res) <- NULL
                     structure (list (res), class = c ("XY", "POLYGON", "sfg")) })
    attr (x, "n_empty") <- 0
    attr(x, "precision") <- 0.0
    class(x) <- c ("sfc_POLYGON", "sfc")
    attr(x, "bbox") <- bb
    attr(x, "crs") <- crs
    x
}
