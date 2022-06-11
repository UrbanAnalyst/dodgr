context ("dodgr_dists")

test_all <- (identical (Sys.getenv ("MPADGE_LOCAL"), "true") |
    identical (Sys.getenv ("GITHUB_WORKFLOW"), "test-coverage"))

if (!test_all) {
    RcppParallel::setThreadOptions (numThreads = 2)
}

test_that ("dists", {
    expect_silent (graph <- weight_streetnet (hampi))
    nf <- 100
    nt <- 50
    set.seed (1)
    from <- sample (graph$from_id, size = nf)
    to <- sample (graph$to_id, size = nt)
    expect_silent (d <- dodgr_distances (graph, from = from, to = to))
    expect_equal (nrow (d), nf)
    expect_equal (ncol (d), nt)
    expect_true (all (d [!is.na (d)] >= 0))
    expect_message (
        d2 <- dodgr_dists (graph, from = from, to = to, quiet = FALSE),
        "Calculating shortest paths ..."
    )
    expect_identical (d, d2)

    from [1] <- "not_a_vertex_id"
    expect_error (
        d <- dodgr_distances (graph, from = from, to = to),
        "from/to are not numeric yet can not be matched onto graph vertices"
    )

    # dists from coordinates:
    bb <- attr (hampi$geometry, "bbox")
    fromx <- bb [1] + (bb [3] - bb [1]) * runif (nf)
    fromy <- bb [2] + (bb [4] - bb [2]) * runif (nf)
    tox <- bb [1] + (bb [3] - bb [1]) * runif (nt)
    toy <- bb [2] + (bb [4] - bb [2]) * runif (nt)
    from <- data.frame (x = fromx, y = fromy, id = paste0 ("f", 1:nf))
    to <- data.frame (x = tox, y = toy, id = paste0 ("t", 1:nt))
    d <- dodgr_dists (graph, from = from, to = to)
    expect_equal (nrow (d), nf)
    expect_equal (ncol (d), nt)
    expect_true (all (d [!is.na (d)] >= 0))

    # from as vector
    from <- c (as.numeric (from [1, 1:2]), 1)
    names (from) <- c ("x", "y", "id")
    expect_silent (d <- dodgr_dists (graph, from = from, to = to))
    from <- as.numeric (from [1:2])
    expect_silent (d <- dodgr_dists (graph, from = from, to = to))

    # from as matrix
    from <- cbind (fromx, fromy, 1:nf)
    colnames (from) <- c ("x", "y", "id")
    d <- dodgr_dists (graph, from = from, to = to)
    from <- from [, 1:2]
    expect_silent (d <- dodgr_dists (graph, from = from, to = to))
    rownames (from) <- 1:nf
    expect_silent (d <- dodgr_dists (graph, from = from, to = to))

    from <- data.frame (x = fromx, y = fromy, id = paste0 ("f", 1:nf))
    to <- data.frame (x = tox, y = toy, id = paste0 ("t", 1:nt))
    from <- cbind (from, "x2" = from$x)
    expect_error (
        d <- dodgr_dists (graph, from = from, to = to),
        "Unable to determine geographical coordinates in from/to"
    )

    # from <- sample (graph$from_id, size = nf)
    # to <- sample (graph$to_id, size = nt)
    # graph0 <- graph
    # graph <- graph0
    # graph$from_id <- graph$to_id <- NULL
    # find_spatial_cols (graph)
})

test_that ("dists-pairwise", {
    expect_silent (graph <- weight_streetnet (hampi))
    n <- 50
    set.seed (1)
    from <- sample (graph$from_id, size = n)
    to <- sample (graph$to_id, size = n)
    expect_silent (d <- dodgr_distances (graph, from = from, to = to))
    expect_equal (dim (d), c (n, n))
    expect_silent (d <- dodgr_distances (graph,
        from = from, to = to,
        pairwise = TRUE
    ))
    expect_equal (dim (d), c (50, 1))
})

test_that ("times", {
    graph <- weight_streetnet (hampi)
    nf <- 100
    nt <- 50
    set.seed (1)
    from <- sample (graph$from_id, size = nf)
    to <- sample (graph$to_id, size = nt)
    expect_silent (d0 <- dodgr_dists (graph,
        from = from, to = to,
        shortest = TRUE
    )) # default
    expect_silent (d1 <- dodgr_dists (graph,
        from = from, to = to,
        shortest = FALSE
    ))
    expect_silent (t0 <- dodgr_times (graph,
        from = from, to = to,
        shortest = TRUE
    ))
    expect_silent (t1 <- dodgr_times (graph,
        from = from, to = to,
        shortest = FALSE
    )) # default

    expect_true (!identical (d0, d1))
    expect_true (!identical (d0, t0))
    expect_true (!identical (d0, t1))
    expect_true (!identical (d1, t0))
    expect_true (!identical (d1, t1))
    expect_true (!identical (t0, t1))

    # times are just dists using different columns:
    grapht <- graph
    grapht$d <- grapht$time
    grapht$d_weighted <- grapht$time_weighted
    grapht$time_weighted <- NULL
    expect_silent (t0 <- dodgr_times (grapht,
        from = from, to = to,
        shortest = TRUE
    ))
    expect_error (
        t0 <- dodgr_times (grapht,
            from = from, to = to,
            shortest = FALSE
        ),
        "Graph does not contain a weighted time column"
    )
    expect_error (
        d0 <- dodgr_dists (grapht,
            from = from, to = to,
            shortest = FALSE
        ),
        "Graph does not contain a weighted time column"
    )
    expect_silent (d1 <- dodgr_dists (grapht,
        from = from, to = to,
        shortest = TRUE
    ))
    expect_identical (t0, d1)

    grapht$time <- NULL
    expect_error (
        t0 <- dodgr_times (grapht,
            from = from, to = to,
            shortest = TRUE
        ),
        "graph has no time column"
    )
    expect_error (
        t0 <- dodgr_times (grapht,
            from = from, to = to,
            shortest = FALSE
        ),
        "graph has no time column"
    )
    expect_error (
        t0 <- dodgr_dists (grapht,
            from = from, to = to,
            shortest = FALSE
        ),
        "Graph does not contain a weighted time column"
    )
    expect_silent (t2 <- dodgr_dists (grapht,
        from = from, to = to,
        shortest = TRUE
    ))

    expect_identical (t2, t0)
})

test_that ("all dists", {
    graph <- weight_streetnet (hampi)
    graph <- graph [graph$component == 2, ]
    expect_silent (d <- dodgr_dists (graph))
    v <- dodgr_vertices (graph)
    expect_equal (nrow (d), ncol (d))
    expect_equal (nrow (d), nrow (v))
})

test_that ("to-from-cols", {
    graph <- weight_streetnet (hampi)
    nf <- 100
    nt <- 50
    set.seed (1)
    v <- dodgr_vertices (graph)
    from <- sample (v$id, size = nf)
    to <- sample (v$id, size = nt)
    expect_silent (d0 <- dodgr_dists (graph, from = from, to = to))

    fromi <- match (from, v$id)
    toi <- match (to, v$id)
    expect_silent (d1 <- dodgr_dists (graph, from = fromi, to = toi))
    expect_identical (d0, d1)

    fromm <- as.matrix (fromi, ncol = 1)
    tom <- as.matrix (toi, ncol = 1)
    expect_silent (d2 <- dodgr_dists (graph, from = fromm, to = tom))
    expect_identical (d0, d2)

    fromm [1] <- -1
    expect_error (
        d2 <- dodgr_dists (graph, from = fromm, to = tom),
        "points exceed numbers of vertices"
    )

    fromf <- as.factor (fromi)
    expect_error (
        d2 <- dodgr_dists (graph, from = fromf, to = toi),
        paste0 (
            "routing points are of unknown form; ",
            "must be either character, matrix, or integer"
        )
    )

    from <- sample (nrow (v), size = nf)
    to <- sample (nrow (v), size = nt)
    to [1] <- nrow (v) + 1L
    expect_error (
        d2 <- dodgr_dists (graph, from = from, to = to),
        "Unable to match all routing points to graph vertices"
    )

    to <- sample (nrow (v), size = nt)
    graph$from_id <- graph$from_lon <- NULL
    expect_error (
        d3 <- dodgr_dists (graph, from = from, to = to),
        "Graph appears to be spatial yet unable to extract coordinates"
    )
})

test_that ("dists with no edge ids", {
    graph <- weight_streetnet (hampi)
    nf <- 100
    nt <- 50
    set.seed (1)
    from <- sample (graph$from_id, size = nf)
    to <- sample (graph$to_id, size = nt)
    expect_silent (d0 <- dodgr_distances (graph, from = from, to = to))

    # from/to as coordinates only:
    v <- dodgr_vertices (graph)
    from <- v [match (from, v$id), c ("x", "y")]
    to <- v [match (to, v$id), c ("x", "y")]
    expect_silent (d1 <- dodgr_distances (graph, from = from, to = to))
    expect_identical (as.vector (d0), as.vector (d1))

    # remove from_id/to_id from graph. Now coordinates will be matched on to
    # **first** occurrence in dodgr_vertices, which may not be actual one, so
    # distances won't necessarily be equal
    graph$from_id <- graph$to_id <- NULL
    expect_silent (d2 <- dodgr_distances (graph, from = from, to = to))
    expect_identical (as.vector (d0), as.vector (d2))
})

test_that ("heaps", {
    graph <- weight_streetnet (hampi)
    nf <- 100
    nt <- 50
    from <- sample (graph$from_id, size = nf)
    to <- sample (graph$to_id, size = nt)
    expect_error (
        dodgr_dists (graph, from = from, to = to, heap = "wrong heap"),
        "'arg' should be one of"
    )
    expect_silent (d0 <- dodgr_dists (graph,
        from = from,
        to = to,
        heap = "BHeap"
    ))
    expect_silent (d1 <- dodgr_dists (graph,
        from = from,
        to = to,
        heap = "FHeap"
    ))
    expect_silent (d3 <- dodgr_dists (graph,
        from = from,
        to = to,
        heap = "TriHeap"
    ))
    expect_silent (d4 <- dodgr_dists (graph,
        from = from,
        to = to,
        heap = "TriHeapExt"
    ))
    # This is a compound message that starts "Calculating shortest paths ..."
    # and then "Extended TriHeaps can not be calculated in parallel
    # That can't be tested, so just generic expect_message here
    expect_message (d4a <- dodgr_dists (graph,
        from = from, to = to,
        heap = "TriHeapExt", quiet = FALSE
    ))
    expect_silent (d5 <- dodgr_dists (graph, from = from, to = to, heap = "Heap23"))

    d4 <- dodgr_dists (graph,
        from = from,
        to = to,
        heap = "TriHeapExt",
        quiet = FALSE
    )

    expect_identical (d0, d1)
    expect_identical (d0, d3)
    expect_identical (d0, d4)
    expect_identical (d0, d5)

    # std::set is only applied to non-spatial graphs:
    graph$from_lon <- graph$from_lat <- graph$to_lon <- graph$to_lat <- NULL
    expect_silent (d6 <- dodgr_dists (graph, from = from, to = to, heap = "set"))
    expect_silent (d7 <- dodgr_dists (graph, from = from, to = to, heap = "BHeap"))

    # expect_identical (d0, d6)
    expect_identical (d0, d7)
})

test_that ("graph columns", {
    expect_silent (graph <- weight_streetnet (hampi))
    nf <- 100
    nt <- 50
    set.seed (1)
    v <- dodgr_vertices (graph)
    index_f <- sample (nrow (v), size = nf)
    index_t <- sample (nrow (v), size = nt)
    from <- v$id [index_f]
    to <- v$id [index_t]
    expect_silent (d0 <- dodgr_distances (graph, from = from, to = to))

    from <- v [index_f, c ("x", "y")]
    to <- v [index_t, c ("x", "y")]
    expect_silent (d1 <- dodgr_distances (graph, from = from, to = to))
    colnames (d0) <- colnames (d1) <- rownames (d0) <- rownames (d1) <- NULL
    expect_identical (d0, d1)

    graph$from_lon <- NULL
    expect_error (
        d2 <- dodgr_distances (graph, from = from, to = to),
        "Cannot determine geographical coordinates against which to match pts"
    )

    expect_silent (graph <- weight_streetnet (hampi))
    graph$d_weighted <- graph$d
    expect_silent (d3 <- dodgr_distances (graph, from = from, to = to))
    expect_false (identical (d0, d3))
})

test_that ("negative weights", {
    expect_silent (graph <- weight_streetnet (hampi))
    nf <- 100
    nt <- 50
    set.seed (1)
    from <- sample (graph$from_id, size = nf)
    to <- sample (graph$to_id, size = nt)
    expect_silent (d0 <- dodgr_distances (graph, from = from, to = to))

    nneg <- 100
    graph$d_weighted [sample (nrow (graph), nneg)] <- -runif (nneg)
    expect_silent (d1 <- dodgr_distances (graph, from = from, to = to))
})
