context("dodgr_dists")

test_all <- (identical (Sys.getenv ("MPADGE_LOCAL"), "true") |
             identical (Sys.getenv ("TRAVIS"), "true"))

test_that("dists", {
    expect_silent (graph <- weight_streetnet (hampi))
    nf <- 100
    nt <- 50
    from <- sample (graph$from_id, size = nf)
    to <- sample (graph$from_id, size = nt)
    expect_silent (d <- dodgr_distances (graph, from = from, to = to))
    expect_equal (nrow (d), nf)
    expect_equal (ncol (d), nt)
    expect_true (all (d [!is.na (d)] >= 0))
    expect_message (d2 <- dodgr_dists (graph, from = from, to = to, quiet = FALSE),
                    "Calculating shortest paths ...")
    expect_identical (d, d2)

    from [1] <- "not_a_vertex_id"
    expect_error (d <- dodgr_distances (graph, from = from, to = to),
                  "from/to are not numeric yet can not be matched onto graph vertices")

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
    from <- as.numeric (from [1, ])
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
    expect_error (d <- dodgr_dists (graph, from = from, to = to),
                  "Unable to determine geographical coordinates in from/to")

    #from <- sample (graph$from_id, size = nf)
    #to <- sample (graph$from_id, size = nt)
    #graph0 <- graph
    #graph <- graph0
    #graph$from_id <- graph$to_id <- NULL
    #find_spatial_cols (graph)

})

test_that("all dists", {
    graph <- weight_streetnet (hampi)
    graph <- graph [graph$component == 2, ]
    expect_silent (d <- dodgr_dists (graph))
    v <- dodgr_vertices (graph)
    expect_equal (nrow (d), ncol (d))
    expect_equal (nrow (d), nrow (v))
})

test_that("dists with no edge ids", {
    graph <- weight_streetnet (hampi)
    nf <- 100
    nt <- 50
    from <- sample (graph$from_id, size = nf)
    to <- sample (graph$from_id, size = nt)
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
    to <- sample (graph$from_id, size = nt)
    expect_error (dodgr_dists (graph, from = from, to = to, heap = "wrong heap"),
                  "'arg' should be one of")
    expect_silent (d0 <- dodgr_dists (graph, from = from, to = to, heap = "BHeap"))
    expect_silent (d1 <- dodgr_dists (graph, from = from, to = to, heap = "FHeap"))
    expect_message (d2 <- dodgr_dists (graph, from = from, to = to, heap = "Radix"),
                    "RadixHeap can only be implemented for integer weights")
    expect_silent (d3 <- dodgr_dists (graph, from = from, to = to, heap = "TriHeap"))
    expect_silent (d4 <- dodgr_dists (graph, from = from, to = to, heap = "TriHeapExt"))
    # This is a compound message that starts "Calculating shortest paths ..."
    # and then "Extended TriHeaps can not be calculated in parallel
    # That can't be tested, so just generic expect_message here
    expect_message (d4a <- dodgr_dists (graph, from = from, to = to,
                                        heap = "TriHeapExt", quiet = FALSE))
    expect_silent (d5 <- dodgr_dists (graph, from = from, to = to, heap = "Heap23"))

    d4 <- dodgr_dists (graph, from = from, to = to, heap = "TriHeapExt", quiet = FALSE)

    expect_identical (d0, d1)
    expect_false (identical (d0, d2))
    expect_identical (d0, d3)
    expect_identical (d0, d4)
    expect_identical (d0, d5)

    # std::set is only applied to non-spatial graphs:
    graph$from_lon <- graph$from_lat <- graph$to_lon <- graph$to_lat <- NULL
    expect_silent (d6 <- dodgr_dists (graph, from = from, to = to, heap = "set"))
    expect_silent (d7 <- dodgr_dists (graph, from = from, to = to, heap = "BHeap"))

    expect_identical (d0, d6)
    expect_identical (d0, d7)
})

test_that("graph columns", {
    expect_silent (graph <- weight_streetnet (hampi))
    nf <- 100
    nt <- 50
    from <- sample (graph$from_id, size = nf)
    to <- sample (graph$from_id, size = nt)
    expect_silent (d0 <- dodgr_distances (graph, from = from, to = to))

    graph$d_weighted <- graph$d
    expect_silent (d1 <- dodgr_distances (graph, from = from, to = to))
    expect_false (identical (d0, d1))
})

test_that ("negative weights", {
    expect_silent (graph <- weight_streetnet (hampi))
    nf <- 100
    nt <- 50
    from <- sample (graph$from_id, size = nf)
    to <- sample (graph$from_id, size = nt)
    expect_silent (d0 <- dodgr_distances (graph, from = from, to = to))

    nneg <- 100
    graph$d_weighted [sample (nrow (graph), nneg)] <- -runif (nneg)
    expect_silent (d1 <- dodgr_distances (graph, from = from, to = to))
})
