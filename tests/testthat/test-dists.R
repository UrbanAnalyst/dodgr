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

test_that ("heaps", {
    graph <- weight_streetnet (hampi)
    nf <- 100
    nt <- 50
    from <- sample (graph$from_id, size = nf)
    to <- sample (graph$from_id, size = nt)
    expect_silent (d0 <- dodgr_dists (graph, from = from, to = to, heap = "BHeap"))
    expect_silent (d1 <- dodgr_dists (graph, from = from, to = to, heap = "FHeap"))
    expect_message (d2 <- dodgr_dists (graph, from = from, to = to, heap = "Radix"),
                    "RadixHeap can only be implemented for integer weights")
    expect_silent (d3 <- dodgr_dists (graph, from = from, to = to, heap = "TriHeap"))
    expect_silent (d4 <- dodgr_dists (graph, from = from, to = to, heap = "TriHeapExt"))
    expect_silent (d5 <- dodgr_dists (graph, from = from, to = to, heap = "Heap23"))

    expect_identical (d0, d1)
    expect_false (identical (d0, d2))
    expect_identical (d0, d3)
    expect_identical (d0, d4)
    expect_identical (d0, d5)

    # std::set is only applied to non-spatial graphs:
    graph$from_lon <- graph$from_lat <- graph$to_lon <- graph$to_lat <- NULL
    expect_silent (d6 <- dodgr_dists (graph, from = from, to = to, heap = "set"))


    expect_identical (d0, d6)
})

test_that("graph columns", {
    expect_silent (graph <- weight_streetnet (hampi))
    nf <- 100
    nt <- 50
    from <- sample (graph$from_id, size = nf)
    to <- sample (graph$from_id, size = nt)
    expect_silent (d0 <- dodgr_distances (graph, from = from, to = to))

    graph <- data.frame (graph) # remove "dodgr_streetnet" class
    expect_silent (d1 <- dodgr_distances (graph, from = from, to = to))
    expect_identical (d0, d1)
})
