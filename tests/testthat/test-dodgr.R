context("dodgr")

test_all <- (identical (Sys.getenv ("MPADGE_LOCAL"), "true") |
             identical (Sys.getenv ("TRAVIS"), "true"))

test_that("dists", {
    graph <- weight_streetnet (hampi)
    from <- sample (graph$from_id, size = 100)
    to <- sample (graph$from_id, size = 50)
    d <- dodgr_dists (graph, from = from, to = to)
    expect_equal (nrow (d), 100)
    expect_equal (ncol (d), 50)
})

test_that("paths", {
    graph <- weight_streetnet (hampi)
    from <- sample (graph$from_id, size = 100)
    to <- sample (graph$from_id, size = 50)
    to <- to [!to %in% from]
    dp <- dodgr_paths (graph, from = from, to = to)
    expect_is (dp, "list")
    expect_equal (length (dp), 100)
    expect_equal (unique (sapply (dp, length)), length (to))
    expect_is (dp [[1]] [[1]], "character")
    lens <- unlist (lapply (dp, function (i) lapply (i, length)))
    dp <- dodgr_paths (graph, from = from, to = to, vertices = FALSE)
    expect_is (dp, "list")
    expect_equal (length (dp), 100)
    expect_equal (unique (sapply (dp, length)), length (to))
    lens2 <- unlist (lapply (dp, function (i) lapply (i, length)))
    # edge lists should all have one less item than vertex lists
    lens2 <- lens2 [which (lens > 0)]
    lens <- lens [which (lens > 0)]
    expect_true (all (lens - lens2 == 1))
})

test_that("flows", {
    graph <- weight_streetnet (hampi)
    from <- sample (graph$from_id, size = 10)
    to <- sample (graph$to_id, size = 5)
    to <- to [!to %in% from]
    flows <- matrix (10 * runif (length (from) * length (to)),
                     nrow = length (from))
    graph2 <- dodgr_flows (graph, from = from, to = to, flows = flows)
    expect_equal (ncol (graph2) - ncol (graph), 1)
    expect_true (mean (graph2$flow) > 0)
})

test_that("sample graph", {
    graph <- weight_streetnet (hampi)
    graph_s <- dodgr_sample (graph, nverts = 100)
    expect_true (nrow (graph_s) < nrow (graph))
    v <- dodgr_vertices (graph_s)
    expect_true (nrow (v) == 100)
})

test_that("components", {
    graph <- weight_streetnet (hampi)
    comp <- graph$component
    graph$component <- NULL
    graph <- dodgr_components (graph)
    expect_identical (comp, graph$component)
})

test_that("contract graph", {
    graph <- weight_streetnet (hampi)
    graph_c <- dodgr_contract_graph (graph)$graph
    expect_true (nrow (graph_c) < nrow (graph))
})

test_that("compare heaps", {
    graph <- weight_streetnet (hampi)
    ch <- compare_heaps (graph, nverts = 100, replications = 1)
    expect_equal (nrow (ch), 11)
    # Test that all dodgr calculations are faster than igraph:
    igr <- which (grepl ("igraph", ch$test))
    #expect_true (ch$elapsed [igr] == max (ch$elapsed))
    # This actually fails on some machines (R oldrel on Windows) so:
    if (test_all)
    {
        expect_true (ch$elapsed [igr] > min (ch$elapsed))
    }
})
