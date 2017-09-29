context("dodgr")

test_that("dists", {
    graph <- weight_streetnet (hampi)
    from <- sample (graph$from_id, size = 100)
    to <- sample (graph$from_id, size = 50)
    d <- dodgr_dists (graph, from = from, to = to)
    expect_equal (nrow (d), 100)
    expect_equal (ncol (d), 50)
})

test_that("convert graph", {
    graph <- weight_streetnet (hampi)
    graph2 <- dodgr_convert_graph (graph)$graph
    expect_equal (ncol (graph2), 6)
    expect_true (ncol (graph2) < ncol (graph))
})

test_that("sample graph", {
    graph <- weight_streetnet (hampi)
    graph_s <- dodgr_sample (graph, nverts = 100)
    expect_true (nrow (graph_s) < nrow (graph))
    v <- dodgr_vertices (graph_s)
    expect_equal (nrow (v), 100)
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
    graph_c <- dodgr_contract_graph (graph)
    expect_true (nrow (graph_c) < nrow (graph))
})

test_that("compare heaps", {
    graph <- weight_streetnet (hampi)
    ch <- compare_heaps (graph, nverts = 100, replications = 1)
    expect_equal (nrow (ch), 11)
    # Test that all dodgr calculations are faster than igraph:
    igr <- which (grepl ("igraph", ch$test))
    expect_true (ch$elapsed [igr] == max (ch$elapsed))
})
