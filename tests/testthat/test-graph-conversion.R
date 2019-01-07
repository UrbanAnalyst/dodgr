context("dodgr graph conversion")

test_that("igraph", {
    graph <- weight_streetnet (hampi)
    graph_i <- dodgr_to_igraph (graph)
    expect_equal (nrow (graph), igraph::ecount (graph_i))
})

test_that("tidyraph", {
    graph <- weight_streetnet (hampi)
    graph_t <- dodgr_to_tidygraph (graph)
    expect_equal (nrow (graph), igraph::ecount (graph_t))
})

