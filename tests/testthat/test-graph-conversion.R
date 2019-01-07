context("dodgr graph conversion")

test_that("igraph", {
    graph <- weight_streetnet (hampi)
    graph_i <- dodgr_to_igraph (graph)
    expect_equal (nrow (graph), igraph::ecount (graph_i))
})

