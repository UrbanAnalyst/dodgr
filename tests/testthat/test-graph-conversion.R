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

test_that("sf", {
    graph <- weight_streetnet (hampi)
    expect_silent (gsfc <- dodgr_to_sfc (graph))
    expect_is (gsfc, "list")
    expect_equal (length (gsfc), 2)

    gsf1 <- dodgr_to_sf (graph)
    gsf2 <- sf::st_sf (gsfc$dat, geometry = gsfc$geometry, crs = 4326)
    expect_identical (gsf1, gsf2)
})
