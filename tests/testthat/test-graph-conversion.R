context ("dodgr graph conversion")

test_that ("igraph", {
    graph <- weight_streetnet (hampi)
    graph_i <- dodgr_to_igraph (graph)
    expect_equal (nrow (graph), igraph::ecount (graph_i))

    graph2 <- igraph_to_dodgr (graph_i)
    expect_true (!identical (graph, graph2))
    expect_equal (nrow (graph), nrow (graph2))

    expect_error (
        graph_i <- dodgr_to_igraph (graph,
            weight_column = "does_not_exist"
        ),
        "graph contains no column named 'does_not_exist'"
    )

    graph$from_id <- graph$to_id <- NULL
    expect_error (
        graph_i <- dodgr_to_igraph (graph),
        "Some vertex names in edge list are not listed in vertex data frame"
    )
})

test_that ("tidyraph", {
    graph <- weight_streetnet (hampi)
    graph_t <- dodgr_to_tidygraph (graph)
    expect_equal (nrow (graph), igraph::ecount (graph_t))
})

test_that ("sf", {
    graph <- weight_streetnet (hampi)
    gsfc <- dodgr_to_sfc (graph)
    expect_is (gsfc, "list")
    expect_equal (length (gsfc), 2)

    gsf1 <- dodgr_to_sf (graph)
    gsf2 <- sf::st_sf (gsfc$dat, geometry = gsfc$geometry, crs = 4326)
    expect_identical (gsf1, gsf2)

    gc <- dodgr_contract_graph (graph)
    msg <- "Calling on a contracted graph will result in loss of information"
    expect_warning (gsf <- dodgr_to_sf (gc), msg)
})
