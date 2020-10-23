context("centrality")

test_all <- (identical (Sys.getenv ("MPADGE_LOCAL"), "true"))

testthat::skip_on_cran ()

test_that("centrality", {
    graph_full <- weight_streetnet (hampi)
    graph <- dodgr_contract_graph (graph_full)

    graphc <- dodgr_centrality (graph)
    expect_equal (nrow (graph), nrow (graphc))
    expect_equal (ncol (graph) + 1, ncol (graphc))
    expect_true ("centrality" %in% names (graphc))
    expect_false ("centrality" %in% names (graph))

    v <- dodgr_vertices (graph)
    vc <- dodgr_centrality (graph, edges = FALSE)
    expect_equal (nrow (v), nrow (vc))
    expect_equal (ncol (v) + 1, ncol (vc))
    expect_true ("centrality" %in% names (vc))
    expect_false ("centrality" %in% names (v))
})

test_that("estimate time", {
    graph <- weight_streetnet (hampi)
    expect_message (x <- estimate_centrality_time (graph),
                    "Estimated time to calculate centrality for full graph is")
    if (test_all)
        expect_identical (x, "00:00:00")
})
