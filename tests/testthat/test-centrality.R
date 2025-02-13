test_all <- (identical (Sys.getenv ("MPADGE_LOCAL"), "true") ||
    identical (Sys.getenv ("GITHUB_JOB"), "test-coverage"))

# testthat::skip_on_cran ()
testthat::skip_if (!test_all)

test_that ("centrality", {
    graph_full <- weight_streetnet (hampi)
    graph <- dodgr_contract_graph (graph_full)

    graphc <- dodgr_centrality (graph, check_graph = FALSE)
    expect_equal (nrow (graph), nrow (graphc))
    expect_equal (ncol (graph) + 1, ncol (graphc))
    expect_true ("centrality" %in% names (graphc))
    expect_false ("centrality" %in% names (graph))

    v <- dodgr_vertices (graph)
    vc <- dodgr_centrality (graph, edges = FALSE, check_graph = FALSE)
    expect_equal (nrow (v), nrow (vc))
    expect_equal (ncol (v) + 1, ncol (vc))
    expect_true ("centrality" %in% names (vc))
    expect_false ("centrality" %in% names (v))
})

test_that ("weighted centrality", {
    graph_full <- weight_streetnet (hampi)
    graph <- dodgr_contract_graph (graph_full)
    graphc0 <- dodgr_centrality (graph, check_graph = FALSE)
    expect_equal (nrow (graphc0), nrow (graph))

    v <- dodgr_vertices (graph)
    set.seed (1)
    vert_wts <- runif (nrow (v))
    graphc1 <- dodgr_centrality (graph, vert_wts = vert_wts, check_graph = FALSE)
    expect_equal (nrow (graphc1), nrow (graph))
    expect_true (!all (graphc1$centrality == graphc0$centrality))

    v0 <- dodgr_centrality (graph, edges = FALSE, check_graph = FALSE)
    expect_equal (nrow (dodgr_vertices (graph)), nrow (v0))

    vert_wts <- runif (nrow (v))
    v1 <- dodgr_centrality (graph, vert_wts = vert_wts, edges = FALSE, check_graph = FALSE)
    expect_equal (nrow (v1), nrow (v0))
    expect_true (!all (v1$centrality == v0$centrality))

    vert_wts <- NULL
    expect_silent (
        vx <- dodgr_centrality (graph, vert_wts = vert_wts, edges = FALSE, check_graph = FALSE)
    )

    vert_wts <- runif (nrow (v)) [-1]
    expect_error (
        vx <- dodgr_centrality (graph, vert_wts = vert_wts, edges = FALSE, check_graph = FALSE),
        "vert_wts must be a vector of same length"
    )

    vert_wts <- "A"
    expect_error (
        vx <- dodgr_centrality (graph, vert_wts = vert_wts, edges = FALSE, check_graph = FALSE),
        "vert_wts must be a vector of same length"
    )
})

test_that ("estimate time", {
    graph <- weight_streetnet (hampi)
    expect_message (
        x <- estimate_centrality_time (graph),
        "Estimated time to calculate centrality for full graph is"
    )
    # Convert `x` as H:M:S to integer
    x <- as.integer (strsplit (x, ":") [[1]])
    x <- sum (x * c (3600, 60, 1))
    if (test_all) {
        expect_true (x <= 1L)
    }
})
