test_all <- (identical (Sys.getenv ("MPADGE_LOCAL"), "true") ||
    identical (Sys.getenv ("GITHUB_JOB"), "test-coverage"))

# testthat::skip_if (!test_all)

test_that ("deduplicate", {

    graph <- weight_streetnet (hampi)
    n0 <- nrow (graph)
    # Duplicate some edges:
    set.seed (1L)
    index <- sample (nrow (graph), size = 10)
    graph <- graph [c (seq (nrow (graph)), index), ]
    n1 <- nrow (graph)

    # expect_message (
    #    duplicated_edge_check (graph, proceed = TRUE),
    #    "Graph has duplicated edges"
    # )

    graph_dedup <- dodgr_deduplicate_graph (graph)
    expect_equal (nrow (graph_dedup), n0)
})
