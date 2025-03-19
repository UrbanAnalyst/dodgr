test_all <- (identical (Sys.getenv ("MPADGE_LOCAL"), "true") ||
    identical (Sys.getenv ("GITHUB_JOB"), "test-coverage"))

test_that ("geodist measure", {

    # Reset any measures stored in options:
    options ("dodgr_dist_measure" = NULL)

    graph <- weight_streetnet (hampi)
    st0 <- system.time (
        m <- get_geodist_measure (graph)
    )
    st1 <- system.time ( # should use option, not calculate
        m <- get_geodist_measure (graph)
    )
    expect_equal (m, "cheap")
    if (test_all) {
        expect_true (st1 [3] < st0 [3])
    }

    op <- getOption ("dodgr_dist_measure")
    expect_true (length (op) > 0L)
    hash <- attr (graph, "hash")
    expect_true (hash %in% names (op))
})
