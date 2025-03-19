test_that ("geodist measure", {

    expect_null (getOption ("dodgr_dist_measure"))

    graph <- weight_streetnet (hampi)
    st0 <- system.time (
        m <- get_geodist_measure (graph)
    )
    st1 <- system.time ( # should use option, not calculate
        m <- get_geodist_measure (graph)
    )
    expect_equal (m, "cheap")
    expect_true (st1 [3] < st0 [3])

    op <- getOption ("dodgr_dist_measure")
    expect_true (length (op) > 0L)
    hash <- attr (graph, "hash")
    expect_true (hash %in% names (op))
})
