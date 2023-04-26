test_all <- (identical (Sys.getenv ("MPADGE_LOCAL"), "true") |
    identical (Sys.getenv ("GITHUB_WORKFLOW"), "test-coverage"))

if (!test_all) {
    RcppParallel::setThreadOptions (numThreads = 2)
}

test_that ("categorical nearest dists", {

    expect_silent (graph <- weight_streetnet (hampi))

    nf <- 50
    nt <- 100
    set.seed (1)
    from <- sample (graph$from_id, size = nf)
    to <- sample (graph$to_id, size = nt)

    expect_silent (
        d <- dodgr_dists_nearest (graph, from, to)
    )

    expect_s3_class (d, "data.frame")
    expect_equal (ncol (d), 3L)
    expect_equal (nrow (d), length (from))
    expect_identical (names (d), c ("from", "to", "d"))
    expect_type (d$from, "character")
    expect_type (d$to, "character")
    expect_type (d$d, "double")

    expect_identical (d$from, from)
    expect_false (all (to %in% d$to))
    expect_true (all (d$to [which (!is.na (d$to))] %in% to))
    expect_true (min (d$d, na.rm = TRUE) >= 0.0)
})
