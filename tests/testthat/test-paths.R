context("dodgr_paths")

test_all <- (identical (Sys.getenv ("MPADGE_LOCAL"), "true") |
             identical (Sys.getenv ("TRAVIS"), "true"))

test_that("paths", {
    graph <- weight_streetnet (hampi)
    from <- sample (graph$from_id, size = 100)
    to <- sample (graph$from_id, size = 50)
    to <- to [!to %in% from]
    dp <- dodgr_paths (graph, from = from, to = to)
    expect_is (dp, "list")
    expect_equal (length (dp), 100)
    expect_equal (unique (sapply (dp, length)), length (to))
    expect_is (dp [[1]] [[1]], "character")
    lens <- unlist (lapply (dp, function (i) lapply (i, length)))
    dp <- dodgr_paths (graph, from = from, to = to, vertices = FALSE)
    expect_is (dp, "list")
    expect_equal (length (dp), 100)
    expect_equal (unique (sapply (dp, length)), length (to))
    lens2 <- unlist (lapply (dp, function (i) lapply (i, length)))
    # edge lists should all have one less item than vertex lists
    lens2 <- lens2 [which (lens > 0)]
    lens <- lens [which (lens > 0)]
    expect_true (all (abs (lens - lens2) == 1))
})

test_that("pairwise paths", {
    graph <- weight_streetnet (hampi)
    from <- sample (graph$from_id, size = 10)
    to <- sample (graph$from_id, size = 10)
    indx <- which (!to %in% from)
    to <- to [indx]
    from <- from [indx]
    n <- length (indx)
    dp <- dodgr_paths (graph, from = from, to = to, pairwise = TRUE)
    expect_is (dp, "list")
    expect_equal (length (dp), n)
    expect_true (all (lapply (dp, length) == 1))
})
