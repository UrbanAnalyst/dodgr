context ("dodgr_dists_categorical")

test_all <- (identical (Sys.getenv ("MPADGE_LOCAL"), "true") |
    identical (Sys.getenv ("GITHUB_WORKFLOW"), "test-coverage"))

if (!test_all) {
    RcppParallel::setThreadOptions (numThreads = 2)
}

test_that ("categorical dists", {

    expect_silent (graph <- weight_streetnet (hampi))

    nf <- 100
    nt <- 50
    set.seed (1)
    from <- sample (graph$from_id, size = nf)
    to <- sample (graph$to_id, size = nt)

    expect_error (
        d <- dodgr_dists_categorical (
            graph,
            from,
            to
        ),
        "graph must have a column named 'edge_type'"
    )

    graph$edge_type <- 1L
    graph$edge_type [1] <- 0L
    expect_error (
        d <- dodgr_dists_categorical (
            graph,
            from,
            to
        ),
        "graphs with integer edge_type columns may not contain 0s"
    )

    graph$edge_type <- graph$highway
    expect_silent (d <- dodgr_dists_categorical (graph, from = from, to = to))
    expect_type (d, "list")
    ntypes <- length (unique (graph$highway))
    expect_length (d, ntypes + 1L)
    expect_true (all (unique (graph$highway) %in% names (d)))

    # All dimensions should be equal:
    dims <- vapply (d, dim, integer (2))
    expect_identical (sum (diff (dims [1, ])), 0L)
    expect_true (all (dims [1, ] == nf))
    expect_true (all (dims [2, ] == nt))

    expect_message (
        d2 <- dodgr_dists_categorical (
            graph,
            from = from,
            to = to,
            quiet = FALSE
        ),
        "Calculating shortest paths ..."
    )
    expect_identical (d, d2)
})

test_that ("categorical dists results", {

    net <- weight_streetnet (hampi, wt_profile = "foot")
    v <- dodgr_vertices (net)
    set.seed (1L)
    from <- sample (v$id, 20)
    to <- sample (v$id, 20)

    net$edge_type <- net$highway
    d0 <- dodgr_dists_categorical (net, from, to, proportions_only = FALSE, pairwise = FALSE)
    # Sums of all "type" matrices should equal main "distances" matrix:
    types <- which (names (d0) != "distances")
    d_types <- lapply (types, function (i) colSums (d0 [[i]]))
    d_types <- colSums (do.call (rbind, d_types))
    d_total <- colSums (d0$distances)
    expect_true (all (abs (d_total - d_types) < 1e-6))

    d1 <- dodgr_dists_categorical (net, from, to, proportions_only = FALSE, pairwise = TRUE)
    index_rows <- which (d1 [, 1] > 0)
    index_cols <- which (colnames (d1) != "total")
    d_total <- d1 [index_rows, 1L]
    d_types <- rowSums (d1 [index_rows, index_cols])
    expect_true (all (abs (d_total - d_types) < 1e-6))

    p0 <- dodgr_dists_categorical (net, from, to, proportions_only = TRUE)
    dp <- vapply (d0, function (i) sum (colSums (i)), numeric (1L))
    p1 <- dp [-1] / dp [1]
    expect_true (all (abs (p0 - p1) < 1e-6))
})

test_that ("categorical dists summary", {

    expect_silent (graph <- weight_streetnet (hampi))
    graph <- graph [graph$component == 1, ]

    v <- dodgr_vertices (graph)
    from <- v$id
    to <- v$id

    graph$edge_type <- graph$highway

    expect_silent (d <- dodgr_dists_categorical (graph, from, to))

    expect_message (ds <- summary (d))
    expect_is (ds, "numeric")
    expect_equal (length (ds), length (unique (graph$edge_type)))
    expect_true (all (ds < 1.0))
})

test_that ("proportions only", {

    expect_silent (graph <- weight_streetnet (hampi))
    graph <- graph [graph$component == 1, ]

    v <- dodgr_vertices (graph)
    from <- v$id
    to <- v$id

    graph$edge_type <- graph$highway

    expect_silent (d <- dodgr_dists_categorical (graph, from, to))

    d0 <- d$distances
    d <- d [-1]
    dtypes <- vapply (
        d, function (i) sum (colSums (i)),
        numeric (1)
    )
    d1 <- dtypes / sum (colSums (d0))

    expect_silent (d2 <- dodgr_dists_categorical (graph, from, to,
        proportions_only = TRUE
    ))
    # These will only be approximately equal:
    expect_true (mean (abs (d1 - d2)) > 0)
})

test_that ("categorical threshold", {

    expect_silent (graph <- weight_streetnet (hampi))
    graph <- graph [graph$component == 1, ]

    v <- dodgr_vertices (graph)
    from <- v$id

    graph$edge_type <- graph$highway

    expect_silent (d <- dodgr_dists_categorical (graph, from, dlimit = 2000))

    expect_is (d, "data.frame")
    expect_equal (nrow (d), length (from))
    expect_equal (ncol (d), length (unique (graph$edge_type)) + 1L)

    # d [, 1] is distance; all other cols are edge-type-specific
    for (i in 2:ncol (d)) {
        expect_true (all (d [, i] <= d [, 1]))
    }
})

test_that ("categorical pairwise", {

    expect_silent (graph <- weight_streetnet (hampi))
    graph <- graph [graph$component == 1, ]

    v <- dodgr_vertices (graph)
    set.seed (1L)
    from <- sample (v$id, size = 100)
    to <- sample (v$id, size = 100)

    graph$edge_type <- graph$highway

    expect_silent (
        d <- dodgr_dists_categorical (graph, from, to, pairwise = TRUE)
    )

    expect_is (d, "matrix")
    expect_equal (nrow (d), length (from))
    expect_equal (ncol (d), length (unique (graph$edge_type)) + 1L)

    # d [, 1] is distance; all other cols are edge-type-specific
    for (i in 2:ncol (d)) {
        expect_true (all (d [, i] <= d [, 1]))
    }

    # rownames are "<from>-<to>":
    from_id <- gsub ("\\-.*$", "", rownames (d))
    to_id <- gsub ("^.*\\-", "", rownames (d))
    expect_identical (from_id, from)
    expect_identical (to_id, to)
})
