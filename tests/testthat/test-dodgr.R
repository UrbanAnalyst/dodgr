context("dodgr")

test_all <- (identical (Sys.getenv ("MPADGE_LOCAL"), "true") |
             identical (Sys.getenv ("TRAVIS"), "true"))

test_that("dists", {
    graph <- weight_streetnet (hampi)
    nf <- 100
    nt <- 50
    from <- sample (graph$from_id, size = nf)
    to <- sample (graph$from_id, size = nt)
    d <- dodgr_dists (graph, from = from, to = to)
    expect_equal (nrow (d), nf)
    expect_equal (ncol (d), nt)
    expect_true (all (d [!is.na (d)] >= 0))

    # dists from coordinates:
    bb <- attr (hampi$geometry, "bbox")
    fromx <- bb [1] + (bb [3] - bb [1]) * runif (nf)
    fromy <- bb [2] + (bb [4] - bb [2]) * runif (nf)
    tox <- bb [1] + (bb [3] - bb [1]) * runif (nt)
    toy <- bb [2] + (bb [4] - bb [2]) * runif (nt)
    from <- data.frame (x = fromx, y = fromy, id = paste0 ("f", 1:nf))
    to <- data.frame (x = tox, y = toy, id = paste0 ("t", 1:nt))
    d <- dodgr_dists (graph, from = from, to = to)
    expect_equal (nrow (d), nf)
    expect_equal (ncol (d), nt)
    expect_true (all (d [!is.na (d)] >= 0))
})

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

test_that("flows aggregate", {
    graph <- weight_streetnet (hampi)
    from <- sample (graph$from_id, size = 10)
    to <- sample (graph$to_id, size = 5)
    to <- to [!to %in% from]
    flows <- matrix (10 * runif (length (from) * length (to)),
                     nrow = length (from))

    graph2 <- dodgr_flows_aggregate (graph, from = from, to = to, flows = flows)
    expect_equal (ncol (graph2) - ncol (graph), 1)
    expect_true (mean (graph2$flow) > 0)

    graph3 <- dodgr_flows_aggregate (graph, from = from, to = to, flows = flows,
                                     contract = TRUE)
    #expect_true (all (graph3$flow == graph2$flow))
    # That's not true because routing points can be betweeen junctions
})

test_that ("flows disperse", {
    graph <- weight_streetnet (hampi)
    from <- sample (graph$from_id, size = 10)
    dens <- runif (length (from))

    graph2 <- dodgr_flows_disperse (graph, from = from, dens = dens)
    expect_equal (ncol (graph2) - ncol (graph), 1)
    expect_true (mean (graph2$flow) > 0)

    graph3 <- dodgr_flows_disperse (graph, from = from, dens = dens,
                                     contract = TRUE)
    # Dispersed flows calculated on contracted graph should **NOT** equal those
    # calculated on full graph
    expect_false (all (graph3$flow == graph2$flow))
})

test_that("sample graph", {
    graph <- weight_streetnet (hampi)
    graph_s <- dodgr_sample (graph, nverts = 100)
    expect_true (nrow (graph_s) < nrow (graph))
    v <- dodgr_vertices (graph_s)
    expect_true (nrow (v) == 100)
})

test_that("components", {
    graph <- weight_streetnet (hampi)
    comp <- graph$component
    graph$component <- NULL
    graph <- dodgr_components (graph)
    expect_identical (comp, graph$component)
})

test_that("contract graph", {
    graph <- weight_streetnet (hampi)
    graph_c <- dodgr_contract_graph (graph)$graph
    expect_true (nrow (graph_c) < nrow (graph))
})

test_that("compare heaps", {
    graph <- weight_streetnet (hampi)
    ch <- compare_heaps (graph, nverts = 100, replications = 1)
    expect_equal (nrow (ch), 11)
    # Test that all dodgr calculations are faster than igraph:
    igr <- which (grepl ("igraph", ch$test))
    #expect_true (ch$elapsed [igr] == max (ch$elapsed))
    # This actually fails on some machines (R oldrel on Windows) so:
    if (test_all)
    {
        expect_true (ch$elapsed [igr] > min (ch$elapsed))
    }
})

test_that("dodgr2sf", {
    hw <- weight_streetnet (hampi)
    expect_warning (dodgr_to_sf (hw), "'dodgr_to_sf' is deprecated")
    y <- dodgr_to_sfc (hw)
    # y should have more linestrings than the original sf object:
    expect_true (length (y) > length (hw$geometry))
})

test_that ("flowmap", {
    graph <- weight_streetnet (hampi)
    from <- sample (graph$from_id, size = 10)
    to <- sample (graph$to_id, size = 5)
    to <- to [!to %in% from]
    flows <- matrix (10 * runif (length (from) * length (to)),
                     nrow = length (from))
    graph <- dodgr_flows_aggregate (graph, from = from, to = to, flows = flows)
    graph_undir <- merge_directed_flows (graph)
    if (nrow (graph_undir) > 0)
    {
        # just test that is produces a plot
        png (filename = "junk.png")
        dodgr_flowmap (graph_undir)
        a <- dev.off (which = dev.cur ())
        expect_true (file.remove ("junk.png")) # false if no plot
    }
})

test_that ("weight_profiles", {
    graph0 <- weight_streetnet (hampi, wt_profile = "foot")
    graph1 <- weight_streetnet (hampi, wt_profile = 1)
    # hampi has some undefined highway types which should be removed from graph,
    # giving less rows in streetnet
    expect_true (nrow (graph0) < nrow (graph1))
    #expect_identical (graph0$d, graph1$d)
    #expect_true (!identical (graph0$d_weighted, graph1$d_weighted))
    wtp <- dodgr::weighting_profiles [dodgr::weighting_profiles == "foot", ]
    graph3 <- weight_streetnet (hampi, wt_profile = wtp)
    expect_identical (graph0$d_weighted, graph3$d_weighted)

    wtp$value [wtp$way == "path"] <- 0.9
    graph4 <- weight_streetnet (hampi, wt_profile = wtp)
    expect_true (!identical (graph0$d_weighted, graph4$d_weighted))

    names (wtp) [3] <- "Value"
    expect_error (graph4 <- weight_streetnet (hampi, wt_profile = wtp),
                  "Weighting profiles must have")
})
