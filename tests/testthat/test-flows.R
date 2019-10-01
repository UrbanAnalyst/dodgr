context("dodgr_flows")

test_all <- (identical (Sys.getenv ("MPADGE_LOCAL"), "true") |
             identical (Sys.getenv ("TRAVIS"), "true"))

testthat::skip_on_cran ()

test_that("flows aggregate", {
    graph <- weight_streetnet (hampi)
    # get routing points from contracted graph:
    graphc <- dodgr_contract_graph (graph)
    from <- sample (graphc$from_id, size = 10)
    to <- sample (graphc$to_id, size = 5)
    to <- to [!to %in% from]
    flows <- matrix (10 * runif (length (from) * length (to)),
                     nrow = length (from))

    expect_message (graph2 <- dodgr_flows_aggregate (graph, from = from, 
                                                     to = to, flows = flows,
                                                     quiet = FALSE),
                    "Aggregating flows ...")
    expect_equal (ncol (graph2) - ncol (graph), 1)
    if (test_all) # fails on CRAN
        expect_true (mean (graph2$flow) > 0)

    flows [1, 2] <- NA
    graph3 <- dodgr_flows_aggregate (graph, from = from, to = to, flows = flows)
    if (test_all)
        #expect_true (max (graph3$flow) <= max (graph2$flow))

    graph4 <- dodgr_flows_aggregate (graph, from = from, to = to, flows = flows,
                                     contract = TRUE)
    if (test_all)
        expect_true (all ((graph4$flow - graph3$flow) < 1e-6))

    expect_warning (graph4 <- dodgr_flows_aggregate (graph3, from = from,
                                                     to = to, flows = flows),
                    "graph already has a 'flow' column; this will be overwritten")

    flowsv <- as.vector (flows)
    graph5 <- dodgr_flows_aggregate (graph, from = from, to = to, flows = flowsv)
    expect_identical (graph5$flow, graph4$flow)
})

test_that("flow points", {
    graph <- weight_streetnet (hampi)
    v <- dodgr_vertices (graph)
    from <- v [sample (nrow (v), size = 10), c ("x", "y")]
    to <- v [sample (nrow (v), size = 10), c ("x", "y")]
    flows <- matrix (10 * runif (length (from) * length (to)),
                     nrow = length (from))

    expect_silent (graph2 <- dodgr_flows_aggregate (graph, from = from,
                                                    to = to, flows = flows))
    expect_true ("flow" %in% names (graph2))
    expect_true (ncol (graph2) == (ncol (graph) + 1))
})

test_that ("flows disperse", {
    graph <- weight_streetnet (hampi)
    from <- sample (graph$from_id, size = 10)
    dens <- runif (length (from))

    expect_message (graph2 <- dodgr_flows_disperse (graph, from = from,
                                                    dens = dens, quiet = FALSE),
                    "Aggregating flows ...")
    expect_equal (ncol (graph2) - ncol (graph), 1)
    if (test_all) # fails on CRAN
        expect_true (mean (graph2$flow) > 0)

    expect_silent (graph3a <- dodgr_flows_disperse (graph, from = from,
                                                    k = 500, dens = dens))
    k <- rep (500, length (from))
    expect_silent (graph3b <- dodgr_flows_disperse (graph, from = from,
                                                    k = k, dens = dens))
    expect_identical (graph3a$flows, graph3b$flows)
    k <- c (k, 1)
    expect_error (graph3b <- dodgr_flows_disperse (graph, from = from,
                                                    k = k, dens = dens),
                  "'k' must be either single value or vector of same")

    expect_warning (graph3 <- dodgr_flows_disperse (graph2, from = from,
                                                    dens = dens),
                    "graph already has a 'flow' column; this will be overwritten")
    #expect_identical (graph3, graph2)
    # flow values are not identical, but
    r2 <- summary (lm (graph3$flow ~ graph2$flow))$r.squared
    if (test_all) # fails on CRAN
        expect_true (r2 > 0.9)

    graph3 <- dodgr_flows_disperse (graph, from = from, dens = dens,
                                    contract = TRUE)
    # Dispersed flows calculated on contracted graph should **NOT** equal those
    # calculated on full graph
    if (test_all) # fails on CRAN
        expect_false (all (graph3$flow == graph2$flow))

    dens [1] <- NA
    graph4 <- dodgr_flows_disperse (graph, from = from, dens = dens)
    # graph4 values are on contracted graph, so flows should generally be less
    # than those on full graph, but may be every so maginally greater
    expect_true (max (graph4$flow - graph2$flow) < 0.1)
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
        expect_silent (dodgr_flowmap (graph_undir))
        a <- dev.off (which = dev.cur ())
        expect_true (file.remove ("junk.png")) # false if no plot

        graph_undir$flow <- NULL
        png (filename = "junk.png")
        expect_silent (dodgr_flowmap (graph_undir))
        a <- dev.off (which = dev.cur ())
        expect_true (file.remove ("junk.png")) # false if no plot
    }

    graph$flow <- NULL
    expect_error (graph_undir <- merge_directed_flows (graph),
                  "graph does not have any flows to merge")
})

