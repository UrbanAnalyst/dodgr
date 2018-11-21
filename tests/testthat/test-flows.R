context("dodgr_flows")

test_all <- (identical (Sys.getenv ("MPADGE_LOCAL"), "true") |
             identical (Sys.getenv ("TRAVIS"), "true"))

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

