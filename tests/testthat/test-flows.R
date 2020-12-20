context("dodgr_flows")

test_all <- (identical (Sys.getenv ("MPADGE_LOCAL"), "true"))

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
    #if (test_all)
    #    expect_true (max (graph3$flow) <= max (graph2$flow))

    graph4 <- dodgr_flows_aggregate (graph, from = from, to = to, flows = flows,
                                     contract = TRUE)
    # this test is not longer true with aggregated flows normalised via #121:
    #if (test_all)
    #    expect_true (all ((graph4$flow - graph3$flow) < 1e-3))

    expect_warning (graph4 <- dodgr_flows_aggregate (graph3, from = from,
                                                     to = to, flows = flows),
                "graph already has a 'flow' column; this will be overwritten")

    flowsv <- as.vector (flows)
    graph5 <- dodgr_flows_aggregate (graph,
                                     from = from,
                                     to = to,
                                     flows = flowsv)
    expect_equal (graph5$flow, graph4$flow)
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
    #expect_equal (graph3a$flow, graph3b$flow)
    flow_diff <- abs (graph3a$flow - graph3b$flow)
    expect_true (max (flow_diff) < 1e-4)

    k <- c (500, 1000)
    expect_silent (graph3b <- dodgr_flows_disperse (graph, from = from,
                                                    k = k, dens = dens))
    expect_true (all (c ("flow1", "flow2") %in% names (graph3b)))
    expect_equal (graph3a$flow, graph3b$flow1)

    expect_silent (graph4 <- dodgr_flows_disperse (graph,
                                                   from = from,
                                                   dens = dens,
                                                   contract = TRUE))
    # Dispersed flows calculated on contracted graph should **NOT** equal those
    # calculated on full graph
    if (test_all) # fails on CRAN
        expect_false (all (graph4$flow == graph2$flow))

    dens [1] <- NA
    graph5 <- dodgr_flows_disperse (graph, from = from, dens = dens)
    # graph4 values are on contracted graph, so flows should generally be less
    # than those on full graph, but may be every so maginally greater
    expect_true (max (graph5$flow - graph2$flow) < 0.1)
})

test_that ("flows_si", {
               graph <- weight_streetnet (hampi, wt_profile = "foot") %>%
                   dodgr_contract_graph ()
               v <- dodgr_vertices (graph)
               nf <- 100
               nt <- nrow (v)
               from <- sample (v$id, nf)
               to <- v$id

               k <- 500 + 10 * rnorm (nf)
               dens_from <- 100 * runif (nf)
               dens_to <- 100 * runif (nt)

               # calculation via explicit matrix and flows_aggregate:
               d <- dodgr_distances (graph, from = from, to = to)
               d_from <- array (dens_from, dim = c (nf, nt))
               d_to <- t (array (dens_to, dim = c (nt, nf)))
               kmat <- array (k, dim = c (nf, nt))
               fmat <- d_to * exp (-d / kmat)
               fmat [is.na (fmat)] <- 0
               csmat <- array (rowSums (fmat), dim = c (nf, nt))
               fmat <- d_from * fmat / csmat
               netf <- dodgr_flows_aggregate (graph, from = from, to = to,
                                              flows = fmat)


               # calculation via flows_si:
               netf_si <- dodgr_flows_si (graph,
                                          from = from,
                                          to = to,
                                          k = k,
                                          dens_from = dens_from,
                                          dens_to = dens_to)
               expect_identical (dim (netf), dim (netf_si))
               expect_identical (names (netf), names (netf_si))
               r2 <- cor (netf$flow, netf_si$flow) ^ 2
               if (test_all)
                   expect_true (r2 > 0.5) # sometimes < 0.9
})

test_that ("flowmap", {
    graph <- weight_streetnet (hampi)
    from <- sample (graph$from_id, size = 10)
    to <- sample (graph$to_id, size = 5)
    to <- to [!to %in% from]
    flows <- matrix (10 * runif (length (from) * length (to)),
                     nrow = length (from))
    graph <- dodgr_flows_aggregate (graph, from = from, to = to, flows = flows)
    graph_undir <- merge_directed_graph (graph)

    if (nrow (graph_undir) > 0) {

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
    expect_error (graph_undir <- merge_directed_graph (graph),
                  "col_names \\[flow\\] do not match columns in graph")
})
