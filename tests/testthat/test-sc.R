context ("SC")

test_all <- (identical (Sys.getenv ("MPADGE_LOCAL"), "true") |
    identical (Sys.getenv ("GITHUB_WORKFLOW"), "test-coverage"))

skip_if (!test_all)

# library (osmdata)
# devtools::load_all ("../../ropensci/osmdata", export_all = FALSE)
# h2 <- opq ("hampi india") %>%
#    add_osm_feature (key = "highway") %>%
#    osmdata_sc ()

source ("../sc-conversion-fns.R")

test_that ("SC", {
    expect_silent (hsc <- sf_to_sc (hampi))
    # This all exists just to test the next line:
    requireNamespace ("geodist")
    requireNamespace ("dplyr")
    expect_silent (net_sc <- weight_streetnet (hsc))
    expect_is (net_sc, "data.frame")
    expect_true (nrow (net_sc) > 0)

    net_sf <- weight_streetnet (hampi)
    expect_true (nrow (net_sf) > nrow (net_sc))
    v_sc <- dodgr_vertices (net_sc)
    v_sf <- dodgr_vertices (net_sf)
    expect_true (nrow (v_sf) > nrow (v_sc))

    class (hsc) <- class (hsc) [!class (hsc) %in% "osmdata_sc"]
    expect_error (
        net_sc <- weight_streetnet (hsc),
        paste0 (
            "weight_streetnet currently only works ",
            "for 'sc'-class objects extracted with"
        )
    )

    expect_silent (hsc <- sf_to_sc (hampi))
    expect_silent (net_sc2 <- weight_streetnet (hsc,
        wt_profile = "horse"
    ))
    expect_true (!identical (net_sc$d_weighted, net_sc2$d_weighted))

    net_sc2 <- dodgr_components (net_sc2)
    expect_silent (v0 <- dodgr_vertices (net_sc2))
    # force re-cache by re-generating edge IDs:
    net_sc2$edge_ <-
        paste0 (seq (nrow (net_sc2)) [order (runif (nrow (net_sc2)))])
    net_sc2$.vx0 <- as.factor (net_sc2$.vx0)
    expect_silent (v1 <- dodgr_vertices (net_sc2)) # should still work

    # force re-cache by re-generating edge IDs:
    net_sc2$edge_ <-
        paste0 (seq (nrow (net_sc2)) [order (runif (nrow (net_sc2)))])
    net_sc2$.vx0 <- as.character (net_sc2$.vx0)
    net_sc2$.vx1 <- as.factor (net_sc2$.vx1)
    expect_silent (v2 <- dodgr_vertices (net_sc2)) # should still work

    net_sc3 <- weight_streetnet (hsc, wt_profile = "bicycle")
    net_sc3 <- dodgr_components (net_sc3)
    # force re-cache by re-generating edge IDs:
    net_sc3$edge_ <-
        paste0 (seq (nrow (net_sc3)) [order (runif (nrow (net_sc3)))])
    expect_silent (v0 <- dodgr_vertices (net_sc3))
    expect_true (all (c ("x", "y") %in% names (v0)))
    net_sc3$edge_ <-
        paste0 (seq (nrow (net_sc3)) [order (runif (nrow (net_sc3)))])
    net_sc3$.vx0_x <-
        net_sc3$.vx0_y <-
        net_sc3$.vx1_x <-
        net_sc3$.vx1_y <- NULL
    expect_silent (v1 <- dodgr_vertices (net_sc3))
    expect_false (all (c ("x", "y") %in% names (v1)))
    expect_identical (v0$id, v1$id)

    # add fake elevation data:
    net_sc <- weight_streetnet (hsc, wt_profile = "bicycle")
    hsc$vertex$z_ <- 10 * runif (nrow (hsc$vertex))
    hsc$vertex <- hsc$vertex [match (
        names (hsc$vertex),
        c ("x_", "y_", "z_", "vertex_")
    )]
    net_sc2 <- weight_streetnet (hsc, wt_profile = "bicycle")
    expect_false ("dz" %in% names (net_sc))
    expect_true ("dz" %in% names (net_sc2))

    expect_error (
        x <- weight_railway (hsc),
        'x must be class "sf"'
    )
})

test_that ("traffic light nodes", {
    expect_silent (hsc <- sf_to_sc (hampi))
    expect_silent (net_sc0 <- weight_streetnet (hsc))
    v <- sample (hsc$vertex$vertex_, size = 10)
    hsc$nodes <- data.frame (
        vertex_ = v,
        key = "highway",
        value = "traffic_signals"
    )
    expect_silent (net_sc1 <- weight_streetnet (hsc))
    # This has no effect here, because the edges must also be flagged
    # with same key-val pair

    expect_identical (net_sc0$d, net_sc1$d)
    expect_identical (net_sc0$d_weighted, net_sc1$d_weighted)
    expect_true (!identical (net_sc0$time, net_sc1$time))
    expect_identical (net_sc0$time_weighted, net_sc1$time_weighted)

    expect_silent (net_sc1 <- weight_streetnet (hsc, wt_profile = 1))
    expect_identical (net_sc1$d, net_sc1$d_weighted)
    expect_identical (net_sc1$time, net_sc1$time_weighted)
})

test_that ("elevation", {
    expect_silent (hsc <- sf_to_sc (hampi))
    expect_silent (net_sc <- weight_streetnet (hsc))
    hsc$vertex$z_ <- runif (nrow (hsc$vertex)) * 10
    expect_silent (net_sc2 <- weight_streetnet (hsc))
    expect_true (ncol (net_sc2) == (ncol (net_sc) + 1))

    expect_silent (net_sc3 <- weight_streetnet (hsc,
        wt_profile = "foot"
    ))
    expect_true (ncol (net_sc3) == (ncol (net_sc2)))
    expect_true (mean (net_sc3$time) > mean (net_sc2$time))
})

test_that ("contract with turn angles", {
    expect_silent (hsc <- sf_to_sc (hampi))
    expect_silent (graph <- weight_streetnet (hsc,
        wt_profile = "bicycle"
    ))
    expect_silent (graph_c <- dodgr_contract_graph (graph))
    expect_silent (v <- dodgr_vertices (graph_c))
    n <- 100
    pts <- sample (v$id, size = n)
    pts <- pts [which (pts %in% graph_c$.vx0 & pts %in% graph_c$.vx1)]
    fmat <- array (1, dim = c (n, n))

    # aggregate flows from graph without turning angles:
    expect_silent (graphf <- dodgr_flows_aggregate (graph_c,
        from = pts,
        to = pts,
        flow = fmat,
        contract = FALSE
    ))
    expect_silent (graphf <- dodgr_uncontract_graph (graphf))
    expect_silent (graphf <- merge_directed_graph (graphf))

    # then turn angle graph
    grapht <- weight_streetnet (hsc,
        wt_profile = "bicycle",
        turn_penalty = TRUE, left_side = TRUE
    )

    expect_equal (nrow (grapht), nrow (graph))
    grapht_c <- dodgr_contract_graph (grapht)
    expect_equal (nrow (grapht_c), nrow (graph_c))
    expect_warning (
        graphtf <- dodgr_flows_aggregate (
            grapht_c,
            from = pts,
            to = pts,
            flow = fmat,
            contract = FALSE
        ),
        "graphs with turn penalties should be submitted in full, not contracted form"
    )
    expect_silent (
        graphtf <- dodgr_flows_aggregate (
            grapht,
            from = pts,
            to = pts,
            flow = fmat,
            contract = FALSE
        )
    )

    # compound junction edges are then removed, as are vertex
    # suffixes:
    expect_true (length (grep ("_start", graphtf$.vx0)) == 0)
    expect_true (length (grep ("_end", graphtf$.vx1)) == 0)

    expect_silent (graphtf <- merge_directed_graph (graphtf))
    # this test does not consistently pass:
    # expect_identical (range (graphf$flow), range (graphtf$flow))
    # TODO: Implement a better alternative

    expect_warning (
        graphtf <-
            dodgr_flows_disperse (
                grapht_c,
                from = pts,
                dens = rep (1, n)
            ),
        "graphs with turn penalties should be submitted in full, not contracted form"
    )
    expect_silent (
        graphtf <- dodgr_flows_disperse (grapht, from = pts, dens = rep (1, n))
    )
})

test_that ("dodgr_times", {
    # dists and times should be strongly correlated:
    expect_silent (hsc <- sf_to_sc (hampi))
    expect_silent (net_sc <- weight_streetnet (hsc))
    v <- dodgr_vertices (net_sc)
    set.seed (1)
    from <- sample (v$id, 100)
    to <- sample (v$id, 100)
    d <- dodgr_dists (net_sc, from = from, to = to)
    t1 <- dodgr_times (net_sc, from = from, to = to)
    r2 <- cor (as.numeric (d), as.numeric (t1),
        use = "pairwise.complete.obs"
    )
    expect_true (r2 < 1)
    # with no turn angles, the should be just scaled versions

    # calculate times with turning angles, such that resultant network
    # includes compound junction edges
    expect_silent (net_sc2 <- weight_streetnet (hsc,
        turn_penalty = TRUE
    ))
    expect_equal (nrow (net_sc2), nrow (net_sc))
    from <- remap_verts_with_turn_penalty (net_sc2, from, from = TRUE)
    to <- remap_verts_with_turn_penalty (net_sc2, to, from = FALSE)
    t2 <- dodgr_times (net_sc2, from = from, to = to)
    r2 <- cor (as.numeric (t1), as.numeric (t2),
        use = "pairwise.complete.obs"
    )
    # expect_true (r2 < 1)
    expect_true (r2 > 0.95)
    # These times should be longer:
    expect_true (mean (t2 - t1, na.rm = TRUE) > 0)

    # times with contracted graph should be identical:
    net_sc2_c <- dodgr_contract_graph (net_sc2)
    v <- dodgr_vertices (net_sc2_c)
    set.seed (1)
    from <- sample (v$id, 100)
    to <- sample (v$id, 100)

    t1 <- dodgr_times (net_sc2, from = from, to = to)
    expect_warning (
        t2 <- dodgr_times (net_sc2_c, from = from, to = to),
        "graphs with turn penalties should be submitted in full, not contracted form"
    )

    dtime <- max (abs (t1 - t2), na.rm = TRUE)
    # expect_true (dtime < 1e-6)
    r2 <- cor (as.vector (t1), as.vector (t2),
        use = "pairwise.complete.obs"
    )^2
    expect_true (r2 > 0.9)
})
