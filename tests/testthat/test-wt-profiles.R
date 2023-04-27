context ("weighting profiles")

test_that ("wp", {
    f <- file.path (tempdir (), "wp")
    expect_false (file.exists (paste0 (f, ".json")))
    # expect_silent ( # pkg startup msgs on some systems
    write_dodgr_wt_profile (f)
    # )
    expect_true (file.exists (paste0 (f, ".json")))
    w <- read_dodgr_wt_profile (f)
    expect_identical (w, dodgr::weighting_profiles)
})

test_that ("local wt_profile", {
    expect_error (
        write_dodgr_wt_profile (),
        "file name must be given"
    )

    f <- file.path (tempdir (), "wp")
    # expect_silent ( # produces namespace messages on some systems
    write_dodgr_wt_profile (f)
    # )
    n0 <- weight_streetnet (hampi, wt_profile = "foot")
    n1 <- weight_streetnet (hampi,
        wt_profile = "foot",
        wt_profile_file = f
    )
    attr (n0, "px") <- NULL
    attr (n1, "px") <- NULL
    expect_identical (n0, n1)

    w <- dodgr::weighting_profiles
    w$weighting_profiles$max_speed [w$weighting_profiles$name == "foot" &
        w$weighting_profiles$max_speed == 5] <- 8

    f <- paste0 (tools::file_path_sans_ext (f), ".json")
    con <- file (f, open = "wt")
    wpj <- jsonlite::toJSON (w, pretty = TRUE)
    writeLines (wpj, con)
    close (con)

    n2 <- weight_streetnet (hampi,
        wt_profile = "foot",
        wt_profile_file = f
    )
    expect_true (mean (n2$time) < mean (n0$time))
    expect_true (mean (n2$time_weighted) < mean (n0$time_weighted))
    expect_identical (n2$d, n0$d)
    expect_identical (n2$d_weighted, n0$d_weighted)

    expect_error (
        weight_streetnet (hampi, wt_profile = 1:2),
        "wt_profile can only be one element"
    )
})

test_that ("weight_streetnet wt_profile_file", {
    expect_silent (graph <- weight_streetnet (hampi, wt_profile = 1))
    expect_identical (graph$d, graph$d_weighted)
})

test_that ("weight_profile structure", {

    graph0 <- weight_streetnet (hampi, wt_profile = "foot")
    graph1 <- weight_streetnet (hampi, wt_profile = 1)
    expect_equal (nrow (graph0), nrow (graph1))
    expect_identical (graph0$d, graph1$d)
    expect_true (!identical (graph0$d_weighted, graph1$d_weighted))

    wp <- dodgr::weighting_profiles$weighting_profiles
    wpf <- wp [wp$name == "foot", ]
    graph3 <- weight_streetnet (hampi, wt_profile = wpf)
    expect_identical (graph0$d_weighted, graph3$d_weighted)

    wpf$value [wpf$way == "path"] <- 0.9
    graph4 <- weight_streetnet (hampi, wt_profile = wpf)
    expect_true (!identical (graph0$d_weighted, graph4$d_weighted))

    names (wpf) [3] <- "Value"
    expect_error (
        graph4 <- weight_streetnet (hampi, wt_profile = wpf),
        "Weighting profiles must have"
    )

    g0 <- weight_streetnet (hampi)
    if ("px" %in% names (attributes (g0))) {
        while (attr (g0, "px")$is_alive ()) {
            attr (g0, "px")$wait ()
        }
    }
    hampi2 <- hampi
    names (hampi2) [grep ("highway", names (hampi2))] <- "waytype"
    expect_error (
        net <- weight_streetnet (hampi2),
        "Please specify type_col to be used for weighting streetnet"
    )

    g1 <- weight_streetnet (hampi2, type_col = "waytype")
    if ("px" %in% names (attributes (g1))) {
        while (attr (g1, "px")$is_alive ()) {
            attr (g1, "px")$wait ()
        }
    }
    attr (g0, "px") <- NULL
    attr (g1, "px") <- NULL

    expect_identical (g0, g1)
    names (hampi2) [grep ("osm_id", names (hampi2))] <- "key"
    expect_message (
        net <- weight_streetnet (hampi2, type_col = "waytype"),
        "x appears to have no ID column"
    )
    if ("px" %in% names (attributes (net))) {
        while (attr (net, "px")$is_alive ()) {
            attr (net, "px")$wait ()
        }
    }

    names (hampi2) [grep ("key", names (hampi2))] <- "id"
    expect_message (
        net <- weight_streetnet (hampi2, type_col = "waytype"),
        "Using column id as ID column for edges"
    )
    if ("px" %in% names (attributes (net))) {
        while (attr (net, "px")$is_alive ()) {
            attr (net, "px")$wait ()
        }
    }

    names (hampi2) [grep ("width", names (hampi2))] <- "w"
    expect_message (
        g1 <- weight_streetnet (hampi2, type_col = "waytype"),
        "Using column id as ID column for edges"
    )
    if ("px" %in% names (attributes (g1))) {
        while (attr (g1, "px")$is_alive ()) {
            attr (g1, "px")$wait ()
        }
    }
    attr (g0, "px") <- NULL
    attr (g1, "px") <- NULL

    expect_identical (g0, g1)
})


test_that ("railway", {
    expect_error (
        g0 <- weight_streetnet (hampi, wt_profile = "rail"),
        "Please use the weight_railway function for railway routing"
    )
    expect_error (
        g0 <- weight_railway (hampi),
        "Please specify type_col to be used for weighting railway"
    )

    # expect_message (
    #     g0 <- weight_railway (hampi, type_col = "highway"),
    #     "Data has no columns named maxspeed"
    # )
    expect_silent (g0 <- weight_railway (hampi,
        type_col = "highway",
        keep_cols = NULL
    ))

    expect_identical (g0$d, g0$d_weighted)
})
