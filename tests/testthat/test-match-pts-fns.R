test_all <- (identical (Sys.getenv ("MPADGE_LOCAL"), "true") |
    identical (Sys.getenv ("GITHUB_WORKFLOW"), "test-coverage"))

skip_if (!test_all)

dodgr_cache_off ()
clear_dodgr_cache ()

test_that ("points to verts", {

    bb <- attr (hampi$geometry, "bbox")
    n <- 100
    x <- bb [1] + (bb [3] - bb [1]) * runif (n)
    y <- bb [2] + (bb [4] - bb [2]) * runif (n)
    pts <- data.frame (x = x, y = y)
    net <- weight_streetnet (hampi)
    expect_message (
        index1 <- match_pts_to_verts (net, pts),
        paste0 (
            "First argument to match_pts_to_verts should ",
            "be result of dodgr_vertices"
        )
    )

    v <- dodgr_vertices (net)
    expect_silent (index2 <- match_pts_to_verts (v, pts))
    expect_identical (index1, index2)

    colnames (pts) <- NULL
    expect_message (
        index3 <- match_pts_to_verts (v, pts),
        "xy has no named columns; assuming order is x then y"
    )
    expect_identical (index1, index3)

    pts <- data.frame (x = x, y = y, x2 = x)
    expect_error (
        index4 <- match_pts_to_verts (v, list (pts)),
        "xy must be a matrix or data.frame"
    )
    expect_error (
        index4 <- match_pts_to_verts (v, pts),
        "xy must have only two columns"
    )

    pts <- data.frame (x = x, y = y)
    expect_silent (index4 <- match_pts_to_verts (v, pts, connected = TRUE))
    expect_true (!identical (index1, index4))

    class (pts) <- c (class (pts), "tbl")
    expect_silent (index5 <- match_pts_to_verts (v, pts, connected = TRUE))
    expect_identical (index4, index5)

    pts <- sf::st_as_sf (pts, coords = c (1, 2), crs = 4326)
    expect_silent (index6 <- match_pts_to_verts (v, pts, connected = TRUE))
    expect_identical (index4, index6)
    expect_silent (index7 <- match_pts_to_verts (v, pts, connected = TRUE))
    expect_identical (index4, index7)

    pts <- hampi [1, ]
    expect_error (index7 <- match_pts_to_verts (v, pts))
    # error is "xy$geometry must be a collection of sfc_POINT objects", but
    # expect_error does not match on the "$" symbo, but expect_error does not
    # match on the "$" symbol

})

test_that ("points to graph", {

    bb <- attr (hampi$geometry, "bbox")
    n <- 100
    x <- bb [1] + (bb [3] - bb [1]) * runif (n)
    y <- bb [2] + (bb [4] - bb [2]) * runif (n)
    pts <- data.frame (x = x, y = y)
    net <- weight_streetnet (hampi)
    verts <- dodgr_vertices (net)

    expect_error (
        match_pts_to_graph (verts, pts),
        "Points may only be matched to spatial graphs."
    )

    expect_silent (index1 <- match_pts_to_graph (net, pts))

    colnames (pts) <- NULL
    expect_message (
        index2 <- match_pts_to_graph (net, pts),
        "xy has no named columns; assuming order is x then y"
    )
    expect_identical (index1, index2)

    pts <- data.frame (x = x, y = y, x2 = x)
    expect_error (
        match_pts_to_graph (net, list (pts)),
        "xy must be a matrix or data.frame"
    )
    expect_error (
        match_pts_to_graph (net, pts),
        "xy must have only two columns"
    )

    pts <- data.frame (x = x, y = y)
    expect_silent (index4 <- match_pts_to_graph (net, pts, connected = TRUE))
    expect_true (!identical (index1, index4))

    class (pts) <- c (class (pts), "tbl")
    expect_silent (index5 <- match_pts_to_graph (net, pts, connected = TRUE))
    expect_identical (index4, index5)

    pts <- sf::st_as_sf (pts, coords = c (1, 2), crs = 4326)
    expect_silent (index6 <- match_pts_to_graph (net, pts, connected = TRUE))
    expect_identical (index4, index6)
    expect_silent (index7 <- match_pts_to_graph (net, pts, connected = TRUE))
    expect_identical (index4, index7)

})

test_that ("add nodes to graph", {

    graph0 <- weight_streetnet (hampi, wt_profile = "foot")
    verts <- dodgr_vertices (graph0)
    set.seed (1)
    npts <- 10
    xy <- data.frame (
        x = min (verts$x) + runif (npts) * diff (range (verts$x)),
        y = min (verts$y) + runif (npts) * diff (range (verts$y))
    )

    graph1 <- add_nodes_to_graph (graph0, xy)

    expect_identical (colnames (graph0), colnames (graph1))
    expect_true ((nrow (graph1) - nrow (graph0)) > npts)
    # actually equals 2 * npts when all edges are bi-directional.
})
