test_that ("turn penalty junction without incoming edges errors", {

    edges <- data.frame (
        .vx0 = c ("a", "a", "a", "a", "f"),
        .vx1 = c ("b", "c", "d", "e", "g"),
        .vx0_x = c (4.4760821, 4.4760821, 4.4760821, 4.4760821, 4.4774875),
        .vx0_y = c (51.8901285, 51.8901285, 51.8901285, 51.8901285, 51.8911503),
        .vx1_x = c (4.4749388, 4.4760279, 4.4761105, 4.4761437, 4.4773415),
        .vx1_y = c (51.8901787, 51.8900749, 51.890205, 51.8901351, 51.89111),
        d = c (78.5, 7.0, 8.7, 4.3, 11.0),
        d_weighted = c (157.1, 14.0, 17.4, 8.6, 21.9),
        time = c (9.4, 0.8, 1.0, 0.5, 1.3),
        time_weighted = c (18.8, 1.7, 2.1, 1.0, 2.6),
        stringsAsFactors = FALSE
    )
    edges$edge_ <- paste0 ("e", seq_len (nrow (edges)))
    edges$object_ <- edges$edge_
    edges$highway <- "residential"
    edges$component <- 1L
    attr (edges, "turn_penalty") <- 1
    attr (edges, "wt_profile") <- "motorcar"
    attr (edges, "left_side") <- FALSE
    class (edges) <- c ("dodgr_streetnet", "dodgr_streetnet_sc", "data.frame")

    expect_error (
        dodgr_dists (edges, from = "f", to = "a"),
        "No incoming edges for vertices: a"
    )
    expect_silent (d <- dodgr_dists (edges, from = "a", to = c ("b", "c")))
    expect_false (anyNA (d))
})
