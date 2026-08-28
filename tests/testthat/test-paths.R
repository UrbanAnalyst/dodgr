context ("dodgr_paths")

test_that ("paths", {
    graph <- weight_streetnet (hampi)
    from <- graph$from_id [1:100]
    to <- graph$to_id [100:150]
    to <- to [!to %in% from]
    expect_message (
        dp <- dodgr_paths (graph,
            from = from,
            to = to,
            quiet = FALSE
        ),
        "Calculating shortest paths ..."
    )
    expect_is (dp, "list")
    expect_length (dp, 100)
    expect_identical (unique (vapply (dp, length, integer (1L))), length (to))
    expect_is (dp [[1]] [[1]], "character")
    lens <- unlist (lapply (dp, lengths))
    dp <- dodgr_paths (graph, from = from, to = to, vertices = FALSE)
    expect_is (dp, "list")
    expect_length (dp, 100)
    expect_identical (unique (vapply (dp, length, integer (1L))), length (to))
    lens2 <- unlist (lapply (dp, lengths))
    # edge lists should all have one less item than vertex lists
    lens2 <- lens2 [which (lens > 0)]
    lens <- lens [which (lens > 0)]
    expect_true (all (abs (lens - lens2) == 1))
})

test_that ("pairwise paths", {
    graph <- weight_streetnet (hampi)
    from <- graph$from_id [1:10]
    to <- graph$to_id [100:105]
    indx <- which (!to %in% from)
    to <- to [indx]
    from <- from [indx]
    n <- length (indx)
    dp <- dodgr_paths (graph, from = from, to = to, pairwise = TRUE)
    expect_is (dp, "list")
    expect_length (dp, n)
    expect_true (all (lapply (dp, length) == 1))

    expect_error (
        dp <- dodgr_paths (graph,
            from = from, to = to [-1],
            pairwise = TRUE
        ),
        "pairwise paths require from and to to have same length"
    )
})

test_that ("sort_transitions", {
    graph <- data.frame (
        edge_id = 1:4,
        from = c ("a", "c", "b", "d"),
        to = c ("b", "d", "c", "e"),
        d = 1:4,
        d_weighted = 1:4
    )
    perm <- sort_transitions (graph)
    expect_identical (perm, c (1L, 3L, 2L, 4L))
    expect_identical (graph$from [perm], c ("a", "b", "c", "d"))
    expect_identical (graph$to [perm], c ("b", "c", "d", "e"))
})

test_that ("dodgr_paths_expand", {
    graph <- weight_streetnet (hampi)
    graph_c <- dodgr_contract_graph (graph)
    verts <- dodgr_vertices (graph_c)
    from <- verts$id [1]
    to <- verts$id [50]

    path_c_chr <- dodgr_paths (graph_c, from = from, to = to) [[1]] [[1]]
    expect_true (length (path_c_chr) > 2)

    path <- dodgr_paths_expand (path_c_chr, graph, graph_c)
    expect_s3_class (path, "data.frame")
    expect_identical (names (path), names (graph))

    # expanded path must be contiguous, and trace the same start and end
    # points as the contracted path from which it was expanded
    expect_identical (path$to_id [-nrow (path)], path$from_id [-1])
    expect_identical (path$from_id [1], path_c_chr [1])
    expect_identical (path$to_id [nrow (path)], path_c_chr [length (path_c_chr)])

    # equivalent integer-index path over 'graph_c' gives identical result
    path_c_int <- match (
        paste (
            path_c_chr [-length (path_c_chr)],
            path_c_chr [-1]
        ),
        paste (graph_c$from_id, graph_c$to_id)
    )
    path_int <- dodgr_paths_expand (path_c_int, graph, graph_c)
    expect_identical (path$edge_id, path_int$edge_id)

    # explicitly-supplied edge map gives identical result to internal lookup
    edge_map <- get_edge_map (graph)
    path_em <- dodgr_paths_expand (
        path_c_chr,
        graph,
        graph_c,
        edge_map = edge_map
    )
    expect_identical (path$edge_id, path_em$edge_id)

    expect_error (
        dodgr_paths_expand (path_c_chr [1], graph, graph_c),
        "There are no transitions."
    )
    expect_error (
        dodgr_paths_expand (c ("not_a_vertex", "another_one"), graph, graph_c),
        "Not all transitions were matched to the contracted graph."
    )
    expect_error (
        dodgr_paths_expand (c (1.5, 2.5), graph, graph_c),
        "Path must be provided as integer or character vector."
    )
})
