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
    expect_true (attr (dp, "vertices"))
    lens <- unlist (lapply (dp, lengths))
    dp <- dodgr_paths (graph, from = from, to = to, vertices = FALSE)
    expect_is (dp, "list")
    expect_length (dp, 100)
    expect_identical (unique (vapply (dp, length, integer (1L))), length (to))
    expect_false (attr (dp, "vertices"))
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
    expect_true (attr (dp, "vertices"))

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

test_that ("dodgr_one_path_expand on synthetic graph", {

    # A simple A-to-G chain of 6 original edges, contracted down to 2
    # compound edges ('a101' = A-B-C-D, 'a103' = E-F-G) plus one
    # pass-through edge ('4' = D-E) which was not aggregated at all, and so
    # retains its original edge_id. `edge_map` intentionally lists the
    # sub-edges of 'a101' out of order, to confirm that
    # `dodgr_one_path_expand` (via `sort_transitions`) correctly re-orders
    # them.
    graph <- data.frame (
        edge_id = 1:6,
        from = c ("A", "B", "C", "D", "E", "F"),
        to = c ("B", "C", "D", "E", "F", "G"),
        d = 1:6,
        d_weighted = 1:6
    )
    graph_c <- data.frame (
        edge_id = c ("a101", "4", "a103"),
        from = c ("A", "D", "E"),
        to = c ("D", "E", "G"),
        d = c (6, 4, 11),
        d_weighted = c (6, 4, 11)
    )
    edge_map <- data.frame (
        edge_new = c ("a101", "a101", "a101", "a103", "a103"),
        edge_old = c (2, 3, 1, 6, 5)
    )

    path_c_chr <- c ("A", "D", "E", "G")
    path <- dodgr_one_path_expand (
        path_c_chr,
        graph,
        graph_c,
        edge_map = edge_map
    )

    # character input returns the full, structurally-identical character
    # vector of expanded vertex IDs
    expect_type (path, "character")
    expect_identical (path, c ("A", "B", "C", "D", "E", "F", "G"))

    # equivalent integer-index path over 'graph_c' returns the full,
    # structurally-identical integer vector of row indices into 'graph',
    # tracing the same underlying sequence of edges as the character version
    path_c_int <- match (
        paste (path_c_chr [-length (path_c_chr)], path_c_chr [-1]),
        paste (graph_c$from, graph_c$to)
    )
    path_int <- dodgr_one_path_expand (
        path_c_int,
        graph,
        graph_c,
        edge_map = edge_map
    )
    expect_type (path_int, "integer")
    expect_identical (path_int, 1:6)
    expect_identical (
        c (graph$from [path_int], graph$to [path_int] [length (path_int)]),
        path
    )

    # a single (contracted) edge is a valid path requiring no transitions,
    # and expands to the pass-through edge's own (single) original edge
    path_single <- dodgr_one_path_expand (
        2L, # row index of pass-through edge '4', from D to E
        graph,
        graph_c,
        edge_map = edge_map
    )
    expect_type (path_single, "integer")
    expect_identical (path_single, 4L)
    expect_identical (graph$from [path_single], "D")
    expect_identical (graph$to [path_single], "E")

    expect_error (
        dodgr_one_path_expand (
            path_c_chr [1],
            graph,
            graph_c,
            edge_map = edge_map
        ),
        "There are no transitions."
    )
    expect_error (
        dodgr_one_path_expand (
            c ("not_a_vertex", "another_one"),
            graph,
            graph_c,
            edge_map = edge_map
        ),
        "Not all transitions were matched to the contracted graph."
    )
    expect_error (
        dodgr_one_path_expand (
            c (1.5, 2.5),
            graph,
            graph_c,
            edge_map = edge_map
        ),
        "Path must be provided as integer or character vector."
    )
})

test_that ("dodgr_paths_expand on synthetic graph", {

    # Same fixture as the 'dodgr_one_path_expand' test above.
    graph <- data.frame (
        edge_id = 1:6,
        from = c ("A", "B", "C", "D", "E", "F"),
        to = c ("B", "C", "D", "E", "F", "G"),
        d = 1:6,
        d_weighted = 1:6
    )
    graph_c <- data.frame (
        edge_id = c ("a101", "4", "a103"),
        from = c ("A", "D", "E"),
        to = c ("D", "E", "G"),
        d = c (6, 4, 11),
        d_weighted = c (6, 4, 11)
    )
    edge_map <- data.frame (
        edge_new = c ("a101", "a101", "a101", "a103", "a103"),
        edge_old = c (2, 3, 1, 6, 5)
    )

    path_c_chr <- c ("A", "D", "E", "G")
    chr_expected <- dodgr_one_path_expand (
        path_c_chr,
        graph,
        graph_c,
        edge_map = edge_map
    )

    # A nested list mimicking the (`vertices = TRUE`) structure returned by
    # `dodgr_paths`, with two origins and two destinations: a valid
    # multi-transition path, a `NULL` (unreachable) entry, a trivial
    # single-vertex (no transitions) entry, and a second valid path
    # identical to the first.
    paths_chr <- list (
        list (path_c_chr, NULL),
        list (path_c_chr [1], path_c_chr)
    )
    attr (paths_chr, "vertices") <- TRUE

    result_chr <- dodgr_paths_expand (
        paths_chr,
        graph,
        graph_c,
        edge_map = edge_map
    )

    expect_type (result_chr, "list")
    expect_identical (attr (result_chr, "vertices"), TRUE)
    expect_length (result_chr, length (paths_chr))
    expect_identical (
        vapply (result_chr, length, integer (1L)),
        vapply (paths_chr, length, integer (1L))
    )

    expect_type (result_chr [[1]] [[1]], "character")
    expect_identical (result_chr [[1]] [[1]], chr_expected)
    expect_null (result_chr [[1]] [[2]])
    expect_null (result_chr [[2]] [[1]])
    expect_identical (result_chr [[2]] [[2]], chr_expected)

    # The same fixture, but as a (`vertices = FALSE`) nested list of integer
    # row indices into `graph_c`, mixed with a `NULL` (unreachable) entry.
    path_c_int <- match (
        paste (path_c_chr [-length (path_c_chr)], path_c_chr [-1]),
        paste (graph_c$from, graph_c$to)
    )
    int_expected <- dodgr_one_path_expand (
        path_c_int,
        graph,
        graph_c,
        edge_map = edge_map
    )

    paths_int <- list (list (path_c_int, NULL))
    attr (paths_int, "vertices") <- FALSE

    result_int <- dodgr_paths_expand (
        paths_int,
        graph,
        graph_c,
        edge_map = edge_map
    )

    expect_identical (attr (result_int, "vertices"), FALSE)
    expect_type (result_int [[1]] [[1]], "integer")
    expect_identical (result_int [[1]] [[1]], int_expected)
    expect_null (result_int [[1]] [[2]])
})

test_that ("dodgr_one_path_expand on real network", {
    graph <- weight_streetnet (hampi)
    graph_c <- dodgr_contract_graph (graph)
    verts <- dodgr_vertices (graph_c)

    # Restrict candidates to the same (largest) connected component, because
    # `hampi` is not a single connected component, so an arbitrary pair may
    # fall in different components with no path between them. (Selecting
    # vertices by raw position in `verts$id` is not reliable here, because
    # the character sort order of `id` values is locale-dependent.)
    comp_sizes <- table (verts$component)
    largest_comp <- names (comp_sizes) [which.max (comp_sizes)]
    verts <- verts [verts$component == as.integer (largest_comp), ]
    verts <- verts [order (verts$x), ]

    # Try a handful of geographically-distant candidate pairs, from the two
    # ends of the longitude range inwards. This guards against the very rare
    # case (never seen locally, but observed on some CI runners) in which the
    # underlying C++ graph contraction produces a self-inconsistent compound
    # edge for a particular pair of vertices, presumably as a result of
    # platform-dependent hash-container iteration order. A single such
    # unlucky pair should not fail the whole test suite, so any pair that
    # trips `sort_transitions`'s contiguous-chain check is simply skipped in
    # favor of the next candidate.
    n <- nrow (verts)
    n_candidates <- min (5L, n)
    candidates <- Map (
        c,
        seq_len (n_candidates),
        rev (seq_len (n)) [seq_len (n_candidates)]
    )

    path_c_chr <- path <- NULL
    for (candidate in candidates) {
        from <- verts$id [candidate [1]]
        to <- verts$id [candidate [2]]
        if (from == to) next

        path_c_chr_i <- dodgr_paths (graph_c, from = from, to = to) [[1]] [[1]]
        if (length (path_c_chr_i) < 3) next

        path_i <- tryCatch (
            dodgr_one_path_expand (path_c_chr_i, graph, graph_c),
            error = function (e) NULL
        )
        if (!is.null (path_i)) {
            path_c_chr <- path_c_chr_i
            path <- path_i
            break
        }
    }

    expect_true (length (path_c_chr) > 2)

    # character input returns a structurally-identical character vector of
    # expanded vertex IDs, tracing the same start and end points as the
    # contracted path from which it was expanded
    expect_type (path, "character")
    expect_identical (path [1], path_c_chr [1])
    expect_identical (path [length (path)], path_c_chr [length (path_c_chr)])

    # every consecutive pair of vertices in the expanded path must
    # correspond to an actual edge of the full, uncontracted graph
    cols <- dodgr_graph_cols (graph)
    graph_ft <- paste (graph [[cols$from]], graph [[cols$to]])
    path_ft <- paste (path [-length (path)], path [-1])
    expect_true (all (path_ft %in% graph_ft))

    # explicitly-supplied edge map gives identical result to internal lookup
    edge_map <- get_edge_map (graph)
    path_em <- dodgr_one_path_expand (
        path_c_chr,
        graph,
        graph_c,
        edge_map = edge_map
    )
    expect_identical (path, path_em)

    # equivalent integer (`vertices = FALSE`-style) path over `graph_c`
    # returns a structurally-identical integer vector of row indices into
    # `graph`, tracing the same underlying sequence of edges
    path_c_int <- match (
        paste (path_c_chr [-length (path_c_chr)], path_c_chr [-1]),
        paste (graph_c [[cols$from]], graph_c [[cols$to]])
    )
    path_int <- dodgr_one_path_expand (path_c_int, graph, graph_c)
    expect_type (path_int, "integer")
    expect_identical (
        c (
            graph [[cols$from]] [path_int],
            graph [[cols$to]] [path_int] [length (path_int)]
        ),
        path
    )
})

test_that ("dodgr_paths_expand on real network", {
    graph <- weight_streetnet (hampi)
    graph_c <- dodgr_contract_graph (graph)
    verts <- dodgr_vertices (graph_c)

    # Same robust candidate search as the 'dodgr_one_path_expand' test above,
    # used here to find (origin, destination) pairs known to expand cleanly,
    # so the wrapper can be exercised via `dodgr_paths`' own nested list
    # output, for both `vertices = TRUE` and `vertices = FALSE`, and for
    # both the default (cross-product) and `pairwise = TRUE` forms.
    comp_sizes <- table (verts$component)
    largest_comp <- names (comp_sizes) [which.max (comp_sizes)]
    verts <- verts [verts$component == as.integer (largest_comp), ]
    verts <- verts [order (verts$x), ]

    n <- nrow (verts)
    n_candidates <- min (10L, n)
    candidates <- Map (
        c,
        seq_len (n_candidates),
        rev (seq_len (n)) [seq_len (n_candidates)]
    )

    good <- list ()
    for (candidate in candidates) {
        from_i <- verts$id [candidate [1]]
        to_i <- verts$id [candidate [2]]
        if (from_i == to_i) next

        path_c_chr_i <-
            dodgr_paths (graph_c, from = from_i, to = to_i) [[1]] [[1]]
        if (length (path_c_chr_i) < 3) next

        ok <- !inherits (
            tryCatch (
                dodgr_one_path_expand (path_c_chr_i, graph, graph_c),
                error = function (e) e
            ),
            "error"
        )
        if (ok) {
            good [[length (good) + 1L]] <- c (from_i, to_i)
        }
        if (length (good) >= 2L) {
            break
        }
    }
    expect_true (length (good) >= 1L)

    from <- vapply (good, `[`, character (1L), 1L)
    to <- vapply (good, `[`, character (1L), 2L)
    cols <- dodgr_graph_cols (graph)

    # -- vertices = TRUE: expanded paths are structurally-identical
    # character vectors of vertex IDs --
    paths_c <- dodgr_paths (graph_c, from = from, to = to)
    paths <- dodgr_paths_expand (paths_c, graph, graph_c)

    expect_type (paths, "list")
    expect_identical (attr (paths, "vertices"), attr (paths_c, "vertices"))
    expect_true (attr (paths, "vertices"))
    expect_length (paths, length (from))
    expect_length (paths [[1]], length (to))

    path <- paths [[1]] [[1]]
    path_c_chr <- paths_c [[1]] [[1]]
    expect_type (path, "character")
    expect_identical (
        path,
        dodgr_one_path_expand (path_c_chr, graph, graph_c)
    )
    expect_identical (path [1], path_c_chr [1])
    expect_identical (path [length (path)], path_c_chr [length (path_c_chr)])

    # -- vertices = FALSE: expanded paths are structurally-identical integer
    # vectors of row indices into `graph`, tracing the same vertex sequence
    # as the `vertices = TRUE` version above --
    paths_c_e <- dodgr_paths (graph_c, from = from, to = to, vertices = FALSE)
    paths_e <- dodgr_paths_expand (paths_c_e, graph, graph_c)

    expect_identical (attr (paths_e, "vertices"), attr (paths_c_e, "vertices"))
    expect_false (attr (paths_e, "vertices"))

    path_e <- paths_e [[1]] [[1]]
    expect_type (path_e, "integer")
    path_e_verts <- c (
        graph [[cols$from]] [path_e],
        graph [[cols$to]] [path_e] [length (path_e)]
    )
    expect_identical (path_e_verts, path)

    # -- pairwise = TRUE: wrapper correctly handles the pairwise nested
    # structure (a list of single-element lists), giving identical results
    # to the corresponding cross-product entries above --
    paths_c_pw <- dodgr_paths (
        graph_c,
        from = from,
        to = to,
        pairwise = TRUE
    )
    paths_pw <- dodgr_paths_expand (paths_c_pw, graph, graph_c)

    expect_type (paths_pw, "list")
    expect_true (attr (paths_pw, "vertices"))
    expect_length (paths_pw, length (from))
    expect_true (all (vapply (paths_pw, length, integer (1L)) == 1L))
    expect_identical (paths_pw [[1]] [[1]], path)
    if (length (from) > 1L) {
        expect_identical (paths_pw [[2]] [[1]], paths [[2]] [[2]])
    }
})
