test_all <- (identical (Sys.getenv ("MPADGE_LOCAL"), "true") ||
    identical (Sys.getenv ("GITHUB_JOB"), "test-coverage"))

if (!test_all) {
    RcppParallel::setThreadOptions (numThreads = 2)
}

# A vertex may only be closed when it is taken from the heap. In this graph, 'x'
# is first reached directly from 'a' at cost 10. 'b' (cost 1) is then taken from
# the heap, and its edge to 'x' (total 21) does not improve on 10. 'c' (cost 2)
# is taken next, and its edge to 'x' (total 3) does. The shortest path from 'a'
# to 'x' is therefore a -> c -> x. Closing 'x' after the failed relaxation from
# 'b' used to discard that improvement.
mre_graph <- function () {
    data.frame (
        from = c ("a", "a", "a", "b", "c"),
        to = c ("x", "b", "c", "x", "x"),
        d = c (10, 1, 2, 20, 1)
    )
}

# Convert to a spatial graph, which is routed with A*. Vertices are placed so
# close together that the A* heuristic is always much smaller than any edge, so
# does not alter the order in which vertices are taken from the heap.
spatial_graph <- function (graph) {
    verts <- unique (c (graph$from, graph$to))
    xy <- 1e-5 * seq_along (verts)
    i <- match (graph$from, verts)
    j <- match (graph$to, verts)
    data.frame (
        from_id = graph$from,
        to_id = graph$to,
        from_lon = xy [i],
        from_lat = xy [i],
        to_lon = xy [j],
        to_lat = xy [j],
        d = graph$d,
        d_weighted = graph$d
    )
}

test_that ("dists after a non-improving relaxation", {

    graph <- mre_graph ()
    for (heap in c ("BHeap", "FHeap", "TriHeap", "Heap23", "set")) {
        d <- dodgr_dists (graph, from = "a", to = "x", heap = heap)
        expect_equal (d [1, 1], 3, info = heap)
    }

    # A* on spatial graphs, both for all 'to' points, and pairwise:
    graph <- spatial_graph (graph)
    d <- dodgr_dists (graph, from = "a", to = "x")
    expect_equal (d [1, 1], 3)
    d <- dodgr_dists (graph, from = "a", to = "x", pairwise = TRUE)
    expect_equal (d [1, 1], 3)
})

test_that ("paths after a non-improving relaxation", {

    graph <- mre_graph ()
    p <- dodgr_paths (graph, from = "a", to = "x")
    expect_identical (p [[1]] [[1]], c ("a", "c", "x"))

    graph <- spatial_graph (graph)
    p <- dodgr_paths (graph, from = "a", to = "x", pairwise = TRUE)
    expect_identical (p [[1]] [[1]], c ("a", "c", "x"))
})

test_that ("flows after a non-improving relaxation", {

    testthat::skip_on_cran ()

    graph <- mre_graph ()
    flows <- dodgr_flows_aggregate (
        graph,
        from = "a", to = "x", flows = matrix (1),
        contract = FALSE, norm_sums = FALSE
    )
    # All flow must go along a -> c -> x, which are rows 3 and 5:
    expect_equal (flows$flow, c (0, 0, 1, 0, 1))
})

test_that ("bidirectional A* after a non-improving relaxation", {

    # Pairwise dists and paths on uncontracted graphs of more than 1000 vertices
    # use bidirectional A*. The search backwards from 'a' here meets the same
    # pattern as above, with 'x' reached first directly from 'a'. The search
    # forwards from 's' is held back by the expensive first edge, and so can not
    # make up for a backwards search which closed 'x' too early. The ring of
    # vertices is disconnected from the rest, and only brings the number of
    # vertices above 1000.
    graph <- data.frame (
        from = c ("s", "x", "b", "c", "x", "x"),
        to = c ("x", "a", "a", "a", "b", "c"),
        d = c (100, 10, 1, 2, 20, 1)
    )
    ring <- paste0 ("r", 1:1100)
    graph <- rbind (
        graph,
        data.frame (from = ring, to = ring [c (2:1100, 1)], d = 1)
    )
    graph <- spatial_graph (graph)
    expect_gt (nrow (dodgr_vertices (graph)), 1000L)

    d <- dodgr_dists (graph, from = "s", to = "a", pairwise = TRUE)
    expect_equal (d [1, 1], 103)
    p <- dodgr_paths (graph, from = "s", to = "a", pairwise = TRUE)
    expect_identical (p [[1]] [[1]], c ("s", "x", "c", "a"))
})

test_that ("categorical dists after a non-improving relaxation", {

    # The first five edges are the graph from above; edges from 'b' and 'x' are
    # added so that these are also explored by the 'dlimit' search below.
    graph <- data.frame (
        from = c ("a", "a", "a", "b", "c", "b", "x", "x"),
        to = c ("x", "b", "c", "x", "x", "y", "z", "w"),
        d = c (10, 1, 2, 20, 1, 0.5, 5, 0.4)
    )
    graph <- spatial_graph (graph)
    graph$edge_type <- c ("p", "p", "q", "p", "q", "p", "q", "p")

    # a -> c -> x is along edges of type "q" only
    d <- dodgr_dists_categorical (graph, from = "a", to = "x", pairwise = TRUE)
    expect_equal (d [1, "total"], 3)
    expect_equal (d [1, "p"], 0)
    expect_equal (d [1, "q"], 3)

    # 'dlimit' aggregates distances to the first vertices beyond the limit; here
    # that is 'z', reached along a -> c -> x -> z, again along type "q" only.
    d <- dodgr_dists_categorical (graph, from = "a", dlimit = 3.5)
    expect_equal (d$distance, 8)
    expect_equal (d$p, 0)
    expect_equal (d$q, 8)
})
