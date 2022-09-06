context ("dodgr graph functions")

test_all <- (identical (Sys.getenv ("MPADGE_LOCAL"), "true") |
    identical (Sys.getenv ("GITHUB_WORKFLOW"), "test-coverage"))

skip_if (!test_all)

dodgr_cache_off ()
clear_dodgr_cache ()

test_that ("sample graph", {

    graph <- weight_streetnet (hampi)
    set.seed (1)
    nverts <- 100
    graph_s <- dodgr_sample (graph, nverts = nverts)
    expect_true (nrow (graph_s) < nrow (graph))
    v <- dodgr_vertices (graph_s)
    expect_true (nrow (v) == nverts)

    # that sample is only of largest component, so the subsequent code removing
    # component should generate longer distances
    d <- mean (geodist::geodist (v))

    graph$component <- NULL
    set.seed (1)
    graph_s2 <- dodgr_sample (graph, nverts = nverts)
    expect_true (nrow (graph_s2) < nrow (graph))
    v <- dodgr_vertices (graph_s2)
    expect_true (nrow (v) == nverts)

    d2 <- mean (geodist::geodist (v))
    # expect_true (d2 > d) # that's not reliably true, but almost always

    graph <- weight_streetnet (hampi)
    graph$edge_id <- NULL
    expect_silent (graphs <- dodgr_sample (graph, nverts = nverts))
    expect_is (graphs$edge_id, "integer")
    expect_true (min (graphs$edge_id) >= 1)
    expect_true (max (graphs$edge_id) <= nrow (graph))
})

test_that ("insert_vertex", {
    graph <- weight_streetnet (hampi)
    e1 <- 2256
    v1 <- graph$from_id [e1]
    v2 <- graph$to_id [e1]
    expect_silent (graph2 <- dodgr_insert_vertex (graph, v1 = v1, v2 = v2))
    # graph should have two more rows added:
    expect_equal (nrow (graph2) - 2, nrow (graph))
})

test_that ("components", {
    graph <- weight_streetnet (hampi)
    comp <- graph$component
    graph$component <- NULL
    expect_silent (graph <- dodgr_components (graph))
    expect_identical (comp, graph$component)
    comp <- graph$component
    expect_identical (comp, graph$component)

    expect_message (
        graph2 <- dodgr_components (graph),
        "graph already has a component column"
    )
    expect_identical (graph, graph2)

    graph$edge_id <- NULL
    expect_message (
        graph3 <- dodgr_components (graph),
        "graph already has a component column"
    )
    expect_identical (graph2$component, graph3$component)

    expect_silent (clear_dodgr_cache ())
    expect_message (
        graph4 <- dodgr_components (graph),
        "graph already has a component column"
    )
    expect_identical (graph3, graph4)

    expect_identical (graph2$component, graph4$component)

    graph$component <- NULL
    expect_silent (graph4 <- dodgr_components (graph))
    expect_identical (graph4$component, graph2$component)
})

test_that ("contract graph", {
    graph <- weight_streetnet (hampi)
    expect_silent (graph_c <- dodgr_contract_graph (graph))
    expect_true (nrow (graph_c) < nrow (graph))

    vc <- dodgr_vertices (graph_c)
    v <- dodgr_vertices (graph)
    verts <- sample (v$id [which (!v$id %in% vc$id)], size = 10)
    expect_silent (graph_c2 <- dodgr_contract_graph (graph, verts = verts))
    expect_true (nrow (graph_c2) > nrow (graph_c))

    verts <- as.matrix (verts, ncol = 1)
    expect_error (
        graph_c3 <- dodgr_contract_graph (graph, verts = verts),
        "verts must be a single value or a vector of vertex IDs"
    )

    verts <- as.numeric (verts [, 1])
    expect_silent (graph_c4 <- dodgr_contract_graph (graph, verts = verts))
    expect_identical (graph_c2, graph_c4)
})

test_that ("uncontract graph", {
    clear_dodgr_cache ()
    graph <- weight_streetnet (hampi)
    graph_c <- dodgr_contract_graph (graph)
    graph2 <- dodgr_uncontract_graph (graph_c)
    expect_identical (dim (graph), dim (graph2))
    expect_identical (graph$edge_id, graph2$edge_id)

    # dodgr_contract_graph in that case just calls the cached version. This
    # checks re-contraction:
    graph$edge_id <- seq (nrow (graph))
    graph_c <- dodgr_contract_graph (graph)
    graph2 <- dodgr_uncontract_graph (graph_c)
    expect_identical (dim (graph), dim (graph2))
    expect_identical (graph$edge_id, graph2$edge_id)

    graph_c$edge_id <- seq (nrow (graph_c))
    graph2 <- dodgr_uncontract_graph (graph_c)
    # with no edge ids, graph uncontraction is not possible:
    expect_equal (nrow (graph2), 0L)
})

test_that ("compare heaps", {
    graph <- weight_streetnet (hampi)
    ch <- compare_heaps (graph, nverts = 100)
    expect_equal (nrow (ch), 11L)
    # Test that all dodgr calculations are faster than igraph:
    igr <- which (grepl ("igraph", ch$expression))
    # expect_true (ch$elapsed [igr] == max (ch$elapsed))
    # This actually fails on some machines (R oldrel on Windows) because elapsed
    # times are sometimes all very small *and equal*, so is turned off:
    # expect_true (ch$elapsed [igr] > min (ch$elapsed))
})

test_that ("dodgr2sf", {
    hw <- weight_streetnet (hampi)
    y <- dodgr_to_sfc (hw)
    # y should have more linestrings than the original sf object:
    expect_true (length (y) > length (hw$geometry))
})

test_that ("different geometry columns", {
    h2 <- hampi
    gcol <- grep ("geometry", names (h2))
    names (h2) [gcol] <- "g"
    attr (h2, "sf_column") <- "g" # not necessary here but should always be done
    expect_silent (net <- weight_streetnet (h2))
    h2 <- hampi
    names (h2) [gcol] <- "geoms"
    attr (h2, "sf_column") <- "geoms"
    expect_silent (net <- weight_streetnet (h2))
    h2 <- hampi
    names (h2) [gcol] <- "geometry"
    attr (h2, "sf_column") <- "geometry"
    expect_silent (net <- weight_streetnet (h2))
    h2 <- hampi
    names (h2) [gcol] <- "xxx"
    expect_error (
        net <- weight_streetnet (h2),
        "Unable to determine geometry column"
    )

    h2 <- data.frame (hampi) # remove sf class
    expect_error (net <- weight_streetnet (h2), "Unknown class")
})

test_that ("no geom rownames", {
    hw0 <- weight_streetnet (hampi)
    g0 <- hampi$geometry
    attr (g0, "names") <- NULL # remove way IDs
    for (i in seq (g0)) {
        rownames (g0 [[i]]) <- NULL
    } # remove all node IDs
    h2 <- hampi
    h2$geometry <- g0
    hw1 <- weight_streetnet (h2)
    expect_true (!identical (hw0, hw1))
    expect_equal (ncol (hw0), ncol (hw1))
    expect_true (!identical (hw0$from_id, hw1$from_id))
    expect_true (!identical (hw0$to_id, hw1$to_id))

    indx0 <- which (!names (hw0) %in% c (
        "from_id", "to_id", "way_id",
        "component"
    ))
    indx1 <- which (!names (hw1) %in% c (
        "from_id", "to_id", "way_id",
        "component"
    ))
    expect_identical (hw0 [, indx0], hw1 [, indx1])
    # components are not identical because ones of equal size are assigned
    # random numbers, but all other columns remain identical:
    # geom_num, edge_id, lon/lat values, d, d_weighted, and highway
})

test_that ("keep cols", {
    hw0 <- weight_streetnet (hampi)
    expect_equal (ncol (hw0), 15)
    hw1 <- weight_streetnet (hampi, keep_cols = "foot")
    expect_equal (ncol (hw1), 16)
    expect_true ("foot" %in% names (hw1))
    expect_false ("foot" %in% names (hw0))

    i <- which (names (hampi) == "foot")
    hw2 <- weight_streetnet (hampi, keep_cols = i)
    attr (hw1, "px") <- attr (hw2, "px") <- NULL
    expect_identical (hw1, hw2)

    expect_error (
        hw3 <- weight_streetnet (hampi, keep_cols = list (i)),
        "keep_cols must be either character or numeric"
    )
})

test_that ("graph columns", {
    graph <- data.frame (weight_streetnet (hampi)) # rm dodgr_streetnet class
    nf <- 100
    nt <- 50
    from <- sample (graph$from_id, size = nf)
    to <- sample (graph$from_id, size = nt)
    d <- dodgr_dists (graph, from = from, to = to)
    graph$from_id2 <- graph$from_id
    expect_error (
        d <- dodgr_dists (graph, from = from, to = to),
        "Unable to determine column with ID of from vertices"
    )
    graph$from_id2 <- NULL
    graph$to_id2 <- graph$to_id
    expect_error (
        d <- dodgr_dists (graph, from = from, to = to),
        "Unable to determine column with ID of to vertices"
    )

    graph <- data.frame (weight_streetnet (hampi)) # rm dodgr_streetnet class
    graph$from_lat <- NULL
    expect_error (
        d <- dodgr_dists (graph, from = from, to = to),
        "Unable to determine coordinate columns of graph"
    )

    graph <- data.frame (weight_streetnet (hampi)) # rm dodgr_streetnet class
    graph$d2 <- graph$d
    expect_error (
        d <- dodgr_dists (graph, from = from, to = to),
        "Unable to determine distance column in graph"
    )
    graph$d2 <- NULL
    graph$d_wt <- graph$d_weighted
    expect_error (
        d <- dodgr_dists (graph, from = from, to = to),
        "Unable to determine weight column in graph"
    )

    graph <- data.frame (weight_streetnet (hampi)) # rm dodgr_streetnet class
    graph$from_lon <- paste0 (graph$from_lon)
    expect_error (
        d <- dodgr_dists (graph, from = from, to = to),
        "graph appears to have non-numeric longitudes and latitudes"
    )

    graph <- data.frame (weight_streetnet (hampi)) # rm dodgr_streetnet class
    class (graph) <- c (class (graph), "tbl")
    expect_silent (d <- dodgr_dists (graph, from = from, to = to))
})

test_that ("get_id_cols", {
    n <- 10
    pts <- cbind (runif (n), runif (n))
    expect_null (ids <- get_id_cols (pts))
    rownames (pts) <- seq (n)
    expect_is (get_id_cols (pts), "character")
    expect_length (get_id_cols (pts), n)

    pts <- runif (n)
    expect_null (get_id_cols (pts))
    names (pts) <- seq (n)
    expect_is (get_id_cols (pts), "character")
    expect_length (get_id_cols (pts), n)

    pts <- cbind (runif (n), runif (n), seq (n))
    expect_null (get_id_cols (pts))
    colnames (pts) <- c ("a", "b", "id")
    expect_is (get_id_cols (pts), "numeric")
    expect_length (get_id_cols (pts), n)

    pts <- data.frame (pts)
    expect_is (get_id_cols (pts), "numeric")
    expect_length (get_id_cols (pts), n)

    names (pts) [3] <- "this_is_an_id_here" # "id" is grepped
    expect_is (get_id_cols (pts), "numeric")
    expect_length (get_id_cols (pts), n)
})

test_that ("get_pts_index", {
    graph <- weight_streetnet (hampi)
    gr_cols <- dodgr_graph_cols (graph)
    vert_map <- make_vert_map (graph, gr_cols)

    n <- 10
    pts <- cbind (runif (n), runif (n))
    expect_error (
        get_pts_index (graph, gr_cols, vert_map, pts),
        "Unable to determine geographical coordinates in from/to"
    )
    rownames (pts) <- seq (n)
    expect_error (
        get_pts_index (graph, gr_cols, vert_map, pts),
        "Unable to determine geographical coordinates in from/to"
    )

    pts <- cbind (runif (n), runif (n), seq (n))
    colnames (pts) <- c ("a", "b", "id")
    expect_error (
        get_pts_index (graph, gr_cols, vert_map, pts),
        "Unable to determine geographical coordinates in from/to"
    )
    colnames (pts) <- c ("x", "y", "id")
    expect_is (get_pts_index (graph, gr_cols, vert_map, pts), "numeric")
    expect_length (get_pts_index (graph, gr_cols, vert_map, pts), n)

    pts <- data.frame (pts)
    expect_is (get_pts_index (graph, gr_cols, vert_map, pts), "numeric")
    expect_length (get_pts_index (graph, gr_cols, vert_map, pts), n)
})
