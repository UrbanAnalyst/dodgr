context("dodgr graph functions")

test_all <- (identical (Sys.getenv ("MPADGE_LOCAL"), "true") |
             identical (Sys.getenv ("TRAVIS"), "true"))

test_that("sample graph", {
    graph <- weight_streetnet (hampi)
    graph_s <- dodgr_sample (graph, nverts = 100)
    expect_true (nrow (graph_s) < nrow (graph))
    v <- dodgr_vertices (graph_s)
    expect_true (nrow (v) == 100)
})

test_that("components", {
    graph <- weight_streetnet (hampi)
    comp <- graph$component
    graph$component <- NULL
    graph <- dodgr_components (graph)
    expect_identical (comp, graph$component)
    comp <- graph$component
    expect_identical (comp, graph$component)
})

test_that("contract graph", {
    graph <- weight_streetnet (hampi)
    graph_c <- dodgr_contract_graph (graph)$graph
    expect_true (nrow (graph_c) < nrow (graph))
})

test_that("compare heaps", {
    graph <- weight_streetnet (hampi)
    ch <- compare_heaps (graph, nverts = 100, replications = 1)
    expect_equal (nrow (ch), 12)
    # Test that all dodgr calculations are faster than igraph:
    igr <- which (grepl ("igraph", ch$test))
    #expect_true (ch$elapsed [igr] == max (ch$elapsed))
    # This actually fails on some machines (R oldrel on Windows) because elapsed
    # times are sometimes all very small *and equal*, so is turned off:
    if (test_all)
    {
        #expect_true (ch$elapsed [igr] > min (ch$elapsed))
    }
})

test_that("dodgr2sf", {
    hw <- weight_streetnet (hampi)
    y <- dodgr_to_sfc (hw)
    # y should have more linestrings than the original sf object:
    expect_true (length (y) > length (hw$geometry))
})

test_that("different geometry columns", {
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
    expect_error (net <- weight_streetnet (h2),
                    "Unable to determine geometry column")

    h2 <- data.frame (hampi) # remove sf class
    expect_error (net <- weight_streetnet (h2), "Unknown class")
})

test_that("no geom rownames", {
    hw0 <- weight_streetnet (hampi)
    g0 <- hampi$geometry
    attr (g0, "names") <- NULL # remove way IDs
    for (i in seq (g0))
        rownames (g0 [[i]]) <- NULL # remove all node IDs
    h2 <- hampi
    h2$geometry <- g0
    hw1 <- weight_streetnet (h2)
    expect_true (!identical (h2, hw1))
    expect_equal (ncol (hw0), 15) # has way_id column
    expect_equal (ncol (hw1), 14) # no way_id
    indx0 <- which (!names (hw0) %in% c ("from_id", "to_id", "way_id",
                                         "component", "time", "time_weighted"))
    indx1 <- which (!names (hw1) %in% c ("from_id", "to_id", "way_id",
                                         "component", "time", "time_weighted"))
    expect_identical (hw0 [, indx0], hw1 [, indx1])
    # components are not identical because ones of equal size are assigned
    # random numbers, but all other columns remain identical:
    # geom_num, edge_id, lon/lat values, d, d_weighted, and highway
})

test_that("keep cols", {
    hw0 <- weight_streetnet (hampi)
    expect_equal (ncol (hw0), 15)
    hw1 <- weight_streetnet (hampi, keep_cols = "foot")
    expect_equal (ncol (hw1), 16)
    expect_true ("foot" %in% names (hw1))
    expect_false ("foot" %in% names (hw0))
})

test_that ("weight_profiles", {
    graph0 <- weight_streetnet (hampi, wt_profile = "foot")
    graph1 <- weight_streetnet (hampi, wt_profile = 1)
    expect_equal (nrow (graph0), nrow (graph1))
    #expect_identical (graph0$d, graph1$d)
    #expect_true (!identical (graph0$d_weighted, graph1$d_weighted))
    wp <- dodgr::weighting_profiles$weighting_profiles
    wpf <- wp [wp$name == "foot", ]
    graph3 <- weight_streetnet (hampi, wt_profile = wpf)
    expect_identical (graph0$d_weighted, graph3$d_weighted)

    wpf$value [wpf$way == "path"] <- 0.9
    graph4 <- weight_streetnet (hampi, wt_profile = wpf)
    expect_true (!identical (graph0$d_weighted, graph4$d_weighted))

    names (wpf) [3] <- "Value"
    expect_error (graph4 <- weight_streetnet (hampi, wt_profile = wpf),
                  "Weighting profiles must have")

    g0 <- weight_streetnet (hampi)
    hampi2 <- hampi
    names (hampi2) [grep ("highway", names (hampi2))] <- "waytype"
    expect_error (weight_streetnet (hampi2),
                  "Please specify type_col to be used for weighting streetnet")
    g1 <- weight_streetnet (hampi2, type_col = "waytype")
    expect_identical (g0, g1)
    names (hampi2) [grep ("osm_id", names (hampi2))] <- "key"
    expect_message (weight_streetnet (hampi2, type_col = "waytype"),
                  "Using column width as ID column for edges")
    names (hampi2) [grep ("key", names (hampi2))] <- "id"
    expect_error (weight_streetnet (hampi2, type_col = "waytype"),
                  "Multiple potential ID columns")
    names (hampi2) [grep ("width", names (hampi2))] <- "w"
    expect_message (g1 <- weight_streetnet (hampi2, type_col = "waytype"),
                  "Using column id as ID column for edges")
    expect_identical (g0, g1)
})


test_that ("railway", {
    expect_error (g0 <- weight_streetnet (hampi, wt_profile = "rail"),
                  "Please use the weight_railway function for railway routing")
    expect_error (g0 <- weight_railway (hampi),
                  "Please specify type_col to be used for weighting railway")
    expect_message (g0 <- weight_railway (hampi, type_col = "highway"),
                    "Data has no columns named maxspeed")
    expect_silent (g0 <- weight_railway (hampi, type_col = "highway",
                                         keep_cols = NULL))
    expect_identical (g0$d, g0$d_weighted)
})

test_that ("points to graph", {
    bb <- attr (hampi$geometry, "bbox")
    n <- 100
    x <- bb [1] + (bb [3] - bb [1]) * runif (n)
    y <- bb [2] + (bb [4] - bb [2]) * runif (n)
    pts <- data.frame (x = x, y = y)
    net <- weight_streetnet (hampi)
    expect_message (index1 <- match_pts_to_graph (net, pts),
                    "First argument to match_pts_to_graph should be result of dodgr_vertices")

    v <- dodgr_vertices (net)
    expect_silent (index2 <- match_pts_to_graph (v, pts))
    expect_identical (index1, index2)

    colnames (pts) <- NULL
    expect_message (index3 <- match_pts_to_graph (v, pts),
                    "xy has no named columns; assuming order is x then y")
    expect_identical (index1, index3)

    pts <- data.frame (x = x, y = y, x2 = x)
    expect_error (index4 <- match_pts_to_graph (v, list (pts)),
                  "xy must be a matrix or data.frame")
    expect_error (index4 <- match_pts_to_graph (v, pts),
                  "xy must have only two columns")

    pts <- data.frame (x = x, y = y)
    expect_silent (index4 <- match_pts_to_graph (v, pts, connected = TRUE))
    expect_true (!identical (index1, index4))

    class (pts) <- c (class (pts), "tbl")
    expect_silent (index5 <- match_pts_to_graph (v, pts, connected = TRUE))
    expect_identical (index4, index5)

    pts <- sf::st_as_sf (pts, coords = c (1, 2), crs = 4326)
    expect_silent (index6 <- match_pts_to_graph (v, pts, connected = TRUE))
    expect_identical (index4, index6)
    expect_silent (index7 <- match_points_to_graph (v, pts, connected = TRUE))
    expect_identical (index4, index7)
})

test_that ("graph columns", {
    graph <- data.frame (weight_streetnet (hampi)) # rm dodgr_streetnet class
    nf <- 100
    nt <- 50
    from <- sample (graph$from_id, size = nf)
    to <- sample (graph$from_id, size = nt)
    d <- dodgr_dists (graph, from = from, to = to)
    graph$from_id2 <- graph$from_id
    expect_error (d <- dodgr_dists (graph, from = from, to = to),
                  "Unable to determine column with ID of from vertices")
    graph$from_id2 <- NULL
    graph$to_id2 <- graph$to_id
    expect_error (d <- dodgr_dists (graph, from = from, to = to),
                  "Unable to determine column with ID of to vertices")

    graph <- data.frame (weight_streetnet (hampi)) # rm dodgr_streetnet class
    graph$from_lat <- NULL
    expect_error (d <- dodgr_dists (graph, from = from, to = to),
                  "Unable to determine coordinate columns of graph")

    graph <- data.frame (weight_streetnet (hampi)) # rm dodgr_streetnet class
    graph$d2 <- graph$d
    expect_error (d <- dodgr_dists (graph, from = from, to = to),
                  "Unable to determine distance column in graph")
    graph$d2 <- NULL
    graph$d_wt <- graph$d_weighted
    expect_error (d <- dodgr_dists (graph, from = from, to = to),
                  "Unable to determine weight column in graph")
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
    expect_error (get_pts_index (graph, gr_cols, vert_map, pts),
                  "Unable to determine geographical coordinates in from/to")
    rownames (pts) <- seq (n)
    expect_error (get_pts_index (graph, gr_cols, vert_map, pts),
                  "Unable to determine geographical coordinates in from/to")

    pts <- cbind (runif (n), runif (n), seq (n))
    colnames (pts) <- c ("a", "b", "id")
    expect_error (get_pts_index (graph, gr_cols, vert_map, pts),
                  "Unable to determine geographical coordinates in from/to")
    colnames (pts) <- c ("x", "y", "id")
    expect_is (get_pts_index (graph, gr_cols, vert_map, pts), "numeric")
    expect_length (get_pts_index (graph, gr_cols, vert_map, pts), n)

    pts <- data.frame (pts)
    expect_is (get_pts_index (graph, gr_cols, vert_map, pts), "numeric")
    expect_length (get_pts_index (graph, gr_cols, vert_map, pts), n)
})
