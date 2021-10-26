context("save and load")

testthat::skip_on_cran ()

test_all <- (identical (Sys.getenv ("MPADGE_LOCAL"), "true") |
             identical (Sys.getenv ("GITHUB_WORKFLOW"), "test-coverage"))

if (!test_all)
    RcppParallel::setThreadOptions (numThreads = 2)

test_that("save", {
    net <- weight_streetnet (hampi)
    f <- file.path (tempdir (), "junk")
    expect_silent (
        dodgr_save_streetnet (net, f)
        )
    expect_false (file.exists (f))
    f <- paste0 (f, ".Rds")
    expect_true (file.exists (f))

    x <- readRDS (f)
    expect_is (x, "list")

    expect_identical (names (x),
                      c ("graph", "verts", "graph_c",
                         "verts_c", "edge_map", "junctions"))
})
