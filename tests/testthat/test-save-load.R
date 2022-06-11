context ("save and load")

testthat::skip_on_cran ()

test_all <- (identical (Sys.getenv ("MPADGE_LOCAL"), "true") |
    identical (Sys.getenv ("GITHUB_WORKFLOW"), "test-coverage"))

if (!test_all) {
    RcppParallel::setThreadOptions (numThreads = 2)
}

test_that ("save and load", {
    net0 <- weight_streetnet (hampi)
    f <- file.path (tempdir (), "junk")
    expect_silent (
        dodgr_save_streetnet (net0, f)
    )
    expect_false (file.exists (f))
    f <- paste0 (f, ".Rds")
    expect_true (file.exists (f))

    x <- readRDS (f)
    expect_is (x, "list")

    expect_identical (
        names (x),
        c (
            "graph", "verts", "graph_c",
            "verts_c", "edge_map", "junctions"
        )
    )

    flist0 <- list.files (tempdir (), pattern = "^dodgr\\_")
    clear_dodgr_cache ()

    net1 <- dodgr_load_streetnet (f)
    expect_equal (net0, net1)
    flist1 <- list.files (tempdir (), pattern = "^dodgr\\_")

    expect_true (all (flist1 %in% flist0))
})
