context("weighting profiles")

test_that("wp", {
              f <- file.path (tempdir (), "wp")
              expect_false (file.exists (paste0 (f, ".json")))
              write_dodgr_wt_profile (f)
              expect_true (file.exists (paste0 (f, ".json")))
              w <- read_dodgr_wt_profile (f)
              expect_identical (w, dodgr::weighting_profiles)
})
