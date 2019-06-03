context("weighting profiles")

test_that("wp", {
              f <- file.path (tempdir (), "wp")
              expect_false (file.exists (paste0 (f, ".json")))
              expect_silent (write_dodgr_wt_profile (f))
              expect_true (file.exists (paste0 (f, ".json")))
              w <- read_dodgr_wt_profile (f)
              expect_identical (w, dodgr::weighting_profiles)
})

test_that ("local wt_profile", {
        expect_error (write_dodgr_wt_profile (),
                      "file name must be given")

        f <- file.path (tempdir (), "wp")
        expect_silent (write_dodgr_wt_profile (f))
        n0 <- weight_streetnet (hampi, wt_profile = "foot")
        n1 <- weight_streetnet (hampi, wt_profile = "foot",
                              wt_profile_file = f)
        attr (n0, "px") <- NULL
        attr (n1, "px") <- NULL
        expect_identical (n0, n1)

        w <- dodgr::weighting_profiles
        w$weighting_profiles$max_speed [w$weighting_profiles$name == "foot" &
                                        w$weighting_profiles$max_speed == 5] <- 8

        f <- paste0 (tools::file_path_sans_ext (f), ".json")
        con <- file (f, open = "wt")
        wpj <- jsonlite::toJSON (w, pretty = TRUE)
        writeLines (wpj, con)
        close (con)

        n2 <- weight_streetnet (hampi, wt_profile = "foot",
                                wt_profile_file = f)
        expect_true (mean (n2$time) < mean (n0$time))
        expect_true (mean (n2$time_weighted) < mean (n0$time_weighted))
        expect_identical (n2$d, n0$d)
        expect_identical (n2$d_weighted, n0$d_weighted)

        expect_error (weight_streetnet (hampi, wt_profile = 1:2),
                      "wt_profile can only be one element")
})

test_that ("weight_streetnet wt_profile_file", {
        expect_silent (graph <- weight_streetnet (hampi, wt_profile = 1))
        expect_identical (graph$d, graph$d_weighted)
})
