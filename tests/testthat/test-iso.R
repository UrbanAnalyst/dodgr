context("iso")

test_all <- (identical (Sys.getenv ("MPADGE_LOCAL"), "true") |
             identical (Sys.getenv ("TRAVIS"), "true"))

source ("../sc-conversion-fns.R")

test_that("isodists", {
              expect_silent (hsc <- sf_to_sc (hampi))
              # This all exists just to test the next line:
              requireNamespace ("geodist")
              requireNamespace ("dplyr")
              expect_silent (net <- weight_streetnet (hsc,
                                                      wt_profile = "bicycle"))
              npts <- 100
              from <- sample (net$.vx0, size = npts)
              dlim <- c (1, 2, 5, 10, 20) * 100
              expect_silent (d <- dodgr_isodists (net, from = from, dlim))
              expect_is (d, "data.frame")
              expect_equal (ncol (d), 5)
              expect_identical (names (d), c ("from", "dlim", "id", "x", "y"))
              expect_true (nrow (d) > 0)
              # some points may give no iso countours, so the following 2 may
              # not always be equal:
              expect_true ((npts - length (unique (d$from))) >= 0)
              expect_true (length (unique (d$dlim)) <= length (dlim))
})

test_that ("errors", {
              expect_silent (hsc <- sf_to_sc (hampi))
              expect_silent (net <- weight_streetnet (hsc,
                                                      wt_profile = "bicycle"))
              npts <- 100
              from <- sample (net$.vx0, size = npts)
              dlim <- c (1, 2, 5, 10, 20) * 100
              expect_error (d <- dodgr_isodists (net, from = from),
                            "dlim must be specified")
              expect_error (d <- dodgr_isodists (net, from = from, dlim = "a"),
                            "dlim must be numeric")
             })

test_that("isochrones", {
              expect_silent (hsc <- sf_to_sc (hampi))
              # This all exists just to test the next line:
              requireNamespace ("geodist")
              requireNamespace ("dplyr")
              expect_silent (net <- weight_streetnet (hsc,
                                                      wt_profile = "bicycle"))
              npts <- 100
              from <- sample (net$.vx0, size = npts)
              tlim <- c (5, 10, 20, 30, 60) * 60 # times in seconds
              expect_silent (x <- dodgr_isochrones (net, from = from, tlim))
              expect_is (x, "data.frame")
              expect_equal (ncol (x), 5)
              expect_identical (names (x), c ("from", "tlim", "id", "x", "y"))
              expect_true (nrow (x) > 0)
              # some points may give no iso countours, so the following 2 may
              # not always be equal:
              expect_true ((npts - length (unique (x$from))) >= 0)
              expect_true (length (unique (x$tlim)) <= length (tlim))
})
