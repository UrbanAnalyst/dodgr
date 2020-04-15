context("iso")

test_all <- (identical (Sys.getenv ("MPADGE_LOCAL"), "true") |
             identical (Sys.getenv ("TRAVIS"), "true"))

if (!test_all)
    RcppParallel::setThreadOptions (numThreads = 2)

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

skip_on_cran ()
skip_if (!test_all)

test_that ("turn penalty", {
              hsc <- sf_to_sc (hampi)
              net0 <- weight_streetnet (hsc, wt_profile = "bicycle")
              npts <- 100
              from <- sample (net0$.vx0, size = npts)
              dlim <- c (1, 2, 5, 10, 20) * 100
              d0 <- dodgr_isodists (net0, from = from, dlim)

              net <- weight_streetnet (hsc, wt_profile = "bicycle",
                                       turn_penalty = TRUE)
              expect_true (nrow (net) > nrow (net0))
              d <- dodgr_isodists (net, from = from, dlim)
              # d includes compound vertices with "_start" suffix, and routes
              # differently because of turning angles
              expect_false (identical (d0, d))
              #expect_true (length (grep ("_start", d$from)) > 0)
              #expect_false (length (grep ("_start", d0$from)) > 0)
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
              expect_error (d <- dodgr_isoverts (net, from = from,
                                                 dlim = dlim, tlim = 500),
                            "Only one of dlim or tlim can be provided")

              net <- weight_streetnet (hampi)
              expect_error (d <- dodgr_isochrones (net, from = from, tlim = 500),
                            "isochrones can only be calculated from SC-class networks")
              expect_error (d <- dodgr_isoverts (net, from = from, tlim = 500),
                            "isoverts can only be calculated from SC-class networks")
             })

test_that("isochrones", {
              expect_silent (hsc <- sf_to_sc (hampi))
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

test_that("isoverts", {
              expect_silent (hsc <- sf_to_sc (hampi))
              expect_silent (net <- weight_streetnet (hsc,
                                                      wt_profile = "bicycle"))
              npts <- 100
              from <- sample (net$.vx0, size = npts)
              dlim <- c (1, 2, 5, 10, 20) * 100
              expect_silent (dd <- dodgr_isodists (net, from = from, dlim))
              expect_silent (d <- dodgr_isoverts (net, from = from, dlim))
              expect_identical (names (dd), names (d))
              # isodists may not return all from pts, so the following may not
              # always hold:
              #expect_identical (unique (d$from), unique (dd$from))
              # dd has all vertices within isodistance hulls; d has only those
              # on the actual hulls, so far fewer vertices
              expect_true (nrow (d) > nrow (dd))

              tlim <- c (60, 120, 300)
              expect_silent (v <- dodgr_isoverts (net, from = from, tlim = tlim))
              expect_true ("tlim" %in% names (v))
              expect_true ("dlim" %in% names (d))
              expect_true (all (paste0 (tlim) %in% unique (v$tlim)))
})
