context("SC")

#library (osmdata)
#devtools::load_all ("../../ropensci/osmdata", export_all = FALSE)
#h2 <- opq ("hampi india") %>%
#    add_osm_feature (key = "highway") %>%
#    osmdata_sc ()

source ("../sc-conversion-fns.R")

test_that("SC", {
              expect_silent (hsc <- sf_to_sc (hampi))
              # This all exists just to test the next line:
              requireNamespace ("geodist")
              requireNamespace ("dplyr")
              expect_silent (net_sc <- weight_streetnet (hsc))
              expect_is (net_sc, "data.frame")
              expect_true (nrow (net_sc) > 0)

              net_sf <- weight_streetnet (hampi)
              # The SC version is sanitised by, for example, removing duplicate
              # edges, so is a more compact graph:
              expect_true (nrow (net_sf) > nrow (net_sc))
              v_sc <- dodgr_vertices (net_sc)
              v_sf <- dodgr_vertices (net_sf)
              expect_true (nrow (v_sf) > nrow (v_sc))

              class (hsc) <- class (hsc) [!class (hsc) %in% "osmdata_sc"]
              expect_error (net_sc <- weight_streetnet (hsc),
                            paste0 ("weight_streetnet currently only works for ",
                                    "'sc'-class objects extracted with"))
})

test_that ("elevation", {
               expect_silent (hsc <- sf_to_sc (hampi))
               expect_silent (net_sc <- weight_streetnet (hsc))
               hsc$vertex$z_ <- runif (nrow (hsc$vertex)) * 10
               expect_silent (net_sc2 <- weight_streetnet (hsc))
               expect_true (ncol (net_sc2) == (ncol (net_sc) + 1))
})

test_that("dodgr_times", {
              # dists and times should be strongly correlated:
              expect_silent (hsc <- sf_to_sc (hampi))
              expect_silent (net_sc <- weight_streetnet (hsc))
              v <- dodgr_vertices (net_sc)
              set.seed (1)
              from <- sample (v$id, 100)
              to <- sample (v$id, 100)
              d <- dodgr_dists (net_sc, from = from, to = to)
              t1 <- dodgr_times (net_sc, from = from, to = to)
              r2 <- cor (as.numeric (d), as.numeric (t1),
                         use = "pairwise.complete.obs")
              expect_true (r2 < 1)
              # with no turn angles, the should be just scaled versions

              # calculate times with turning angles, such that resultant network
              # includes compound junction edges
              expect_silent (net_sc2 <- weight_streetnet (hsc, turn_angle = TRUE))
              expect_true (nrow (net_sc2) > nrow (net_sc))
              v0 <- net_sc2$.vx0 [grep ("_start", net_sc2$.vx0)]
              v0 <- gsub ("_start", "", v0)
              v1 <- net_sc2$.vx1 [grep ("_end", net_sc2$.vx1)]
              v1 <- gsub ("_end", "", v1)
              from [from %in% v0] <- paste0 (from [from %in% v0], "_start")
              to [to %in% v1] <- paste0 (to [to %in% v1], "_end")
              t2 <- dodgr_times (net_sc2, from = from, to = to)
              r2 <- cor (as.numeric (t1), as.numeric (t2),
                         use = "pairwise.complete.obs")
              expect_true (r2 < 1)
              expect_true (r2 > 0.95)
              # These times should be longer:
              expect_true (mean (t2 - t1, na.rm = TRUE) > 0)

              # times with contracted graph should be identical:
              net_sc2_c <- dodgr_contract_graph (net_sc2)
              v <- dodgr_vertices (net_sc2_c)
              set.seed (1)
              from <- sample (v$id, 100)
              to <- sample (v$id, 100)

              #t1 <- dodgr_times (net_sc2, from = from, to = to)
              #t2 <- dodgr_times (net_sc2_c, from = from, to = to)
              #dtime <- max (abs (t1 - t2), na.rm = TRUE)
              #expect_true (dtime < 1e-6)
})
