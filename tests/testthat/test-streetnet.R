context("dodgr streetnet")

test_that("streetnet bbox", {
              n <- 12
              bbox <- cbind (runif (n), 2 * runif (n))
              bb <- process_bbox (bbox, NULL, 0)
              expect_is (bb, "list")
              expect_length (bb, 2)
              expect_equal (nrow (bb$bbox), 2)
              expect_equal (nrow (bb$bbox_poly), n)

              bbox2 <- as.vector (apply (bbox, 2, range))
              bb2 <- process_bbox (bbox2, NULL, 0)
              expect_identical (bb$bbox, bb2$bbox)

              expect_message (bb2 <- process_bbox (list (bbox), NULL, 0),
                              "selecting the first polygon from bbox")
              expect_identical (bb, bb2)

              bbox <- list (matrix (letters [ceiling (runif (n) * 26)],
                                    ncol = 2))
              expect_error (bb <- process_bbox (bbox, NULL, 0),
                            "bbox is a list, so items must be numeric")
              bbox <- runif (6)
              expect_error (bb <- process_bbox (bbox, NULL, 0),
                            "bbox must have four numeric values")
})

test_that ("streetnet pts", {
              n <- 12
              pts <- cbind (runif (n), 2 * runif (n))
              expect_error (bb <- process_bbox (pts = pts, expand = 0),
                            paste0 ("Can not unambiguously determine ",
                                    "coordinates in graph"))

              colnames (pts) <- c ("x", "y")
              expect_silent (bb <- process_bbox (pts = pts, expand = 0))
              expect_silent (bb2 <- process_bbox (bbox = pts, expand = 0))

              bb2_bb <- c (bb2$bbox [1, 1], bb2$bbox [1, 2],
                           bb2$bbox [2, 1], bb2$bbox [2, 2])
              names (bb2_bb) <- NULL
              # -> as.vector (t (bb2$bbox))
              expect_identical (bb2_bb, bb$bbox)
})

test_that ("streetnet times", {
               expect_error (graph <- weight_streetnet (hampi, turn_angle = TRUE),
                             paste0 ("Turn-angle calculations only currently ",
                                     "implemented for street network data ",
                                     "generated with"))
               expect_silent (graph <- weight_streetnet (hampi))
               h <- hampi
               names (h) [names (h) == "osm_id"] <- "id"
               expect_silent (graph2 <- weight_streetnet (h, id_col = "id"))
               expect_identical (graph, graph2)

               h$id <- NULL
               expect_message (graph3 <- weight_streetnet (h),
                               "Using column width as ID column for edges")

               h <- hampi
               names (h$geometry) <- NULL
               graph4 <- weight_streetnet (h)
               expect_identical (graph$edge_id, seq (nrow (graph)))

               h$oneway_bicycle <- h$oneway
               graph5 <- weight_streetnet (h)
               expect_identical (graph5, graph4)

               expect_error (weight_streetnet (hampi,
                                               wt_profile = list (1)),
                             "Custom named profiles must be vectors")
})
