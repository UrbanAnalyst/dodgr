context("dodgr streetnet")

test_that("streetnet bbox", {
              n <- 24
              bbox <- matrix (runif (n), ncol = 2)
              bb <- process_bbox (bbox, NULL, 0)
              expect_is (bb, "list")
              expect_length (bb, 2)
              expect_equal (nrow (bb$bbox), 2)
              expect_equal (nrow (bb$bbox_poly), n / 2)

              expect_message (bb2 <- process_bbox (list (bbox), NULL, 0),
                              "selecting the first polygon from bbox")
              #expect_identical (bb, bb2)
              # TODO: They should be identical - fix!
})

