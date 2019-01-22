context("fundamental cycles")

test_all <- (identical (Sys.getenv ("MPADGE_LOCAL"), "true") |
             identical (Sys.getenv ("TRAVIS"), "true"))

test_that("dodgr_fundamental_cycles", {
              net <- weight_streetnet (hampi)
              graph <- dodgr_contract_graph (net)$graph
              expect_error (x <- dodgr_fundamental_cycles (),
                            "graph must be provided")
              expect_error (x <- dodgr_fundamental_cycles (graph = "a"),
                            "graph must be a data.frame object")
              expect_silent (x <- dodgr_fundamental_cycles (graph))
              expect_is (x, "list")
              expect_length (x, 43)
             })

test_that("cycles_with_max_graph_size", {
              net <- weight_streetnet (hampi)
              expect_message (
                    x <- dodgr_fundamental_cycles (graph = net,
                                                   graph_max_size = 1000),
                              "Now computing fundamental cycles")
              expect_is (x, "list")
              expect_length (x, 66) # more cycles than before!

              expect_silent (
                    xf <- dodgr_full_cycles (graph = net, graph_max_size = 1000))
              # full_cycles creates the contracted graph, which is < 1000!
              expect_length (xf, 43)
             })

test_that("sflines_to_poly", {
              expect_silent (p <- dodgr_sflines_to_poly (hampi))
              expect_is (hampi$geometry, "sfc_LINESTRING")
              expect_is (p, "sfc_POLYGON")
              expect_equal (length (p), 59)

              net <- weight_streetnet (hampi, wt_profile = 1)

              net1 <- net [net$component == 1, ]
              net1$edge_id <- seq (nrow (net1))
              p1 <- dodgr_full_cycles (net1)

              net2 <- net [net$component == 2, ]
              net2$edge_id <- seq (nrow (net2))
              p2 <- dodgr_full_cycles (net2)

              # `dodgr_sflines_to_poly` analyses components seperately, so p1 +
              # p2 should give same result:
              expect_equal (length (p1) + length (p2), length (p))
             })
