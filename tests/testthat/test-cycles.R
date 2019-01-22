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
             })
