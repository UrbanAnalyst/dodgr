context("cache")

test_all <- (identical (Sys.getenv ("MPADGE_LOCAL"), "true") |
             identical (Sys.getenv ("TRAVIS"), "true"))

#library (osmdata)
#devtools::load_all ("../../ropensci/osmdata", export_all = FALSE)
#h2 <- opq ("hampi india") %>%
#    add_osm_feature (key = "highway") %>%
#    osmdata_sc ()

source ("../sc-conversion-fns.R")

test_that("cache on", {
              expect_silent (hsc <- sf_to_sc (hampi))
              requireNamespace ("geodist")
              requireNamespace ("dplyr")
              expect_silent (graph <- weight_streetnet (hsc))
              expect_message (graph <- dodgr_components (graph),
                              "graph already has a component column")
              expect_silent (v0 <- dodgr_vertices (graph))
              expect_silent (graph_c <- dodgr_contract_graph (graph))
              expect_silent (v <- dodgr_vertices (graph_c))

              n <- 100
              pts <- sample (v$id, size = n)
              pts <- pts [which (pts %in% graph_c$.vx0 & pts %in% graph_c$.vx1)]
              fmat <- array (1, dim = c (n, n))

              # aggregate flows from graph without turning angles:
              expect_silent (graphf <- dodgr_flows_aggregate (graph_c,
                                                              from = pts,
                                                              to = pts,
                                                              flow = fmat))
              expect_silent (graphf <- dodgr_uncontract_graph (graphf))
              expect_silent (graphf <- merge_directed_flows (graphf))

              # then turn angle graph
              grapht <- weight_streetnet (hsc, wt_profile = "bicycle",
                                          turn_penalty = TRUE, left_side = TRUE)
              # grapht has extra compound edges for turning angles:
              expect_true (nrow (grapht) > nrow (graph))
              grapht_c <- dodgr_contract_graph (grapht)
              expect_true (nrow (grapht_c) > nrow (graph_c))
              expect_silent (graphtf <- dodgr_flows_aggregate (grapht_c,
                                                              from = pts,
                                                              to = pts,
                                                              flow = fmat))
              expect_silent (graphtf <- dodgr_uncontract_graph (graphtf))
              expect_silent (graphtf <- merge_directed_flows (graphtf))

              expect_silent (graphtf <- dodgr_flows_disperse (grapht_c,
                                                              from = pts,
                                                              dens = rep (1, n)))

})

test_that("cache off", {
              expect_silent (clear_dodgr_cache ())
              expect_silent (dodgr_cache_off ())
              expect_silent (hsc <- sf_to_sc (hampi))
              expect_silent (graph <- weight_streetnet (hsc))
              expect_message (graph <- dodgr_components (graph),
                              "graph already has a component column")
              expect_silent (v0 <- dodgr_vertices (graph))
              expect_silent (graph_c <- dodgr_contract_graph (graph))
              expect_silent (v <- dodgr_vertices (graph_c))

              n <- 100
              pts <- sample (v$id, size = n)
              pts <- pts [which (pts %in% graph_c$.vx0 & pts %in% graph_c$.vx1)]
              fmat <- array (1, dim = c (n, n))

              # aggregate flows from graph without turning angles:
              expect_silent (graphf <- dodgr_flows_aggregate (graph_c,
                                                              from = pts,
                                                              to = pts,
                                                              flow = fmat))
              expect_silent (graphf <- dodgr_uncontract_graph (graphf))
              expect_silent (graphf <- merge_directed_flows (graphf))

              # then turn angle graph
              expect_silent (grapht <- weight_streetnet (hsc,
                                                         wt_profile = "bicycle",
                                                         turn_penalty = TRUE,
                                                         left_side = TRUE))
              expect_silent (grapht_c <- dodgr_contract_graph (grapht))
              expect_silent (graphtf <- dodgr_flows_aggregate (grapht_c,
                                                              from = pts,
                                                              to = pts,
                                                              flow = fmat))
              expect_silent (graphtf <- dodgr_uncontract_graph (graphtf))
              expect_silent (graphtf <- merge_directed_flows (graphtf))

              expect_silent (graphtf <- dodgr_flows_disperse (grapht_c,
                                                              from = pts,
                                                              dens = rep (1, n)))

              expect_silent (dodgr_cache_on ())

})
