devtools::load_all ("../dodgr")
#streetnet <- dodgr_streetnet ("hampi india")
#nms <- c ("osm_id", "highway", "oneway", "geometry")
#streetnet <- streetnet [names (streetnet) %in% nms, ]
#save (streetnet, file = "hampi.rda")
load ("tests/hampi.rda")
graph <- weight_streetnet (streetnet)
graph <- dodgr_sample (graph, nverts = 100)
d <- dodgr_dists (graph)
