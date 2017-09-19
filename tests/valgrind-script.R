library (dodgr)
graph <- weight_streetnet (hampi)
graph <- dodgr_sample (graph, nverts = 100)
d <- dodgr_dists (graph)
