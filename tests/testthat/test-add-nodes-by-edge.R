library(dplyr)
library(glue)
library(testthat)

# test_all <- (identical (Sys.getenv ("MPADGE_LOCAL"), "true") ||
#                identical (Sys.getenv ("GITHUB_JOB"), "test-coverage"))
# 
# skip_if (!test_all)

library(dplyr)
dodgr_cache_off ()
clear_dodgr_cache ()

test_that ("add_nodes_to_graph_by_edge single point per edge", {
  
  # Load a sample graph
  graph <- weight_streetnet (hampi, wt_profile = "foot")%>%
    mutate(graph_index=1:n())
  verts <- dodgr_vertices (graph)
  
  # Create a small set of points that will likely match to different edges
  set.seed (1)
  npts <- 5
  xy <- data.frame (
    x = min (verts$x) + runif (npts) * diff (range (verts$x)),
    y = min (verts$y) + runif (npts) * diff (range (verts$y))
  )
  # Match points to graph to verify they match to different edges
  pts <- match_pts_to_graph (graph, xy, distances = TRUE)%>%
    mutate(xy_index=1:n())%>%
    group_by(index)%>%
    dplyr::filter(n()==1)

  # Filter to keep only points that match to unique edges
  xy_single <- xy [pts$xy_index, ]
  xy_single <- xy_single[1,]
  # find bidirectional edges
  bi <- graph[pts$index,]%>%
    select(from_id=to_id, to_id=from_id, graph_index)%>%
    inner_join(graph%>%rename(graph_index_bi=graph_index))
  bi_index <- sort(c(bi$graph_index, bi$graph_index_bi))
  
  # Process with both functions
  graph1 <- add_nodes_to_graph (graph, xy_single, intersections_only = TRUE, dist_tol = 0)
  graph2 <- add_nodes_to_graph_by_edge (graph, xy_single, intersections_only = TRUE, dist_tol = 0)
  
  # Compare results
  expect_equal (nrow (graph1), nrow (graph2))
  
  # Compare total distance in the graph
  total_dist1 <- sum (graph1$d)
  total_dist2 <- sum (graph2$d)
  expect_equal (total_dist1, total_dist2, tolerance = 1e-6)
  
  # Compare structure
  expect_equal (ncol (graph1), ncol (graph2))
  
  # Compare number of unique vertices
  verts1 <- unique (c (graph1$from, graph1$to))
  verts2 <- unique (c (graph2$from, graph2$to))
  expect_equal (length (verts1), length (verts2))
  
  
  # Compare with edge to point creation
  
  # Process with both functions
  graph1 <- add_nodes_to_graph (graph, xy_single, intersections_only = FALSE, dist_tol = 0)
  graph2 <- add_nodes_to_graph_by_edge (graph, xy_single, intersections_only = FALSE, dist_tol = 0)
  expect_equal (sum(graph1$d), sum(graph2$d), tolerance = 1e-6)
  expect_equal (sum(graph1$time), sum(graph2$time), tolerance = 1e-6)
  expect_equal (sum(graph1$time_weighted), sum(graph2$time_weighted), tolerance = 1e-6)
  expect_equal (sum(graph1$d_weighted), sum(graph2$d_weighted), tolerance = 1e-6)
  
  
  
  graph2 <- add_nodes_to_graph_by_edge (graph, xy_single, intersections_only = FALSE, dist_tol = 0, wt_profile = "foot", highway="residential")
  expect_equal (sum(graph1$d), sum(graph2$d), tolerance = 1e-6)
  expect_equal (sum(graph1$time), sum(graph2$time), tolerance = 1e-6)
  expect_equal (sum(graph1$time_weighted), sum(graph2$time_weighted), tolerance = 1e-6)
  expect_equal (sum(graph1$d_weighted), sum(graph2$d_weighted), tolerance = 1e-6)
  
  
  
  graph2 <- add_nodes_to_graph_by_edge (graph, xy_single, intersections_only = FALSE, dist_tol = 0, wt_profile = "motorcar", highway="residential")
  expect_equal (sum(graph1$d), sum(graph2$d), tolerance = 1e-6)
  expect_gt (sum(graph1$time), sum(graph2$time))
  expect_gt (sum(graph1$time_weighted), sum(graph2$time_weighted))
})

test_that ("add_nodes_to_graph_by_edge multiple points per edge", {
  
  # Load a sample graph
  graph <- weight_streetnet (hampi, wt_profile = "foot")
  verts <- dodgr_vertices (graph)
  
  # Create a set of points where multiple points will match to the same edge
  set.seed (2)
  npts <- 20 # More points increases chance of multiple matches
  xy <- data.frame (
    x = min (verts$x) + runif (npts) * diff (range (verts$x)),
    y = min (verts$y) + runif (npts) * diff (range (verts$y))
  )

  # Match points to graph
  pts <- match_pts_to_graph (graph, xy, distances = TRUE)
  edge_counts <- table (pts$index)
  
  # Identify edges with multiple points
  multi_point_edges <- as.integer (names (edge_counts [edge_counts > 1]))
  
  # Verify we have at least one edge with multiple points
  expect_true (length (multi_point_edges) > 0)
  
  # Create a dataset with multiple points per edge
  multi_point_indices <- which (pts$index %in% multi_point_edges)
  xy_multi <- xy [multi_point_indices, ]
  # Process with both functions
  graph1 <- add_nodes_to_graph (graph, xy_multi)
  graph2 <- add_nodes_to_graph_by_edge (graph, xy_multi)
  
  # Compare results
  
  # The edge-based approach should be more efficient with multiple points
  # So the total distance should be less or equal
  total_dist1 <- sum (graph1$d)
  total_dist2 <- sum (graph2$d)
  
  # The edge-based approach should create fewer edges
  expect_true (nrow (graph2) <= nrow (graph1))
  
  # The total distance should be less or equal for the edge-based approach
  expect_true (total_dist2 <= total_dist1 )
  
  # Print the difference for diagnostic purposes
  cat ("\nMultiple points per edge comparison:\n")
  cat ("Original graph edges:", nrow (graph), "\n")
  cat ("add_nodes_to_graph edges:", nrow (graph1), "\n")
  cat ("add_nodes_to_graph_by_edge edges:", nrow (graph2), "\n")
  cat ("add_nodes_to_graph total distance:", total_dist1, "\n")
  cat ("add_nodes_to_graph_by_edge total distance:", total_dist2, "\n")
  cat ("Distance ratio (edge/point):", total_dist2 / total_dist1, "\n")
})

test_that ("add_nodes_to_graph_by_edge with mixed point distribution", {
  
  # Load a sample graph
  graph <- weight_streetnet (hampi, wt_profile = "foot")%>%
    mutate(edge_id=as.character(edge_id))%>%
    std_graph()
  verts <- dodgr_vertices (graph)
  
  # Create a larger set of points with mixed distribution
  set.seed (3)
  npts <- 3
  xy <- data.frame (
    x = min (verts$x) + runif (npts) * diff (range (verts$x)),
    y = min (verts$y) + runif (npts) * diff (range (verts$y))
  )
  #xy <- rbind(xy,xy)
  # Process with both functions
  graph1 <- add_nodes_to_graph (graph, xy, dist_tol = 0)
  graph2 <- add_nodes_to_graph_by_edge (graph, xy, dist_tol = 0)
  
  g1 <- graph1%>%
    anti_join(graph)%>%
    filter(grepl("_", edge_id))%>%
    tidyr::separate(edge_id, c("edge_id", "edge_id_seq"), sep="_")%>%
    left_join(graph%>%select(edge_id, highway), by="edge_id")%>%
    filter(highway.x!=highway.y)
  
  g1 <- graph1%>%anti_join(graph)%>%std_graph()
  g2 <- graph2%>%anti_join(graph)%>%std_graph()
  g1%>%anti_join(g2, by=c("from_lat", "from_lon", "to_lat", "to_lon", "d", "d_weighted", "highway", "time", "time_weighted"))%>%left_join(g2, by=c("from_lat", "from_lon", "to_lat", "to_lon"))%>%select(sort(names(.)))%>%View()
  
  # Compare results
  
  # The edge-based approach should generally be more efficient
  total_dist1 <- sum (graph1$d)
  total_dist2 <- sum (graph2$d)
  cat ("Distance ratio (edge/point):", total_dist2 / total_dist1, "\n")
  
  # Print the difference for diagnostic purposes
  cat ("\nMixed point distribution comparison:\n")
  cat ("Original graph edges:", nrow (graph), "\n")
  cat ("add_nodes_to_graph edges:", nrow (graph1), "\n")
  cat ("add_nodes_to_graph_by_edge edges:", nrow (graph2), "\n")
  cat ("add_nodes_to_graph total distance:", total_dist1, "\n")
  cat ("add_nodes_to_graph_by_edge total distance:", total_dist2, "\n")
  cat ("Distance ratio (edge/point):", total_dist2 / total_dist1, "\n")
  
  # Calculate efficiency metrics
  edge_increase1 <- nrow (graph1) - nrow (graph)
  edge_increase2 <- nrow (graph2) - nrow (graph)
  
  cat ("Edge increase (point method):", edge_increase1, "\n")
  cat ("Edge increase (edge method):", edge_increase2, "\n")
  cat ("Edge efficiency ratio:", edge_increase2 / edge_increase1, "\n")
  
  # The edge-based approach should be more efficient in terms of edges added
  expect_true (edge_increase2 <= edge_increase1)
  
  # Verify that both graphs have the same number of unique vertices
  # (excluding the intermediate vertices created during edge splitting)
  verts1 <- unique (c (graph1$from_id, graph1$to_id))
  verts2 <- unique (c (graph2$from_id, graph2$to_id))
  
  # The number of vertices might differ slightly due to different splitting approaches
  # but the difference should be small relative to the total
  vertex_diff_ratio <- abs (length (verts1) - length (verts2)) / length (verts1)
  expect_true (vertex_diff_ratio < 0.1) # Allow up to 10% difference
})

test_that ("add_nodes_to_graph_by_edge preserves edge properties", {
  
  # Load a sample graph
  graph <- weight_streetnet (hampi, wt_profile = "foot")
  
  # Create a small set of test points
  set.seed (4)
  npts <- 10
  verts <- dodgr_vertices (graph)
  xy <- data.frame (
    x = min (verts$x) + runif (npts) * diff (range (verts$x)),
    y = min (verts$y) + runif (npts) * diff (range (verts$y))
  )
  
  # First, identify which edges will be split by finding the matching edges for each point
  pts <- match_pts_to_graph (graph, xy, distances = TRUE)
  
  # Get the original edges that will be split
  edges_to_split <- graph[pts$index, ]
  
  # Calculate the ratios for these specific edges
  orig_d_ratios <- edges_to_split$d_weighted / edges_to_split$d
  orig_time_ratios <- edges_to_split$time_weighted / edges_to_split$time
  
  # Process with both functions
  graph1 <- add_nodes_to_graph (graph, xy)
  graph2 <- add_nodes_to_graph_by_edge (graph, xy)
  
  # Check that all required columns are preserved
  expect_equal (sort (names (graph)), sort (names (graph1)))
  expect_equal (sort (names (graph)), sort (names (graph2)))
  
  # For each original edge that was split, find the corresponding new edges
  for (i in seq_len(nrow(edges_to_split))) {
    edge_id <- edges_to_split$edge_id[i]
    
    # Find new edges in graph1 that replaced this edge
    # These will have edge_ids that start with the original edge_id followed by "_"
    new_edges1 <- graph1[grep(paste0("^", edge_id, "_"), graph1$edge_id), ]
    
    # Find new edges in graph2 that replaced this edge
    new_edges2 <- graph2[grep(paste0("^", edge_id, "_"), graph2$edge_id), ]
    
    # Skip if no matching edges found (could happen if the point was very close to a vertex)
    if (nrow(new_edges1) == 0 || nrow(new_edges2) == 0) next
    
    # Calculate the ratios for the new edges
    new_d_ratios1 <- new_edges1$d_weighted / new_edges1$d
    new_time_ratios1 <- new_edges1$time_weighted / new_edges1$time
    
    new_d_ratios2 <- new_edges2$d_weighted / new_edges2$d
    new_time_ratios2 <- new_edges2$time_weighted / new_edges2$time
    
    # The ratios for the new edges should be similar to the original edge
    # Use mean to account for small variations due to floating point arithmetic
    expect_equal (orig_d_ratios[i], mean(new_d_ratios1), tolerance = 0.01)
    expect_equal (orig_time_ratios[i], mean(new_time_ratios1), tolerance = 0.01)
    
    expect_equal (orig_d_ratios[i], mean(new_d_ratios2), tolerance = 0.01)
    expect_equal (orig_time_ratios[i], mean(new_time_ratios2), tolerance = 0.01)
    
    # Also check that both functions produce similar results
    expect_equal (mean(new_d_ratios1), mean(new_d_ratios2), tolerance = 0.01)
    expect_equal (mean(new_time_ratios1), mean(new_time_ratios2), tolerance = 0.01)
  }
  
  # Print some diagnostic information
  cat("\nEdge property preservation test:\n")
  cat("Original edges to split:", nrow(edges_to_split), "\n")
  cat("Original d_weighted/d ratios:", mean(orig_d_ratios), "\n")
  cat("Original time_weighted/time ratios:", mean(orig_time_ratios), "\n")
})

test_that ("add_nodes_to_graph_by_edge handles dist_tol parameter correctly", {
  
  # Load a sample graph
  graph <- weight_streetnet (hampi, wt_profile = "foot")
  verts <- dodgr_vertices (graph)
  
  # Create a set of points that will be placed very close to vertices
  # to test the dist_tol parameter
  set.seed (5)
  npts <- 10
  
  # Get some random vertices from the graph
  sample_verts <- verts[sample(nrow(verts), npts), ]
  
  # Create points that are very close to these vertices (within 1e-5 units)
  xy <- data.frame(
    x = sample_verts$x + rnorm(npts, 0, 1e-5),
    y = sample_verts$y + rnorm(npts, 0, 1e-5)
  )
  
  # Test with different dist_tol values
  
  # With a small tolerance, points should be treated as separate
  small_tol <- 1e-6
  graph1_small <- add_nodes_to_graph(graph, xy, dist_tol = small_tol)
  graph2_small <- add_nodes_to_graph_by_edge(graph, xy, dist_tol = small_tol)
  
  # With a larger tolerance, points should be merged with existing vertices
  large_tol <- 1
  graph1_large <- add_nodes_to_graph(graph, xy, dist_tol = large_tol)
  graph2_large <- add_nodes_to_graph_by_edge(graph, xy, dist_tol = large_tol)
  
  # Compare results
  
  # With small tolerance, both functions should add more edges
  expect_true(nrow(graph1_small) > nrow(graph))
  expect_true(nrow(graph2_small) > nrow(graph))
  
  # With large tolerance, both functions should add fewer edges
  # compared to the small tolerance case
  expect_true(nrow(graph1_large) <= nrow(graph1_small))
  expect_true(nrow(graph2_large) <= nrow(graph2_small))
  
  # The edge-based approach should be consistent with the original function
  # in how it handles the tolerance parameter
  small_diff_ratio <- abs(nrow(graph1_small) - nrow(graph2_small)) / nrow(graph1_small)
  large_diff_ratio <- abs(nrow(graph1_large) - nrow(graph2_large)) / nrow(graph1_large)
  
  # The difference between the two approaches should be similar regardless of tolerance
  expect_true(small_diff_ratio < 0.2) # Allow up to 20% difference
  expect_true(large_diff_ratio < 0.2)
  
  # Print diagnostic information
  cat("\ndist_tol parameter comparison:\n")
  cat("Original graph edges:", nrow(graph), "\n")
  cat(glue::glue("Small tolerance {small_tol}:\n"))
  cat("  add_nodes_to_graph edges:", nrow(graph1_small), "\n")
  cat("  add_nodes_to_graph_by_edge edges:", nrow(graph2_small), "\n")
  cat(glue::glue("Large tolerance {large_tol}:\n"))
  cat("  add_nodes_to_graph edges:", nrow(graph1_large), "\n")
  cat("  add_nodes_to_graph_by_edge edges:", nrow(graph2_large), "\n")
  cat("Edge difference ratio (small tolerance):", small_diff_ratio, "\n")
  cat("Edge difference ratio (large tolerance):", large_diff_ratio, "\n")
})

test_that("add_nodes_to_graph_by_edge maintains consistent coordinates by ID", {
  
  # Load a sample graph
  graph <- weight_streetnet(hampi, wt_profile = "foot")
  verts <- dodgr_vertices(graph)

  # Create a set of points
  set.seed(42)
  npts <- 2
  xy <- data.frame(
    x = min(verts$x) + runif(npts) * diff(range(verts$x)),
    y = min(verts$y) + runif(npts) * diff(range(verts$y))
  )
  
  # Create custom IDs for the points
  xy_id <- paste0("point_", seq_len(nrow(xy)))
  
  # Process with add_nodes_to_graph_by_edge using custom IDs
  graph_with_nodes <- add_nodes_to_graph_by_edge(graph, xy, xy_id = xy_id)
  
  # Get coordinates for this ID
  inconsistent_coords <- graph_with_nodes %>%
    dplyr::reframe(id=c(from_id, to_id), lon=c(from_lon, to_lon), lat=c(from_lat, to_lat))%>%
    distinct(id, lon, lat)%>%
    count(id)%>%
    filter(n>1)
  expect_equal(nrow(inconsistent_coords),0)
})
