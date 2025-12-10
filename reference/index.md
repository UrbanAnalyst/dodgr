# Package index

## Main distance and time functions

- [`dodgr_distances()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_distances.md)
  : Calculate matrix of pair-wise distances between points.
- [`dodgr_dists()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_dists.md)
  : Calculate matrix of pair-wise distances between points.
- [`dodgr_dists_categorical()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_dists_categorical.md)
  : Cumulative distances along different edge categories
- [`dodgr_dists_nearest()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_dists_nearest.md)
  : Calculate vector of shortest distances from a series of 'from'
  points to nearest one of series of 'to' points.
- [`dodgr_paths()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_paths.md)
  : Calculate lists of pair-wise shortest paths between points.
- [`dodgr_times()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_times.md)
  : Calculate matrix of pair-wise travel times between points.

## Main flow functions

- [`dodgr_flows_aggregate()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_flows_aggregate.md)
  : Aggregate flows throughout a network.
- [`dodgr_flows_disperse()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_flows_disperse.md)
  : Aggregate flows dispersed from each point in a network.
- [`dodgr_flows_si()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_flows_si.md)
  : Aggregate flows throughout a network using a spatial interaction
  model.

## Main iso functions

- [`dodgr_isochrones()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_isochrones.md)
  : Calculate isochrone contours from specified points.
- [`dodgr_isodists()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_isodists.md)
  : Calculate isodistance contours from specified points.
- [`dodgr_isoverts()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_isoverts.md)
  : Calculate isodistance or isochrone vertices from specified points.

## Functions to Obtain Graphs

- [`dodgr_streetnet()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_streetnet.md)
  :

  Extract a street network in sf-format for a given location.

- [`dodgr_streetnet_geodesic()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_streetnet_geodesic.md)
  : Force weight_streetnet to use geodesic distances.

- [`dodgr_streetnet_sc()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_streetnet_sc.md)
  :

  Extract a street network in silicate-format for a given location.

- [`weight_railway()`](https://UrbanAnalyst.github.io/dodgr/reference/weight_railway.md)
  : Weight a network for routing along railways.

- [`weight_streetnet()`](https://UrbanAnalyst.github.io/dodgr/reference/weight_streetnet.md)
  : Weight a street network according to a specified weighting profile.

## Functions to Modify Graphs

- [`dodgr_components()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_components.md)
  : Identify connected components of graph.
- [`dodgr_contract_graph()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_contract_graph.md)
  : Contract graph to junction vertices only.
- [`dodgr_uncontract_graph()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_uncontract_graph.md)
  : Re-expand a contracted graph.

## Functions to Convert Graphs

- [`dodgr_deduplicate_graph()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_deduplicate_graph.md)
  : Deduplicate edges in a graph

- [`dodgr_to_igraph()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_to_igraph.md)
  :

  Convert a `dodgr` graph to an igraph.

- [`dodgr_to_sf()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_to_sf.md)
  :

  Convert a `dodgr` graph into an equivalent sf object.

- [`dodgr_to_sfc()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_to_sfc.md)
  :

  Convert a `dodgr` graph into an equivalent
  [`sf::sfc`](https://r-spatial.github.io/sf/reference/sfc.html) object.

- [`dodgr_to_tidygraph()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_to_tidygraph.md)
  :

  Convert a `dodgr` graph to an tidygraph.

- [`igraph_to_dodgr()`](https://UrbanAnalyst.github.io/dodgr/reference/igraph_to_dodgr.md)
  :

  Convert a igraph network to an equivalent `dodgr` representation.

## Graph Centrality

- [`dodgr_centrality()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_centrality.md)
  : Calculate betweenness centrality for a 'dodgr' network.
- [`estimate_centrality_threshold()`](https://UrbanAnalyst.github.io/dodgr/reference/estimate_centrality_threshold.md)
  : Estimate a value for the 'dist_threshold' parameter of the
  dodgr_centrality function.
- [`estimate_centrality_time()`](https://UrbanAnalyst.github.io/dodgr/reference/estimate_centrality_time.md)
  : Estimate time required for a planned centrality calculation.

## Matching Points to Graphs

- [`add_nodes_to_graph()`](https://UrbanAnalyst.github.io/dodgr/reference/add_nodes_to_graph.md)
  : Insert new nodes into a graph, breaking edges at point of nearest
  intersection.
- [`match_points_to_graph()`](https://UrbanAnalyst.github.io/dodgr/reference/match_points_to_graph.md)
  : Alias for match_pts_to_graph
- [`match_points_to_verts()`](https://UrbanAnalyst.github.io/dodgr/reference/match_points_to_verts.md)
  : Alias for match_pts_to_verts
- [`match_pts_to_graph()`](https://UrbanAnalyst.github.io/dodgr/reference/match_pts_to_graph.md)
  : Match spatial points to the edges of a spatial graph.
- [`match_pts_to_verts()`](https://UrbanAnalyst.github.io/dodgr/reference/match_pts_to_verts.md)
  : Match spatial points to the vertices of a spatial graph.

## Miscellaneous Functions

- [`compare_heaps()`](https://UrbanAnalyst.github.io/dodgr/reference/compare_heaps.md)
  : Compare timings of different sort heaps for a given input graph.

- [`dodgr_flowmap()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_flowmap.md)
  :

  Create a map of `dodgr` flows.

- [`dodgr_full_cycles()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_full_cycles.md)
  : Calculate fundamental cycles on a FULL (that is, non-contracted)
  graph.

- [`dodgr_fundamental_cycles()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_fundamental_cycles.md)
  : Calculate fundamental cycles in a graph.

- [`dodgr_insert_vertex()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_insert_vertex.md)
  : Insert a new node or vertex into a network

- [`dodgr_sample()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_sample.md)
  : Sample a random but connected sub-component of a graph

- [`dodgr_sflines_to_poly()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_sflines_to_poly.md)
  :

  Convert sf `LINESTRING` objects to `POLYGON` objects representing all
  fundamental cycles within the `LINESTRING` objects.

- [`dodgr_vertices()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_vertices.md)
  : Extract vertices of graph, including spatial coordinates if
  included.

- [`merge_directed_graph()`](https://UrbanAnalyst.github.io/dodgr/reference/merge_directed_graph.md)
  : Merge directed edges into equivalent undirected edges.

- [`summary(`*`<dodgr_dists_categorical>`*`)`](https://UrbanAnalyst.github.io/dodgr/reference/summary.dodgr_dists_categorical.md)
  : Transform a result from dodgr_dists_categorical to summary
  statistics

- [`write_dodgr_wt_profile()`](https://UrbanAnalyst.github.io/dodgr/reference/write_dodgr_wt_profile.md)
  :

  Write `dodgr` weighting profiles to local file.

## Save, Load, and Caching

- [`clear_dodgr_cache()`](https://UrbanAnalyst.github.io/dodgr/reference/clear_dodgr_cache.md)
  :

  Remove cached versions of `dodgr` graphs.

- [`dodgr_cache_off()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_cache_off.md)
  : Turn off all dodgr caching in current session.

- [`dodgr_cache_on()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_cache_on.md)
  : Turn on all dodgr caching in current session.

- [`dodgr_load_streetnet()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_load_streetnet.md)
  : Load a street network previously saved with dodgr_save_streetnet.

- [`dodgr_save_streetnet()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_save_streetnet.md)
  : Save a weighted streetnet to a local file

## Package Data

- [`hampi`](https://UrbanAnalyst.github.io/dodgr/reference/hampi.md) :
  Sample street network from Hampi, India.
- [`os_roads_bristol`](https://UrbanAnalyst.github.io/dodgr/reference/os_roads_bristol.md)
  : Sample street network from Bristol, U.K.
- [`weighting_profiles`](https://UrbanAnalyst.github.io/dodgr/reference/weighting_profiles.md)
  : Weighting profiles used to route different modes of transport.

## The ‘dodgr’ package

- [`dodgr-package`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr.md)
  [`dodgr`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr.md) :
  Distances On Directed GRaphs ("dodgr")
