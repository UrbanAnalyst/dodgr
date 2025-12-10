# Alias for [match_pts_to_verts](https://UrbanAnalyst.github.io/dodgr/reference/match_pts_to_verts.md)

The
[match_pts_to_graph](https://UrbanAnalyst.github.io/dodgr/reference/match_pts_to_graph.md)
function matches points to the nearest edge based on geometric
intersections; this function only matches to the nearest vertex based on
point-to-point distances.

## Usage

``` r
match_points_to_verts(verts, xy, connected = FALSE)
```

## Arguments

- verts:

  A `data.frame` of vertices obtained from `dodgr_vertices(graph)`.

- xy:

  coordinates of points to be matched to the vertices, either as matrix
  or sf-formatted `data.frame`.

- connected:

  Should points be matched to the same (largest) connected component of
  graph? If `FALSE` and these points are to be used for a `dodgr`
  routing routine
  ([dodgr_dists](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_dists.md),
  [dodgr_paths](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_paths.md),
  or
  [dodgr_flows_aggregate](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_flows_aggregate.md)),
  then results may not be returned if points are not part of the same
  connected component. On the other hand, forcing them to be part of the
  same connected component may decrease the spatial accuracy of
  matching.

## Value

A vector index into verts

## See also

Other match:
[`add_nodes_to_graph()`](https://UrbanAnalyst.github.io/dodgr/reference/add_nodes_to_graph.md),
[`match_points_to_graph()`](https://UrbanAnalyst.github.io/dodgr/reference/match_points_to_graph.md),
[`match_pts_to_graph()`](https://UrbanAnalyst.github.io/dodgr/reference/match_pts_to_graph.md),
[`match_pts_to_verts()`](https://UrbanAnalyst.github.io/dodgr/reference/match_pts_to_verts.md)

## Examples

``` r
net <- weight_streetnet (hampi, wt_profile = "foot")
verts <- dodgr_vertices (net)
# Then generate some random points to match to graph
npts <- 10
xy <- data.frame (
    x = min (verts$x) + runif (npts) * diff (range (verts$x)),
    y = min (verts$y) + runif (npts) * diff (range (verts$y))
)
pts <- match_pts_to_verts (verts, xy)
pts # an index into verts
#>  [1] 1704 1174 3055 1450 1704 2349 3053 3053   22  654
pts <- verts$id [pts]
pts # names of those vertices
#>  [1] "2588119056" "7799711008" "1388481924" "1398748033" "2588119056"
#>  [6] "2398957595" "7769190961" "7769190961" "2398957995" "5351515769"
```
