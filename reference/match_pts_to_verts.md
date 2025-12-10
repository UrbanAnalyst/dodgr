# Match spatial points to the vertices of a spatial graph.

The
[match_pts_to_graph](https://UrbanAnalyst.github.io/dodgr/reference/match_pts_to_graph.md)
function matches points to the nearest edge based on geometric
intersections; this function only matches to the nearest vertex based on
point-to-point distances.

## Usage

``` r
match_pts_to_verts(verts, xy, connected = FALSE)
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
[`match_points_to_verts()`](https://UrbanAnalyst.github.io/dodgr/reference/match_points_to_verts.md),
[`match_pts_to_graph()`](https://UrbanAnalyst.github.io/dodgr/reference/match_pts_to_graph.md)

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
#>  [1] 2926  389 3076 1704 2569 1460 3052 1044 1813 2926
pts <- verts$id [pts]
pts # names of those vertices
#>  [1] "676635868"  "571423425"  "1204772868" "2588119056" "3921522991"
#>  [6] "1427080597" "7769271419" "1376768725" "7794286099" "676635868" 
```
