# Match spatial points to the edges of a spatial graph.

Match spatial points to the edges of a spatial graph, through finding
the edge with the closest perpendicular intersection. NOTE:
Intersections are calculated geometrically, and presume planar geometry.
It is up to users of projected geometrical data, such as those within a
`dodgr_streetnet` object, to ensure that either: (i) Data span an
sufficiently small area that errors from presuming planar geometry may
be ignored; or (ii) Data are re-projected to an equivalent planar
geometry prior to calling this routine.

## Usage

``` r
match_pts_to_graph(graph, xy, connected = FALSE, distances = FALSE)
```

## Arguments

- graph:

  A `dodgr` graph with spatial coordinates, such as a `dodgr_streetnet`
  object.

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

- distances:

  If `TRUE`, return a 'data.frame' object with 'index' column as
  described in return value; and additional columns with perpendicular
  distance to nearest edge in graph, and coordinates of points of
  intersection. See description of return value for details.

## Value

For `distances = FALSE` (default), a vector index matching the `xy`
coordinates to nearest edges. For bi-directional edges, only one match
is returned, and it is up to the user to identify and suitably process
matching edge pairs. For 'distances = TRUE', a 'data.frame' of four
columns:

- "index" The index of closest edges in "graph", as described above.

- "d_signed" The perpendicular distance from ech point to the nearest
  edge, with negative distances denoting points to the left of edges,
  and positive distances denoting points to the right. Distances of zero
  denote points lying precisely on the line of an edge (potentially
  including cases where nearest point of bisection lies beyond the
  actual edge).

- "x" The x-coordinate of the point of intersection.

- "y" The y-coordinate of the point of intersection.

## See also

Other match:
[`add_nodes_to_graph()`](https://UrbanAnalyst.github.io/dodgr/reference/add_nodes_to_graph.md),
[`match_points_to_graph()`](https://UrbanAnalyst.github.io/dodgr/reference/match_points_to_graph.md),
[`match_points_to_verts()`](https://UrbanAnalyst.github.io/dodgr/reference/match_points_to_verts.md),
[`match_pts_to_verts()`](https://UrbanAnalyst.github.io/dodgr/reference/match_pts_to_verts.md)

## Examples

``` r
graph <- weight_streetnet (hampi, wt_profile = "foot")
# Then generate some random points to match to graph
verts <- dodgr_vertices (graph)
npts <- 10
xy <- data.frame (
    x = min (verts$x) + runif (npts) * diff (range (verts$x)),
    y = min (verts$y) + runif (npts) * diff (range (verts$y))
)
edges <- match_pts_to_graph (graph, xy)
graph [edges, ] # The edges of the graph closest to `xy`
#>      geom_num edge_id    from_id from_lon from_lat      to_id   to_lon   to_lat
#> 3395       89    3395 2588119056 76.42341 15.31717 2588146107 76.42349 15.31746
#> 3553       90    3553 2588155738 76.43808 15.33034 7794286115 76.43804 15.33042
#> 6211      203    6211 1204772662 76.38895 15.34748 1204772772 76.38938 15.34775
#> 2887       73    2887 1427080652 76.46669 15.31934 1427080653 76.46640 15.31955
#> 2653       70    2653 4474520513 76.48875 15.33092 4474520512 76.48875 15.33085
#> 2641       70    2641 4474520519 76.48811 15.33201 4474520518 76.48828 15.33176
#> 3605       90    3605 7794286102 76.43465 15.32939 7794286101 76.43451 15.32927
#> 6285      203    6285 3921522973 76.44042 15.34841 1388482448 76.44155 15.34693
#> 5877      183    5877 1128374399 76.47924 15.31509  313796427 76.47893 15.31474
#> 6203      203    6203 7769190961 76.38203 15.34708 1204772804 76.38329 15.34736
#>               d d_weighted      highway    way_id component       time
#> 3395  33.008446  41.260558 unclassified 252786290         1  23.766081
#> 3553  10.408489  10.956305        track 252787544         1   7.494112
#> 6211  55.472027 110.944054      primary 835018468         2  39.939860
#> 2887  38.142290  38.142290         path 129323700         1  27.462449
#> 2653   8.057293   8.057293         path 123463598         1   5.801251
#> 2641  33.672123  33.672123         path 123463598         1  24.243928
#> 3605  19.575233  20.605508        track 252787544         1  14.094167
#> 6285 203.231629 406.463259      primary 835018468         2 146.326773
#> 5877  51.465621 102.931241      primary 652570479         1  37.055247
#> 6203 138.666644 277.333288      primary 835018468         2  99.839984
#>      time_weighted
#> 3395     29.707601
#> 3553      7.888539
#> 6211     79.879719
#> 2887     27.462449
#> 2653      5.801251
#> 2641     24.243928
#> 3605     14.835966
#> 6285    292.653546
#> 5877     74.110494
#> 6203    199.679967
```
