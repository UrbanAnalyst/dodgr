# Alias for [match_pts_to_graph](https://UrbanAnalyst.github.io/dodgr/reference/match_pts_to_graph.md)

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
match_points_to_graph(graph, xy, connected = FALSE, distances = FALSE)
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
[`match_points_to_verts()`](https://UrbanAnalyst.github.io/dodgr/reference/match_points_to_verts.md),
[`match_pts_to_graph()`](https://UrbanAnalyst.github.io/dodgr/reference/match_pts_to_graph.md),
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
#>        geom_num edge_id    from_id from_lon from_lat      to_id   to_lon
#> 3403         89    3403 7793366198 76.42424 15.31863 2588146138 76.42439
#> 5867        183    5867 1376768565 76.48143 15.31815 2398957701 76.48100
#> 1713         50    1713 2632626796 76.46953 15.34602 2632626792 76.46943
#> 3395         89    3395 2588119056 76.42341 15.31717 2588146107 76.42349
#> 3395.1       89    3395 2588119056 76.42341 15.31717 2588146107 76.42349
#> 6251        203    6251 6597298435 76.42010 15.35136 6597300510 76.42343
#> 187           2     187 2398957685 76.47764 15.31653 2398957688 76.47753
#> 1097         35    1097 1204772807 76.44527 15.34382 1204772888 76.44661
#> 4155        100    4155 2627461399 76.46968 15.34330 2627461395 76.46979
#> 4253        107    4253 2627486261 76.47163 15.35948 2627486260 76.47137
#>          to_lat         d d_weighted      highway    way_id component
#> 3403   15.31886  31.37821   39.22276 unclassified 252786290         1
#> 5867   15.31770  68.04148  136.08296      primary 652570479         1
#> 1713   15.34568  38.77995   43.08884  residential  84014148         2
#> 3395   15.31746  33.00845   41.26056 unclassified 252786290         1
#> 3395.1 15.31746  33.00845   41.26056 unclassified 252786290         1
#> 6251   15.35101 359.82586  719.65173      primary 835018468         2
#> 187    15.31697  50.09514   50.09514         path  30643853         1
#> 1097   15.34379 144.87782  289.75564      primary  53658844         2
#> 4155   15.34326  12.19799   12.19799         path 257144184         2
#> 4253   15.35960  31.22245   32.86573        track 257147579         2
#>              time time_weighted
#> 3403    22.592310     28.240387
#> 5867    48.989866     97.979732
#> 1713    27.921568     31.023964
#> 3395    23.766081     29.707601
#> 3395.1  23.766081     29.707601
#> 6251   259.074622    518.149244
#> 187     90.171255     90.171255
#> 1097   104.312029    208.624057
#> 4155     8.782554      8.782554
#> 4253    22.480162     23.663329
```
