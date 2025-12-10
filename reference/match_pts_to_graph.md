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
#>        geom_num edge_id    from_id from_lon from_lat      to_id   to_lon
#> 6201        203    6201 7769271419 76.37263 15.34499 7769190961 76.38203
#> 3395         89    3395 2588119056 76.42341 15.31717 2588146107 76.42349
#> 3617         90    3617 2588155770 76.43412 15.32851 7794286097 76.43416
#> 3395.1       89    3395 2588119056 76.42341 15.31717 2588146107 76.42349
#> 3395.2       89    3395 2588119056 76.42341 15.31717 2588146107 76.42349
#> 503          18     503  339574188 76.44575 15.33279  339574189 76.44510
#> 5919        183    5919 2398957539 76.47550 15.30880  676875267 76.47492
#> 3395.3       89    3395 2588119056 76.42341 15.31717 2588146107 76.42349
#> 6201.1      203    6201 7769271419 76.37263 15.34499 7769190961 76.38203
#> 4715        133    4715  338905103 76.46961 15.31386  338905105 76.47080
#>          to_lat          d d_weighted      highway    way_id component
#> 6201   15.34708 1035.68141 2071.36283      primary 835018468         2
#> 3395   15.31746   33.00845   41.26056 unclassified 252786290         1
#> 3617   15.32842   10.78597   11.35366        track 252787544         1
#> 3395.1 15.31746   33.00845   41.26056 unclassified 252786290         1
#> 3395.2 15.31746   33.00845   41.26056 unclassified 252786290         1
#> 503    15.33287   70.55725   74.27078        track  30704678         1
#> 5919   15.30807  102.42454  204.84908      primary 652570479         1
#> 3395.3 15.31746   33.00845   41.26056 unclassified 252786290         1
#> 6201.1 15.34708 1035.68141 2071.36283      primary 835018468         2
#> 4715   15.31338  138.14499  230.24166    secondary 327102382         1
#>              time time_weighted
#> 6201   745.690618   1491.381236
#> 3395    23.766081     29.707601
#> 3617     7.765901      8.174633
#> 3395.1  23.766081     29.707601
#> 3395.2  23.766081     29.707601
#> 503     50.801217     53.474965
#> 5919    73.745670    147.491340
#> 3395.3  23.766081     29.707601
#> 6201.1 745.690618   1491.381236
#> 4715    99.464395    165.773992
```
