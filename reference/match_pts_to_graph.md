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
#> 6309      205    6309 7799710769 76.47825 15.32601 7799710767 76.47857 15.32603
#> 213         2     213 2398957706 76.47670 15.31802 2398957709 76.47663 15.31828
#> 3689       90    3689 2588155728 76.44036 15.32559 2588155748 76.44051 15.32549
#> 2447       70    2447 1376768866 76.47733 15.33116 1376768379 76.47788 15.33140
#> 4253      107    4253 2627486261 76.47163 15.35948 2627486260 76.47137 15.35960
#> 3109       82    3109 2195425012 76.46200 15.35325 2195425016 76.46177 15.35348
#> 3467       89    3467 5351719173 76.43085 15.31804 5351719172 76.43136 15.31796
#> 3469       89    3469 5351719172 76.43136 15.31796 2588146136 76.43231 15.31795
#> 6201      203    6201 7769271419 76.37263 15.34499 7769190961 76.38203 15.34708
#> 5585      170    5585 5974426284 76.44659 15.31570 5974426290 76.44659 15.31617
#>               d d_weighted      highway    way_id component      time
#> 6309   34.53666   36.35438        track 835627549         3  24.86639
#> 213    30.08419   30.08419         path  30643853         1  54.15155
#> 3689   19.99728   21.04977        track 252787544         1  14.39804
#> 2447   64.78795   64.78795         path 123463598         1  46.64732
#> 4253   31.22245   32.86573        track 257147579         2  22.48016
#> 3109   36.15702   45.19627 unclassified 209318354         2  26.03305
#> 3467   55.71128   69.63911 unclassified 252786290         1  40.11212
#> 3469  102.32710  127.90887 unclassified 252786290         1  73.67551
#> 6201 1035.68141 2071.36283      primary 835018468         2 745.69062
#> 5585   51.80895   54.53574        track 554572318         1  37.30244
#>      time_weighted
#> 6309      26.17515
#> 213       54.15155
#> 3689      15.15583
#> 2447      46.64732
#> 4253      23.66333
#> 3109      32.54131
#> 3467      50.14016
#> 3469      92.09439
#> 6201    1491.38124
#> 5585      39.26573
```
