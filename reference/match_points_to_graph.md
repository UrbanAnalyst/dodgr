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
#>      geom_num edge_id    from_id from_lon from_lat      to_id   to_lon   to_lat
#> 5247      157    5247 5351564647 76.43545 15.34650 3921522971 76.43542 15.34657
#> 2835       72    2835 1398747966 76.44911 15.35290 1398747940 76.44880 15.35303
#> 3405       89    3405 2588146138 76.42439 15.31886 2588146016 76.42477 15.31911
#> 3395       89    3395 2588119056 76.42341 15.31717 2588146107 76.42349 15.31746
#> 3583       90    3583 7794286109 76.43558 15.33028 2588155764 76.43554 15.33031
#> 3179       82    3179 2195425010 76.46729 15.35314 2195425005 76.46774 15.35302
#> 3781       90    3781 7794251342 76.44447 15.32376 7794251341 76.44445 15.32392
#> 3017       82    3017 2195424965 76.45517 15.35171 2195424960 76.45557 15.35149
#> 3569       90    3569 2588155753 76.43665 15.33008 2588155735 76.43640 15.33009
#> 247         2     247 2398957753 76.47679 15.32222 2398957756 76.47659 15.32304
#>              d d_weighted      highway    way_id component       time
#> 5247  9.413905   9.909373        track 388983220         2   6.778011
#> 2835 36.189595  45.236994 unclassified 126094049         2  26.056508
#> 3405 48.824353  61.030441 unclassified 252786290         1  35.153534
#> 3395 33.008446  41.260558 unclassified 252786290         1  23.766081
#> 3583  4.767803   5.018740        track 252787544         1   3.432818
#> 3179 50.815069  63.518836 unclassified 209318354         2  36.586849
#> 3781 18.415651  19.384896        track 252787544         1  13.259269
#> 3017 49.280905  61.601131 unclassified 209318354         2  35.482251
#> 3569 27.752699  29.213367        track 252787544         1  19.981943
#> 247  93.004062  93.004062         path  30643853         1 167.407312
#>      time_weighted
#> 5247      7.134749
#> 2835     32.570636
#> 3405     43.941918
#> 3395     29.707601
#> 3583      3.613493
#> 3179     45.733562
#> 3781     13.957125
#> 3017     44.352814
#> 3569     21.033625
#> 247     167.407312
```
