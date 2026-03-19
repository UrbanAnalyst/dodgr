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
#> 5283      158    5283 7794286139 76.44353 15.34281 7794286140 76.44359 15.34272
#> 1287       35    1287 5351515769 76.48789 15.35570 1206252310 76.48878 15.35523
#> 3395       89    3395 2588119056 76.42341 15.31717 2588146107 76.42349 15.31746
#> 3547       89    3547 2588146120 76.44917 15.32131 2588146091 76.45042 15.32141
#> 909        34     909 5974503437 76.45142 15.31516 2398957668 76.45161 15.31592
#> 2839       72    2839 1398748037 76.44852 15.35307 1398748016 76.44688 15.35315
#> 6199      203    6199 6025347259 76.37261 15.34499 7769271419 76.37263 15.34499
#> 6247      203    6247 1204772780 76.40951 15.35236 1204772868 76.41809 15.35157
#> 6225      203    6225 1204772661 76.39213 15.35169 1204772759 76.39828 15.35260
#> 6255      203    6255 1388482473 76.42580 15.35076 1204772830 76.42729 15.35071
#>               d  d_weighted      highway    way_id component       time
#> 5283  11.992240   12.623410        track 388983289         2   8.634413
#> 1287 108.839929  217.679858      primary  53658844         2  78.364749
#> 3395  33.008446   41.260558 unclassified 252786290         1  23.766081
#> 3547 134.664487  168.330608 unclassified 252786290         1  96.958430
#> 909   86.572237  144.287062    secondary  53626074         1  62.332011
#> 2839 176.004015  220.005019 unclassified 126094049         2 126.722891
#> 6199   2.646247    5.292493      primary 835018468         2   1.905298
#> 6247 925.873570 1851.747140      primary 835018468         2 666.628970
#> 6225 667.359739 1334.719477      primary 835018468         2 480.499012
#> 6255 160.538138  321.076277      primary 835018468         2 115.587460
#>      time_weighted
#> 5283      9.088856
#> 1287    156.729498
#> 3395     29.707601
#> 3547    121.198038
#> 909     103.886684
#> 2839    158.403613
#> 6199      3.810595
#> 6247   1333.257941
#> 6225    960.998024
#> 6255    231.174919
```
