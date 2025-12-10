# Insert new nodes into a graph, breaking edges at point of nearest intersection.

Note that this routine presumes graphs to be `dodgr_streetnet` object,
with geographical coordinates.

## Usage

``` r
add_nodes_to_graph(graph, xy, dist_tol = 0.000001, intersections_only = FALSE)
```

## Arguments

- graph:

  A `dodgr` graph with spatial coordinates, such as a `dodgr_streetnet`
  object.

- xy:

  coordinates of points to be matched to the vertices, either as matrix
  or sf-formatted `data.frame`.

- dist_tol:

  Only insert new nodes if they are further from existing nodes than
  this distance, expressed in units of the distance column of `graph`.

- intersections_only:

  If `FALSE`

## Value

A modified version of `graph`, with additional edges formed by breaking
previous edges at nearest perpendicular intersections with the points,
`xy`.

## Details

This inserts new nodes by extending lines from each input point to the
edge with the closest point of perpendicular intersection. That edge is
then split at that point of intersection, creating two new edges (or
four for directed edges). If `intersections_only = FALSE` (default),
then additional edges are inserted from those intersection points to the
input points. If `intersections_only = TRUE`, then nodes are added by
splitting graph edges at points of nearest perpendicular intersection,
without adding additional edges out to the actual input points.

In the former case, the properties of those new edges, such as distance
and time weightings, are inherited from the edges which are intersected,
and may need to be manually modified after calling this function.

## See also

Other match:
[`match_points_to_graph()`](https://UrbanAnalyst.github.io/dodgr/reference/match_points_to_graph.md),
[`match_points_to_verts()`](https://UrbanAnalyst.github.io/dodgr/reference/match_points_to_verts.md),
[`match_pts_to_graph()`](https://UrbanAnalyst.github.io/dodgr/reference/match_pts_to_graph.md),
[`match_pts_to_verts()`](https://UrbanAnalyst.github.io/dodgr/reference/match_pts_to_verts.md)

## Examples

``` r
graph <- weight_streetnet (hampi, wt_profile = "foot")
dim (graph)
#> [1] 6813   15

verts <- dodgr_vertices (graph)
set.seed (2)
npts <- 10
xy <- data.frame (
    x = min (verts$x) + runif (npts) * diff (range (verts$x)),
    y = min (verts$y) + runif (npts) * diff (range (verts$y))
)

graph <- add_nodes_to_graph (graph, xy)
dim (graph) # more edges than original
#> [1] 6863   15
```
