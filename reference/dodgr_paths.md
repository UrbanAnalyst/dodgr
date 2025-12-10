# Calculate lists of pair-wise shortest paths between points.

Calculate lists of pair-wise shortest paths between points.

## Usage

``` r
dodgr_paths(
  graph,
  from,
  to,
  vertices = TRUE,
  pairwise = FALSE,
  heap = "BHeap",
  quiet = TRUE
)
```

## Arguments

- graph:

  `data.frame` or equivalent object representing the network graph (see
  Details)

- from:

  Vector or matrix of points **from** which route paths are to be
  calculated (see Details)

- to:

  Vector or matrix of points **to** which route paths are to be
  calculated (see Details)

- vertices:

  If `TRUE`, return lists of lists of vertices for each path, otherwise
  return corresponding lists of edge numbers from `graph`.

- pairwise:

  If `TRUE`, calculate paths only between the ordered pairs of `from`
  and `to`. In this case, each of these must be the same length, and the
  output will contain paths the i-th members of each, and thus also be
  of that length.

- heap:

  Type of heap to use in priority queue. Options include Fibonacci Heap
  (default; `FHeap`), Binary Heap (`BHeap`), `Radix`, Trinomial Heap
  (`TriHeap`), Extended Trinomial Heap (`TriHeapExt`, and 2-3 Heap
  (`Heap23`).

- quiet:

  If `FALSE`, display progress messages on screen.

## Value

List of list of paths tracing all connections between nodes such that if
`x <- dodgr_paths (graph, from, to)`, then the path between `from[i]`
and `to[j]` is `x [[i]] [[j]]`. Each individual path is then a vector of
integers indexing into the rows of `graph` if `vertices = FALSE`, or
into the rows of `dodgr_vertices (graph)` if `vertices = TRUE`.

## Note

`graph` must minimally contain four columns of `from`, `to`, `dist`. If
an additional column named `weight` or `wt` is present, shortest paths
are calculated according to values specified in that column; otherwise
according to `dist` values. Either way, final distances between `from`
and `to` points are calculated according to values of `dist`. That is,
paths between any pair of points will be calculated according to the
minimal total sum of `weight` values (if present), while reported
distances will be total sums of `dist` values.

The `from` and `to` columns of `graph` may be either single columns of
numeric or character values specifying the numbers or names of graph
vertices, or combinations to two columns specifying geographical
(longitude and latitude) coordinates. In the latter case, almost any
sensible combination of names will be accepted (for example,
`fromx, fromy`, `from_x, from_y`, or `fr_lat, fr_lon`.)

`from` and `to` values can be either two-column matrices of equivalent
of longitude and latitude coordinates, or else single columns precisely
matching node numbers or names given in `graph$from` or `graph$to`. If
`to` is missing, pairwise distances are calculated between all points
specified in `from`. If neither `from` nor `to` are specified, pairwise
distances are calculated between all nodes in `graph`.

## See also

Other distances:
[`dodgr_distances()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_distances.md),
[`dodgr_dists()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_dists.md),
[`dodgr_dists_categorical()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_dists_categorical.md),
[`dodgr_dists_nearest()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_dists_nearest.md),
[`dodgr_times()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_times.md)

## Examples

``` r
graph <- weight_streetnet (hampi)
from <- sample (graph$from_id, size = 100)
to <- sample (graph$to_id, size = 50)
dp <- dodgr_paths (graph, from = from, to = to)
# dp is a list with 100 items, and each of those 100 items has 30 items, each
# of which is a single path listing all vertiex IDs as taken from `graph`.

# it is also possible to calculate paths between pairwise start and end
# points
from <- sample (graph$from_id, size = 5)
to <- sample (graph$to_id, size = 5)
dp <- dodgr_paths (graph, from = from, to = to, pairwise = TRUE)
# dp is a list of 5 items, each of which just has a single path between each
# pairwise from and to point.
```
