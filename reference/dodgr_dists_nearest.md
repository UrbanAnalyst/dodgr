# Calculate vector of shortest distances from a series of 'from' points to nearest one of series of 'to' points.

Calculate vector of shortest distances from a series of 'from' points to
nearest one of series of 'to' points.

## Usage

``` r
dodgr_dists_nearest(
  graph,
  from = NULL,
  to = NULL,
  shortest = TRUE,
  heap = "BHeap",
  quiet = TRUE
)
```

## Arguments

- graph:

  `data.frame` or equivalent object representing the network graph (see
  Notes). For `dodgr` street networks, this may be a network derived
  from either sf or silicate ("sc") data, generated with
  [weight_streetnet](https://UrbanAnalyst.github.io/dodgr/reference/weight_streetnet.md).

  The `from` and `to` columns of `graph` may be either single columns of
  numeric or character values specifying the numbers or names of graph
  vertices, or combinations to two columns specifying geographical
  (longitude and latitude,) coordinates. In the latter case, almost any
  sensible combination of names will be accepted (for example,
  `fromx, fromy`, `from_x, from_y`, or `fr_lat, fr_lon`.)

  Note that longitude and latitude values are always interpreted in
  'dodgr' to be in EPSG:4326 / WSG84 coordinates. Any other kinds of
  coordinates should first be reprojected to EPSG:4326 before submitting
  to any 'dodgr' routines.

  See further information in Details.

- from:

  Vector or matrix of points **from** which route distances are to be
  calculated, specified as one of the following:

  - Single character vector precisely matching node numbers or names
    given in `graph$from` or `graph$to`.

  - Single vector of integer-ish values, in which case these will be
    presumed to specify indices into
    [dodgr_vertices](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_vertices.md),
    and NOT to correspond to values in the 'from' or 'to' columns of the
    graph. See the example below for a demonstration.

  - Matrix or equivalent of longitude and latitude coordinates, in which
    case these will be matched on to the nearest coordinates of 'from'
    and 'to' points in the graph.

- to:

  Vector or matrix of points **to** which route distances are to be
  calculated. If `to` is `NULL`, pairwise distances will be calculated
  from all `from` points to all other nodes in `graph`. If both `from`
  and `to` are `NULL`, pairwise distances are calculated between all
  nodes in `graph`.

- shortest:

  If `FALSE`, calculate distances along the *fastest* rather than
  shortest routes. For street networks produced with
  [weight_streetnet](https://UrbanAnalyst.github.io/dodgr/reference/weight_streetnet.md),
  distances may also be calculated along the *fastest* routes with the
  `shortest = FALSE` option. Graphs must in this case have columns of
  `time` and `time_weighted`. Note that the fastest routes will only be
  approximate when derived from sf-format data generated with the
  osmdata function `osmdata_sf()`, and will be much more accurate when
  derived from `sc`-format data generated with `osmdata_sc()`. See
  [weight_streetnet](https://UrbanAnalyst.github.io/dodgr/reference/weight_streetnet.md)
  for details.

- heap:

  Type of heap to use in priority queue. Options include Fibonacci Heap
  (default; `FHeap`), Binary Heap (`BHeap`),
  `Trinomial Heap (`TriHeap`), Extended Trinomial Heap (`TriHeapExt`, and 2-3 Heap (`Heap23\`).

- quiet:

  If `FALSE`, display progress messages on screen.

## Value

Vector of distances, one element for each 'from' point giving the
distance to the nearest 'to' point.

## Note

`graph` must minimally contain three columns of `from`, `to`, `dist`. If
an additional column named `weight` or `wt` is present, shortest paths
are calculated according to values specified in that column; otherwise
according to `dist` values. Either way, final distances between `from`
and `to` points are calculated by default according to values of `dist`.
That is, paths between any pair of points will be calculated according
to the minimal total sum of `weight` values (if present), while reported
distances will be total sums of `dist` values.

For street networks produced with
[weight_streetnet](https://UrbanAnalyst.github.io/dodgr/reference/weight_streetnet.md),
distances may also be calculated along the *fastest* routes with the
`shortest = FALSE` option. Graphs must in this case have columns of
`time` and `time_weighted`. Note that the fastest routes will only be
approximate when derived from sf-format data generated with the osmdata
function `osmdata_sf()`, and will be much more accurate when derived
from `sc`-format data generated with `osmdata_sc()`. See
[weight_streetnet](https://UrbanAnalyst.github.io/dodgr/reference/weight_streetnet.md)
for details.

The `from` and `to` columns of `graph` may be either single columns of
numeric or character values specifying the numbers or names of graph
vertices, or combinations to two columns specifying geographical
(longitude and latitude) coordinates. In the latter case, almost any
sensible combination of names will be accepted (for example,
`fromx, fromy`, `from_x, from_y`, or `fr_lat, fr_lon`.)

`from` and `to` values can be either two-column matrices or equivalent
of longitude and latitude coordinates, or else single columns precisely
matching node numbers or names given in `graph$from` or `graph$to`. If
`to` is `NULL`, pairwise distances are calculated from all `from` points
to all other nodes in `graph`. If both `from` and `to` are `NULL`,
pairwise distances are calculated between all nodes in `graph`.

Calculations are always calculated in parallel, using multiple threads.

## See also

Other distances:
[`dodgr_distances()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_distances.md),
[`dodgr_dists()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_dists.md),
[`dodgr_dists_categorical()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_dists_categorical.md),
[`dodgr_paths()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_paths.md),
[`dodgr_times()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_times.md)

## Examples

``` r
# A simple graph
graph <- data.frame (
    from = c ("A", "B", "B", "B", "C", "C", "D", "D"),
    to = c ("B", "A", "C", "D", "B", "D", "C", "A"),
    d = c (1, 2, 1, 3, 2, 1, 2, 1)
)
dodgr_dists (graph)
#>   A B C D
#> A 0 1 2 3
#> B 2 0 1 2
#> C 2 2 0 1
#> D 1 2 2 0

# A larger example from the included [hampi()] data.
graph <- weight_streetnet (hampi)
from <- sample (graph$from_id, size = 100)
to <- sample (graph$to_id, size = 50)
d <- dodgr_dists (graph, from = from, to = to)
# d is a 100-by-50 matrix of distances between `from` and `to`

if (FALSE) { # \dontrun{
# a more complex street network example, thanks to @chrijo; see
# https://github.com/UrbanAnalyst/dodgr/issues/47

xy <- rbind (
    c (7.005994, 51.45774), # limbeckerplatz 1 essen germany
    c (7.012874, 51.45041)
) # hauptbahnhof essen germany
xy <- data.frame (lon = xy [, 1], lat = xy [, 2])
essen <- dodgr_streetnet (pts = xy, expand = 0.2, quiet = FALSE)
graph <- weight_streetnet (essen, wt_profile = "foot")
d <- dodgr_dists (graph, from = xy, to = xy)
# First reason why this does not work is because the graph has multiple,
# disconnected components.
table (graph$component)
# reduce to largest connected component, which is always number 1
graph <- graph [which (graph$component == 1), ]
d <- dodgr_dists (graph, from = xy, to = xy)
# should work, but even then note that
table (essen$level)
# There are parts of the network on different building levels (because of
# shopping malls and the like). These may or may not be connected, so it may
# be necessary to filter out particular levels
index <- which (!(essen$level == "-1" | essen$level == "1")) # for example
library (sf) # needed for following sub-select operation
essen <- essen [index, ]
graph <- weight_streetnet (essen, wt_profile = "foot")
graph <- graph [which (graph$component == 1), ]
d <- dodgr_dists (graph, from = xy, to = xy)
} # }
```
