# Calculate isochrone contours from specified points.

Calculates isochrones from input `data.frame` objects (`graph`), which
must minimally contain three columns of `from`, `to`, and `t` or `time`.
If an additional column named `t_weight` or `t_wt` is present, fastest
paths are calculated according to values specified in that column, while
resultant isochrones are calculated from the `t` or `time` column. That
is, the paths tracing isochrones from any point will be calculated
according to the minimal total sum of `t_weight` values (if present),
while reported isochrones will be total sums of `time` values.

Graphs derived from Open Street Map street networks, via the
[weight_streetnet](https://UrbanAnalyst.github.io/dodgr/reference/weight_streetnet.md)
function, have columns labelled `d`, `d_weighted`, `time`, and
`time_weighted`. For these inputs, isochrones are always routed using
`t_weighted`, while final isochrones are sums of values of `t` - that
is, of un-weighted distances or times - along those paths.

Function is fully vectorized to calculate accept vectors of central
points and vectors defining multiple isochrones. Calculations use by
default parallel computation with the maximal number of available cores
or threads. This number can be reduced by specifying a value via
`RcppParallel::setThreadOptions (numThreads = <desired_number>)`.

## Usage

``` r
dodgr_isochrones(
  graph,
  from = NULL,
  tlim = NULL,
  concavity = 0,
  length_threshold = 0,
  heap = "BHeap"
)
```

## Arguments

- graph:

  `data.frame` or equivalent object representing the network graph. For
  `dodgr` street networks, this must be a network derived from silicate
  ("sc") data, generated with
  [weight_streetnet](https://UrbanAnalyst.github.io/dodgr/reference/weight_streetnet.md).
  This function does not work with networks derived from sf data.

- from:

  Vector or matrix of points **from** which isochrones are to be
  calculated.

- tlim:

  Vector of desired limits of isochrones in seconds

- concavity:

  A value between 0 and 1, with 0 giving (generally smoother but less
  detailed) convex iso-contours and 1 giving highly concave (and
  generally more detailed) contours.

- length_threshold:

  The minimal length of a segment of the iso-contour to be made more
  convex according to the 'concavity\` parameter.. Low values will
  produce highly detailed hulls which may cause problems; if in doubt,
  or if odd results appear, increase this value.

- heap:

  Type of heap to use in priority queue. Options include Fibonacci Heap
  (default; `FHeap`), Binary Heap (`BHeap`), Trinomial Heap (`TriHeap`),
  Extended Trinomial Heap (`TriHeapExt`, and 2-3 Heap (`Heap23`).

## Value

A single `data.frame` of isochrones as points sorted anticlockwise
around each origin (`from`) point, with columns denoting the `from`
points and `tlim` value(s). The isochrones are given as `id` values and
associated coordinates of the series of points from each `from` point at
the specified isochrone times.

Isochrones are calculated by default using parallel computation with the
maximal number of available cores or threads. This number can be reduced
by specifying a value via
`RcppParallel::setThreadOptions (numThreads = <desired_number>)`.

## See also

Other iso:
[`dodgr_isodists()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_isodists.md),
[`dodgr_isoverts()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_isoverts.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Use osmdata package to extract 'SC'-format data:
library (osmdata)
dat <- opq ("hampi india") %>%
    add_osm_feature (key = "highway") %>%
    osmdata_sc ()
graph <- weight_streetnet (dat)
from <- sample (graph$.vx0, size = 100)
tlim <- c (5, 10, 20, 30, 60) * 60 # times in seconds
x <- dodgr_isochrones (graph, from = from, tlim)
} # }
```
