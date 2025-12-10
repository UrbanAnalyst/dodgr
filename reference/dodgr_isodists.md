# Calculate isodistance contours from specified points.

Calculates isodistances from input `data.frame` objects (`graph`), which
must minimally contain three columns of `from`, `to`, and `d` or `dist`.
If an additional column named `weight` or `wt` is present, shortest
paths are calculated according to values specified in that column, while
resultant isodistances are calculated from the `d` or `dist` column.
That is, the paths tracing isodistances from any point will be
calculated according to the minimal total sum of `weight` values (if
present), while reported isodistances will be total sums of `dist`
values.

Graphs derived from Open Street Map street networks, via the
[weight_streetnet](https://UrbanAnalyst.github.io/dodgr/reference/weight_streetnet.md)
function, have columns labelled `d`, `d_weighted`, `time`, and
`time_weighted`. For these inputs, isodistances are always routed using
`d_weighted` (or `t_weighted` for times), while final isodistances are
sums of values of `d` (or `t` for times)- that is, of un-weighted
distances or times - along those paths.

Function is fully vectorized to calculate accept vectors of central
points and vectors defining multiple isodistances. Calculations use by
default parallel computation with the maximal number of available cores
or threads. This number can be reduced by specifying a value via
`RcppParallel::setThreadOptions (numThreads = <desired_number>)`.

## Usage

``` r
dodgr_isodists(
  graph,
  from = NULL,
  dlim = NULL,
  concavity = 0,
  length_threshold = 0,
  contract = TRUE,
  heap = "BHeap"
)
```

## Arguments

- graph:

  `data.frame` or equivalent object representing the network graph. For
  `dodgr` street networks, this may be a network derived from either sf
  or silicate ("sc") data, generated with
  [weight_streetnet](https://UrbanAnalyst.github.io/dodgr/reference/weight_streetnet.md).

- from:

  Vector or matrix of points **from** which isodistances are to be
  calculated.

- dlim:

  Vector of desired limits of isodistances in metres.

- concavity:

  A value between 0 and 1, with 0 giving (generally smoother but less
  detailed) convex iso-contours and 1 giving highly concave (and
  generally more detailed) contours.

- length_threshold:

  The minimal length of a segment of the iso-contour to be made more
  convex according to the 'concavity\` parameter.. Low values will
  produce highly detailed hulls which may cause problems; if in doubt,
  or if odd results appear, increase this value.

- contract:

  If `TRUE`, calculate isodists only to vertices in the contract graph,
  in other words, only to junction vertices.

- heap:

  Type of heap to use in priority queue. Options include Fibonacci Heap
  (default; `FHeap`), Binary Heap (`BHeap`), Trinomial Heap (`TriHeap`),
  Extended Trinomial Heap (`TriHeapExt`, and 2-3 Heap (`Heap23`).

## Value

A single `data.frame` of isodistances as points sorted anticlockwise
around each origin (`from`) point, with columns denoting the `from`
points and `dlim` value(s). The isodistance contours are given as `id`
values and associated coordinates of the series of points from each
`from` point at the specified isodistances.

## See also

Other iso:
[`dodgr_isochrones()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_isochrones.md),
[`dodgr_isoverts()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_isoverts.md)

## Examples

``` r
graph <- weight_streetnet (hampi)
from <- sample (graph$from_id, size = 100)
dlim <- c (1, 2, 5, 10, 20) * 100
d <- dodgr_isodists (graph, from = from, dlim)
```
