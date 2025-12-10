# Calculate isodistance or isochrone vertices from specified points.

Returns lists of all network vertices contained within isodistance or
isochrone contours. Input objects must be `data.frame` objects
(`graph`), which must minimally contain three columns of `from`, `to`,
and `d` or `dist`. If an additional column named `weight` or `wt` is
present, iso contours are evaluate via shortest paths calculated
according to values specified in that column, while resultant values of
iso contours are calculated from the `d` or `dist` column. That is, the
paths tracing iso contours from any point will be calculated according
to the minimal total sum of `weight` values (if present), while reported
iso contours will be total sums of `dist` values.

Graphs derived from Open Street Map street networks, via the
[weight_streetnet](https://UrbanAnalyst.github.io/dodgr/reference/weight_streetnet.md)
function, have columns labelled `d`, `d_weighted`, `time`, and
`time_weighted`. For these inputs, iso contours are always routed using
`d_weighted` (or `t_weighted` for times), while final iso contours
reflect sums of values of `d` (or `t` for times) - that is, of
un-weighted distances or times - along those paths.

Function is fully vectorized to accept vectors of central points and
vectors defining multiple isochrone or isodistance thresholds. Provide
one or more `dlim` values for isodistances, or one or more `tlim` values
for isochrones. Calculations use by default parallel computation with
the maximal number of available cores or threads. This number can be
reduced by specifying a value via
`RcppParallel::setThreadOptions (numThreads = <desired_number>)`.

## Usage

``` r
dodgr_isoverts(graph, from = NULL, dlim = NULL, tlim = NULL, heap = "BHeap")
```

## Arguments

- graph:

  `data.frame` or equivalent object representing the network graph. For
  `dodgr` street networks, this must be a network derived from silicate
  ("sc") data, generated with
  [weight_streetnet](https://UrbanAnalyst.github.io/dodgr/reference/weight_streetnet.md).
  This function does not work with networks derived from sf data.

- from:

  Vector or matrix of points **from** which isodistances or isochrones
  are to be calculated.

- dlim:

  Vector of desired limits of isodistances in metres.

- tlim:

  Vector of desired limits of isochrones in seconds

- heap:

  Type of heap to use in priority queue. Options include Fibonacci Heap
  (default; `FHeap`), Binary Heap (`BHeap`), Trinomial Heap (`TriHeap`),
  Extended Trinomial Heap (`TriHeapExt`, and 2-3 Heap (`Heap23`).

## Value

A single `data.frame` of vertex IDs, with columns denoting the `from`
points and `tlim` value(s). The isochrones are given as `id` values and
associated coordinates of the series of points from each `from` point at
the specified isochrone times.

Isoverts are calculated by default using parallel computation with the
maximal number of available cores or threads. This number can be reduced
by specifying a value via
`RcppParallel::setThreadOptions (numThreads = <desired_number>)`.

## See also

Other iso:
[`dodgr_isochrones()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_isochrones.md),
[`dodgr_isodists()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_isodists.md)

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
x <- dodgr_isoverts (graph, from = from, tlim)
} # }
```
