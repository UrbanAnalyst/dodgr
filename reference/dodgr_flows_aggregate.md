# Aggregate flows throughout a network.

Aggregate flows throughout a network based on an input matrix of flows
between all pairs of `from` and `to` points. Flows are calculated by
default on contracted graphs, via the `contract = TRUE` parameter.
(These are derived by reducing the input graph down to junction vertices
only, by joining all intermediate edges between each junction.) If
changes to the input graph do not prompt changes to resultant flows, and
the default `contract = TRUE` is used, it may be that calculations are
using previously cached versions of the contracted graph. If so, please
use either
[clear_dodgr_cache](https://UrbanAnalyst.github.io/dodgr/reference/clear_dodgr_cache.md)
to remove the cached version, or
[dodgr_cache_off](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_cache_off.md)
prior to initial graph construction to switch the cache off completely.

## Usage

``` r
dodgr_flows_aggregate(
  graph,
  from,
  to,
  flows,
  pairwise = FALSE,
  contract = TRUE,
  heap = "BHeap",
  tol = 0.000000000001,
  norm_sums = TRUE,
  quiet = TRUE
)
```

## Arguments

- graph:

  `data.frame` or equivalent object representing the network graph (see
  Details)

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

- flows:

  Matrix of flows with `nrow(flows)==length(from)` and
  `ncol(flows)==length(to)`.

- pairwise:

  If `TRUE`, aggregate flows only only paths connecting the ordered
  pairs of `from` and `to`. In this case, both `from` and `to` must be
  of the same length, and `flows` must be either a vector of the same
  length, or a matrix with only one column and same number of rows.
  `flows` then quantifies the flows between each pair of `from` and `to`
  points.

- contract:

  If `TRUE` (default), calculate flows on contracted graph before
  mapping them back on to the original full graph (recommended as this
  will generally be much faster). `FALSE` should only be used if the
  `graph` has already been contracted.

- heap:

  Type of heap to use in priority queue. Options include Fibonacci Heap
  (default; `FHeap`), Binary Heap (`BHeap`), Trinomial Heap (`TriHeap`),
  Extended Trinomial Heap (`TriHeapExt`, and 2-3 Heap (`Heap23`).

- tol:

  Relative tolerance below which flows towards `to` vertices are not
  considered. This will generally have no effect, but can provide speed
  gains when flow matrices represent spatial interaction models, in
  which case this parameter effectively reduces the radius from each
  `from` point over which flows are aggregated. To remove any such
  effect, set `tol = 0`.

- norm_sums:

  Standardise sums from all origin points, so sum of flows throughout
  entire network equals sum of densities from all origins (see Note).

- quiet:

  If `FALSE`, display progress messages on screen.

## Value

Modified version of graph with additional `flow` column added.

## Note

The `norm_sums` parameter should be used whenever densities at origins
and destinations are absolute values, and ensures that the sum of
resultant flow values throughout the entire network equals the sum of
densities at all origins. For example, with `norm_sums = TRUE` (the
default), a flow from a single origin with density one to a single
destination along two edges will allocate flows of one half to each of
those edges, such that the sum of flows across the network will equal
one, or the sum of densities from all origins. The `norm_sums = TRUE`
option is appropriate where densities are relative values, and ensures
that each edge maintains relative proportions. In the above example,
flows along each of two edges would equal one, for a network sum of two,
or greater than the sum of densities.

Flows are calculated by default using parallel computation with the
maximal number of available cores or threads. This number can be reduced
by specifying a value via
`RcppParallel::setThreadOptions (numThreads = <desired_number>)`.

## See also

Other flows:
[`dodgr_flows_disperse()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_flows_disperse.md),
[`dodgr_flows_si()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_flows_si.md)

## Examples

``` r
graph <- weight_streetnet (hampi)
from <- sample (graph$from_id, size = 10)
to <- sample (graph$to_id, size = 5)
to <- to [!to %in% from]
flows <- matrix (10 * runif (length (from) * length (to)),
    nrow = length (from)
)
graph <- dodgr_flows_aggregate (graph, from = from, to = to, flows = flows)
# graph then has an additonal 'flows' column of aggregate flows along all
# edges. These flows are directed, and can be aggregated to equivalent
# undirected flows on an equivalent undirected graph with:
graph_undir <- merge_directed_graph (graph)
# This graph will only include those edges having non-zero flows, and so:
nrow (graph)
#> [1] 6813
nrow (graph_undir) # the latter is much smaller
#> [1] 996

# The following code can be used to convert the resultant graph to an `sf`
# object suitable for plotting
if (FALSE) { # \dontrun{
gsf <- dodgr_to_sf (graph_undir)

# example of plotting with the 'mapview' package
library (mapview)
flow <- gsf$flow / max (gsf$flow)
ncols <- 30
cols <- c ("lawngreen", "red")
colranmp <- colorRampPalette (cols) (ncols) [ceiling (ncols * flow)]
mapview (gsf, color = colranmp, lwd = 10 * flow)
} # }

# An example of flow aggregation across a generic (non-OSM) highway,
# represented as the `routes_fast` object of the \pkg{stplanr} package,
# which is a SpatialLinesDataFrame containing commuter densities along
# components of a street network.
if (FALSE) { # \dontrun{
library (stplanr)
# merge all of the 'routes_fast' lines into a single network
r <- overline (routes_fast, attrib = "length", buff_dist = 1)
r <- sf::st_as_sf (r)
# then extract the start and end points of each of the original 'routes_fast'
# lines and use these for routing with `dodgr`
l <- lapply (routes_fast@lines, function (i) {
    c (
        sp::coordinates (i) [[1]] [1, ],
        tail (sp::coordinates (i) [[1]], 1)
    )
})
l <- do.call (rbind, l)
xy_start <- l [, 1:2]
xy_end <- l [, 3:4]
# Then just specify a generic OD matrix with uniform values of 1:
flows <- matrix (1, nrow = nrow (l), ncol = nrow (l))
# We need to specify both a `type` and `id` column for the
# \link{weight_streetnet} function.
r$type <- 1
r$id <- seq (nrow (r))
graph <- weight_streetnet (
    r,
    type_col = "type",
    id_col = "id",
    wt_profile = 1
)
f <- dodgr_flows_aggregate (
    graph,
    from = xy_start,
    to = xy_end,
    flows = flows
)
# Then merge directed flows and convert to \pkg{sf} for plotting as before:
f <- merge_directed_graph (f)
geoms <- dodgr_to_sfc (f)
gc <- dodgr_contract_graph (f)
gsf <- sf::st_sf (geoms)
gsf$flow <- gc$flow
# sf plot:
plot (gsf ["flow"])
} # }
```
