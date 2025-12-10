# Aggregate flows dispersed from each point in a network.

Disperse flows throughout a network based on a input vectors of origin
points and associated densities. Dispersal is implemented as an
exponential decay, controlled by a parameter, `k`, so that flows decay
with `exp(-d / k)`, where `d` is distance. The algorithm allows for
efficient fitting of multiple dispersal models for different
coefficients to be fitted with a single call. Values of the dispersal
coefficients, `k`, may take one of the following forms:

- A single numeric value (\> 0), with dispersal along all paths
  calculated with that single value. Return object (see below) will then
  have a single additional column named "flow".

- A vector of length equal to the number of `from` points, with
  dispersal from each point then calculated using the corresponding
  value of `k`. Return object has single additional "flow" column.

- A vector of any other length (that is, \> 1 yet different to number of
  `from` points), in which case different dispersal models will be
  fitted for each of the `n` specified values, and the resultant return
  object will have an additional 'n' columns, named 'flow1', 'flow2',
  ... up to 'n'. These columns must be subsequently matched by the user
  back on to the corresponding 'k' values.

- A matrix with number of rows equal to the number of `from` points, and
  any number of columns. Each column will then specify a distinct
  dispersal model, with different values from each row applied to the
  corresponding `from` points. The return value will then be the same as
  the previous version, with an additional `n` columns, "flow1" to
  "flown".

Flows are calculated by default on contracted graphs, via the
`contract = TRUE` parameter. (These are derived by reducing the input
graph down to junction vertices only, by joining all intermediate edges
between each junction.) If changes to the input graph do not prompt
changes to resultant flows, and the default `contract = TRUE` is used,
it may be that calculations are using previously cached versions of the
contracted graph. If so, please use either
[clear_dodgr_cache](https://UrbanAnalyst.github.io/dodgr/reference/clear_dodgr_cache.md)
to remove the cached version, or
[dodgr_cache_off](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_cache_off.md)
prior to initial graph construction to switch the cache off completely.

## Usage

``` r
dodgr_flows_disperse(
  graph,
  from,
  dens,
  k = 500,
  contract = TRUE,
  heap = "BHeap",
  tol = 0.000000000001,
  quiet = TRUE
)
```

## Arguments

- graph:

  `data.frame` or equivalent object representing the network graph (see
  Details)

- from:

  Vector or matrix of points **from** which aggregate dispersed flows
  are to be calculated (see Details)

- dens:

  Vectors of densities corresponding to the `from` points

- k:

  Width coefficient of exponential diffusion function defined as
  `exp(-d/k)`, in units of distance column of `graph` (metres by
  default). Can also be a vector with same length as `from`, giving
  dispersal coefficients from each point. If value of `k<0` is given, a
  standard logistic polynomial will be used.

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

  Relative tolerance below which dispersal is considered to have
  finished. This parameter can generally be ignored; if in doubt, its
  effect can be removed by setting `tol = 0`.

- quiet:

  If `FALSE`, display progress messages on screen.

## Value

Modified version of graph with additional `flow` column added.

## See also

Other flows:
[`dodgr_flows_aggregate()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_flows_aggregate.md),
[`dodgr_flows_si()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_flows_si.md)

## Examples

``` r
# This is generally needed to explore different values of `k` on same graph:
dodgr_cache_off ()

graph <- weight_streetnet (hampi)
from <- sample (graph$from_id, size = 10)
dens <- rep (1, length (from)) # Uniform densities
graph <- dodgr_flows_disperse (graph, from = from, dens = dens)
# graph then has an additonal 'flows` column of aggregate flows along all
# edges. These flows are directed, and can be aggregated to equivalent
# undirected flows on an equivalent undirected graph with:
graph_undir <- merge_directed_graph (graph)

# Remove `flow` column to avoid warning about over-writing values:
graph$flow <- NULL
# One dispersal coefficient for each origin point:
k <- runif (length (from))
graph <- dodgr_flows_disperse (graph, from = from, dens = dens, k = k)
grep ("^flow", names (graph), value = TRUE)
#> [1] "flow"
# single dispersal model; single "flow" column

# Multiple models, muliple dispersal coefficients:
k <- 1:5
graph$flow <- NULL
graph <- dodgr_flows_disperse (graph, from = from, dens = dens, k = k)
grep ("^flow", names (graph), value = TRUE)
#> [1] "flow1" "flow2" "flow3" "flow4" "flow5"
# Rm all flow columns:
graph [grep ("^flow", names (graph), value = TRUE)] <- NULL

# Multiple models with unique coefficient at each origin point:
k <- matrix (runif (length (from) * 5), ncol = 5)
dim (k)
#> [1] 10  5
graph <- dodgr_flows_disperse (graph, from = from, dens = dens, k = k)
grep ("^flow", names (graph), value = TRUE)
#> [1] "flow1" "flow2" "flow3" "flow4" "flow5"
# 5 "flow" columns again, but this time different dispersal coefficients each
# each origin point.
```
