# Optimal allocation of flows from sources to capacity-limited targets

Solve the transportation problem of allocating flows from `N` source
('from') points with associated densities to `M` target ('to') points
with finite capacities, such that total flow from each source is fully
allocated, no target receives more than its capacity, and total
allocation cost (source-to-target network distance) is minimized. The
resultant optimal allocation is then aggregated on to the network with
[dodgr_flows_aggregate](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_flows_aggregate.md),
so this function returns a graph with an additional `flow` column, just
like
[dodgr_flows_aggregate](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_flows_aggregate.md),
[dodgr_flows_disperse](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_flows_disperse.md),
and
[dodgr_flows_si](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_flows_si.md).

## Usage

``` r
dodgr_flows_optalloc(
  graph,
  from,
  to,
  source_densities,
  target_capacities,
  control = list(algorithm = "sinkhorn"),
  contract = TRUE,
  heap = "BHeap",
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

- source_densities:

  Numeric vector of densities at each `from` point, with
  `length(source_densities) == length(from)`.

- target_capacities:

  Numeric vector of capacities at each `to` point, with
  `length(target_capacities) == length(to)`.

- control:

  A named list controlling the allocation algorithm. Must include an
  `algorithm` entry, either `"sinkhorn"` (default) or `"lp"`:

  - `"sinkhorn"` solves an entropic-regularized approximation to the
    optimal allocation via iterative matrix scaling. This is generally
    much faster for large numbers of source/target points, at the cost
    of only approximating the true optimum. Recognizes the following
    additional `control` entries, all optional:

    - `lambda` Entropic regularization strength (default `1`). Smaller
      values approach the exact optimum more closely, at the cost of
      slower and less numerically stable convergence.

    - `tol` Convergence tolerance on row/column marginal errors (default
      `1e-8`).

    - `maxiter` Maximum number of scaling iterations (default `1000`).

  - `"lp"` solves the exact transportation linear program via lpSolve,
    which must be installed (it is only a "Suggested", not "Imported"
    dependency). No additional `control` entries are used.

- contract:

  If `TRUE` (default), calculate flows on contracted graph before
  mapping them back on to the original full graph (recommended as this
  will generally be much faster). `FALSE` should only be used if the
  `graph` has already been contracted.

- heap:

  Type of heap to use in priority queue. Options include Fibonacci Heap
  (default; `FHeap`), Binary Heap (`BHeap`), Trinomial Heap (`TriHeap`),
  Extended Trinomial Heap (`TriHeapExt`, and 2-3 Heap (`Heap23`).

- norm_sums:

  Standardise sums from all origin points, so sum of flows throughout
  entire network equals sum of densities from all origins (see Note).

- quiet:

  If `FALSE`, display progress messages on screen.

## Value

The input `graph`, with an additional `flow` column added, similar to
behaviour of
[dodgr_flows_aggregate](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_flows_aggregate.md).

## Details

This function performs an initial call to
[dodgr_dists](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_dists.md)
to obtain the `N x M` matrix of shortest-path distances between all
`from` and `to` points. The optimal allocation is then obtained by
numerical optimization over that matrix alone (see `control`, below),
with no further path-finding required, before finally calling
[dodgr_flows_aggregate](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_flows_aggregate.md)
with the resultant allocation matrix as its `flows` argument.

Because targets may collectively have more capacity than sources have
density, `sum(source_densities) <= sum(target_capacities)` must hold.
This is checked before any allocation is attempted.

## Note

The `"sinkhorn"` algorithm is generally the faster choice, especially
for large numbers of source/target points, but only approximates the
true optimal allocation, with accuracy controlled by `control$lambda`.
Use `control = list(algorithm = "lp")` for the exact optimum, at the
cost of both requiring the lpSolve package and scaling less well to
large numbers of points.

## See also

Other flows:
[`dodgr_flows_aggregate()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_flows_aggregate.md),
[`dodgr_flows_disperse()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_flows_disperse.md),
[`dodgr_flows_si()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_flows_si.md)

## Examples

``` r
graph <- weight_streetnet (hampi)
graphc <- dodgr_contract_graph (graph)
set.seed (1)
from <- sample (graphc$from_id, size = 10)
to <- sample (graphc$to_id, size = 5)
to <- to [!to %in% from]
source_densities <- runif (length (from))
target_capacities <- runif (length (to))
# scale target_capacities to ensure sum(source_densities) <=
# sum(target_capacities):
target_capacities <- target_capacities *
    1.5 * sum (source_densities) / sum (target_capacities)
graph <- dodgr_flows_optalloc (
    graph,
    from = from,
    to = to,
    source_densities = source_densities,
    target_capacities = target_capacities
)
# graph then has an additional 'flow' column, exactly as for
# 'dodgr_flows_aggregate'

# The exact optimum can be obtained instead with the 'lpSolve' package:
graph <- dodgr_flows_optalloc (
    graph,
    from = from,
    to = to,
    source_densities = source_densities,
    target_capacities = target_capacities,
    control = list (algorithm = "lp")
)
#> Warning: graph already has a 'flow' column; this will be overwritten
```
