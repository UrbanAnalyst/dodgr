# Estimate time required for a planned centrality calculation.

The 'dodgr' centrality functions are designed to be applied to
potentially very large graphs, and may take considerable time to
execute. This helper function estimates how long a centrality function
may take for a given graph and given value of 'dist_threshold' estimated
via the
[estimate_centrality_threshold](https://UrbanAnalyst.github.io/dodgr/reference/estimate_centrality_threshold.md)
function.

## Usage

``` r
estimate_centrality_time(
  graph,
  contract = TRUE,
  edges = TRUE,
  dist_threshold = NULL,
  heap = "BHeap"
)
```

## Arguments

- graph:

  'data.frame' or equivalent object representing the network graph (see
  Details)

- contract:

  If 'TRUE', centrality is calculated on contracted graph before mapping
  back on to the original full graph. Note that for street networks, in
  particular those obtained from the osmdata package, vertex placement
  is effectively arbitrary except at junctions; centrality for such
  graphs should only be calculated between the latter points, and thus
  'contract' should always be 'TRUE'.

- edges:

  If 'TRUE', centrality is calculated for graph edges, returning the
  input 'graph' with an additional 'centrality' column; otherwise
  centrality is calculated for vertices, returning the equivalent of
  'dodgr_vertices(graph)', with an additional vertex-based 'centrality'
  column.

- dist_threshold:

  If not 'NULL', only calculate centrality for each point out to
  specified threshold. Setting values for this will result in
  approximate estimates for centrality, yet with considerable gains in
  computational efficiency. For sufficiently large values,
  approximations will be accurate to within some constant multiplier.
  Appropriate values can be established via the
  [estimate_centrality_threshold](https://UrbanAnalyst.github.io/dodgr/reference/estimate_centrality_threshold.md)
  function.

- heap:

  Type of heap to use in priority queue. Options include Fibonacci Heap
  (default; 'FHeap'), Binary Heap ('BHeap'), Trinomial Heap ('TriHeap'),
  Extended Trinomial Heap ('TriHeapExt', and 2-3 Heap ('Heap23').

## Value

An estimated calculation time for calculating centrality for the given
value of 'dist_threshold'

## Note

This function may take some time to execute. While running, it displays
ongoing information on screen of estimated values of 'dist_threshold'
and associated errors. Thresholds are progressively increased until the
error is reduced below the specified tolerance.

## See also

Other centrality:
[`dodgr_centrality()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_centrality.md),
[`estimate_centrality_threshold()`](https://UrbanAnalyst.github.io/dodgr/reference/estimate_centrality_threshold.md)

## Examples

``` r
graph <- weight_streetnet (hampi, wt_profile = "foot")
estimate_centrality_time (graph)
#> Estimated time to calculate centrality for full graph is 00:00:00
```
