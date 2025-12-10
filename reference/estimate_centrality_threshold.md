# Estimate a value for the 'dist_threshold' parameter of the [dodgr_centrality](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_centrality.md) function.

Providing distance thresholds to this function generally provides
considerably speed gains, and results in approximations of centrality.
This function enables the determination of values of 'dist_threshold'
corresponding to specific degrees of accuracy.

## Usage

``` r
estimate_centrality_threshold(graph, tolerance = 0.001)
```

## Arguments

- graph:

  'data.frame' or equivalent object representing the network graph (see
  Details)

- tolerance:

  Desired maximal degree of inaccuracy in centrality estimates

  - values will be accurate to within this amount, subject to a constant
    scaling factor. Note that threshold values increase non-linearly
    with decreasing values of 'tolerance'

## Value

A single value for 'dist_threshold' giving the required tolerance.

## Note

This function may take some time to execute. While running, it displays
ongoing information on screen of estimated values of 'dist_threshold'
and associated errors. Thresholds are progressively increased until the
error is reduced below the specified tolerance.

## See also

Other centrality:
[`dodgr_centrality()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_centrality.md),
[`estimate_centrality_time()`](https://UrbanAnalyst.github.io/dodgr/reference/estimate_centrality_time.md)

## Examples

``` r
# No threshold estimation possible on this small example graph:
graph <- weight_streetnet (hampi, wt_profile = "foot")
estimate_centrality_threshold (graph)
#> dist_threshold approaches size of graph;
#> Recommended value of 'dist_threshold' remains 'NULL',
#> with centrality calculated across entire graph
#> NULL
```
