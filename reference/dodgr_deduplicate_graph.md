# Deduplicate edges in a graph

Graph may have duplicated edges, particularly when extracted as
[dodgr_streetnet](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_streetnet.md)
objects. This function de-duplicates any repeated edges, reducing
weighted distances and times to the minimal values from all duplicates.

## Usage

``` r
dodgr_deduplicate_graph(graph)
```

## Arguments

- graph:

  Any 'dodgr' graph or network.

## Value

A potentially modified version of graph, with any formerly duplicated
edges reduces to single rows containing minimal weighted distances and
times.

## See also

Other conversion:
[`dodgr_to_igraph()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_to_igraph.md),
[`dodgr_to_sf()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_to_sf.md),
[`dodgr_to_sfc()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_to_sfc.md),
[`dodgr_to_tidygraph()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_to_tidygraph.md),
[`igraph_to_dodgr()`](https://UrbanAnalyst.github.io/dodgr/reference/igraph_to_dodgr.md)

## Examples

``` r
net0 <- weight_streetnet (hampi, wt_profile = "foot")
nrow (net0)
#> [1] 6813
# Duplicate part of input data:
h2 <- rbind (hampi, hampi [1, ])
net1 <- weight_streetnet (h2, wt_profile = "foot")
nrow (net1) # network then has more edges
#> [1] 6871
net2 <- dodgr_deduplicate_graph (net1)
nrow (net2)
#> [1] 6813
stopifnot (identical (nrow (net0), nrow (net2)))
```
