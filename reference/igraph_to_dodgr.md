# Convert a igraph network to an equivalent `dodgr` representation.

Convert a igraph network to an equivalent `dodgr` representation.

## Usage

``` r
igraph_to_dodgr(graph)
```

## Arguments

- graph:

  An igraph network

## Value

The `dodgr` equivalent of the input.

## See also

[dodgr_to_igraph](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_to_igraph.md)

Other conversion:
[`dodgr_deduplicate_graph()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_deduplicate_graph.md),
[`dodgr_to_igraph()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_to_igraph.md),
[`dodgr_to_sf()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_to_sf.md),
[`dodgr_to_sfc()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_to_sfc.md),
[`dodgr_to_tidygraph()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_to_tidygraph.md)

## Examples

``` r
graph <- weight_streetnet (hampi)
graphi <- dodgr_to_igraph (graph)
graph2 <- igraph_to_dodgr (graphi)
identical (graph2, graph) # FALSE
#> [1] FALSE
```
