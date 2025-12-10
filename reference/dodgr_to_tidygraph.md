# Convert a `dodgr` graph to an tidygraph.

Convert a `dodgr` graph to an tidygraph.

## Usage

``` r
dodgr_to_tidygraph(graph)
```

## Arguments

- graph:

  A `dodgr` graph

## Value

The `tidygraph` equivalent of the input

## See also

Other conversion:
[`dodgr_deduplicate_graph()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_deduplicate_graph.md),
[`dodgr_to_igraph()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_to_igraph.md),
[`dodgr_to_sf()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_to_sf.md),
[`dodgr_to_sfc()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_to_sfc.md),
[`igraph_to_dodgr()`](https://UrbanAnalyst.github.io/dodgr/reference/igraph_to_dodgr.md)

## Examples

``` r
graph <- weight_streetnet (hampi)
grapht <- dodgr_to_tidygraph (graph)
#> Loading required namespace: tidygraph
```
