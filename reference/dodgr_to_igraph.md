# Convert a `dodgr` graph to an igraph.

Convert a `dodgr` graph to an igraph.

## Usage

``` r
dodgr_to_igraph(graph, weight_column = "d")
```

## Arguments

- graph:

  A `dodgr` graph

- weight_column:

  The column of the `dodgr` network to use as the edge weights in the
  `igraph` representation.

## Value

The `igraph` equivalent of the input. Note that this will *not* be a
dual-weighted graph.

## See also

[igraph_to_dodgr](https://UrbanAnalyst.github.io/dodgr/reference/igraph_to_dodgr.md)

Other conversion:
[`dodgr_deduplicate_graph()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_deduplicate_graph.md),
[`dodgr_to_sf()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_to_sf.md),
[`dodgr_to_sfc()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_to_sfc.md),
[`dodgr_to_tidygraph()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_to_tidygraph.md),
[`igraph_to_dodgr()`](https://UrbanAnalyst.github.io/dodgr/reference/igraph_to_dodgr.md)

## Examples

``` r
graph <- weight_streetnet (hampi)
graphi <- dodgr_to_igraph (graph)
```
