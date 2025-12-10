# Convert a `dodgr` graph into an equivalent sf object.

Works by aggregating edges into `LINESTRING` objects representing
longest sequences between all junction nodes. The resultant objects will
generally contain more `LINESTRING` objects than the original sf object,
because the former will be bisected at every junction point.

## Usage

``` r
dodgr_to_sf(graph)
```

## Arguments

- graph:

  A `dodgr` graph

## Value

Equivalent object of class sf.

## Note

Requires the sf package to be installed.

## See also

Other conversion:
[`dodgr_deduplicate_graph()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_deduplicate_graph.md),
[`dodgr_to_igraph()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_to_igraph.md),
[`dodgr_to_sfc()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_to_sfc.md),
[`dodgr_to_tidygraph()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_to_tidygraph.md),
[`igraph_to_dodgr()`](https://UrbanAnalyst.github.io/dodgr/reference/igraph_to_dodgr.md)

## Examples

``` r
hw <- weight_streetnet (hampi)
nrow (hw) # 5,729 edges
#> [1] 6813
xy <- dodgr_to_sf (hw)
dim (xy) # 764 edges; 14 attributes
#> [1] 744  17
```
