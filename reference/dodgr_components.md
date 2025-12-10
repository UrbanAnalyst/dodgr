# Identify connected components of graph.

Identify connected components of graph and add corresponding `component`
column to `data.frame`.

## Usage

``` r
dodgr_components(graph)
```

## Arguments

- graph:

  A `data.frame` of edges

## Value

Equivalent graph with additional `component` column, sequentially
numbered from 1 = largest component.

## See also

Other modification:
[`dodgr_contract_graph()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_contract_graph.md),
[`dodgr_uncontract_graph()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_uncontract_graph.md)

## Examples

``` r
graph <- weight_streetnet (hampi)
graph <- dodgr_components (graph)
#> graph already has a component column
```
