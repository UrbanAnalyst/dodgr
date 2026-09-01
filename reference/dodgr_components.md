# Identify connected components of graph.

Identify connected components of graph and add corresponding `component`
column to `data.frame`.

## Usage

``` r
dodgr_components(graph, strong = FALSE)
```

## Arguments

- graph:

  A `data.frame` of edges

- strong:

  Defaults to
  `FALSE', which may identify components which can only be accessed from a single direction, and therefore not actually used in routing calculations. If `TRUE\`,
  all edges in each identified component are fully connected in both
  directions.

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
