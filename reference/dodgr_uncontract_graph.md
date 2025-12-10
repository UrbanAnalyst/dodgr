# Re-expand a contracted graph.

Revert a contracted graph created with
[dodgr_contract_graph](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_contract_graph.md)
back to a full, uncontracted version. This function is mostly used for
the side effect of mapping any new columns inserted on to the contracted
graph back on to the original graph, as demonstrated in the example.

## Usage

``` r
dodgr_uncontract_graph(graph)
```

## Arguments

- graph:

  A contracted graph created from
  [dodgr_contract_graph](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_contract_graph.md).

## Value

A single `data.frame` representing the equivalent original, uncontracted
graph.

## Details

Note that this function will generally *not* recover original graphs
submitted to
[dodgr_contract_graph](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_contract_graph.md).
Specifically, the sequence
`dodgr_contract_graph(graph) |> dodgr_uncontract_graph()` will generally
produce a graph with fewer edges than the original. This is because
graphs may have multiple paths between a given pair of points.
Contraction will reduce these to the single path with the shortest
weighted distance (or time), and uncontraction will restore only that
single edge with shortest weighted distance, and not any original edges
which may have had longer weighted distances.

## See also

Other modification:
[`dodgr_components()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_components.md),
[`dodgr_contract_graph()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_contract_graph.md)

## Examples

``` r
graph0 <- weight_streetnet (hampi)
nrow (graph0) # 6,813
#> [1] 6813
graph1 <- dodgr_contract_graph (graph0)
nrow (graph1) # 760
#> [1] 744
graph2 <- dodgr_uncontract_graph (graph1)
nrow (graph2) # 6,813
#> [1] 6813

# Insert new data on to the contracted graph and uncontract it:
graph1$new_col <- runif (nrow (graph1))
graph3 <- dodgr_uncontract_graph (graph1)
# graph3 is then the uncontracted graph which includes "new_col" as well
dim (graph0)
#> [1] 6813   16
dim (graph3)
#> [1] 6813   17
```
