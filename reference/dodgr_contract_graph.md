# Contract graph to junction vertices only.

Removes redundant (straight-line) vertices from graph, leaving only
junction vertices.

## Usage

``` r
dodgr_contract_graph(graph, verts = NULL, nocache = FALSE)
```

## Arguments

- graph:

  A flat table of graph edges. Must contain columns labelled `from` and
  `to`, or `start` and `stop`. May also contain similarly labelled
  columns of spatial coordinates (for example `from_x`) or `stop_lon`).

- verts:

  Optional list of vertices to be retained as routing points. These must
  match the `from` and `to` columns of `graph`.

- nocache:

  If `FALSE` (default), load cached version of contracted graph if
  previously calculated and cached. If `TRUE`, then re-contract graph
  even if previously calculated version has been stored in cache.

## Value

A contracted version of the original `graph`, containing the same number
of columns, but with each row representing an edge between two junction
vertices (or between the submitted `verts`, which may or may not be
junctions).

## See also

Other modification:
[`dodgr_components()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_components.md),
[`dodgr_uncontract_graph()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_uncontract_graph.md)

## Examples

``` r
graph <- weight_streetnet (hampi)
nrow (graph) # 5,973
#> [1] 6813
graph <- dodgr_contract_graph (graph)
nrow (graph) # 662
#> [1] 744
```
