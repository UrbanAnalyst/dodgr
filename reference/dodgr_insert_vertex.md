# Insert a new node or vertex into a network

Insert a new node or vertex into a network

## Usage

``` r
dodgr_insert_vertex(graph, v1, v2, x = NULL, y = NULL)
```

## Arguments

- graph:

  A flat table of graph edges. Must contain columns labelled `from` and
  `to`, or `start` and `stop`. May also contain similarly labelled
  columns of spatial coordinates (for example `from_x`) or `stop_lon`).

- v1:

  Vertex defining start of graph edge along which new vertex is to be
  inserted

- v2:

  Vertex defining end of graph edge along which new vertex is to be
  inserted (order of `v1` and `v2` is not important).

- x:

  The `x`-coordinate of new vertex. If not specified, vertex is created
  half-way between `v1` and `v2`.

- y:

  The `y`-coordinate of new vertex. If not specified, vertex is created
  half-way between `v1` and `v2`.

## Value

A modified graph with specified edge between defined start and end
vertices split into two edges either side of new vertex.

## See also

Other misc:
[`compare_heaps()`](https://UrbanAnalyst.github.io/dodgr/reference/compare_heaps.md),
[`dodgr_flowmap()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_flowmap.md),
[`dodgr_full_cycles()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_full_cycles.md),
[`dodgr_fundamental_cycles()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_fundamental_cycles.md),
[`dodgr_sample()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_sample.md),
[`dodgr_sflines_to_poly()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_sflines_to_poly.md),
[`dodgr_vertices()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_vertices.md),
[`merge_directed_graph()`](https://UrbanAnalyst.github.io/dodgr/reference/merge_directed_graph.md),
[`summary.dodgr_dists_categorical()`](https://UrbanAnalyst.github.io/dodgr/reference/summary.dodgr_dists_categorical.md),
[`write_dodgr_wt_profile()`](https://UrbanAnalyst.github.io/dodgr/reference/write_dodgr_wt_profile.md)

## Examples

``` r
graph <- weight_streetnet (hampi)
e1 <- sample (nrow (graph), 1)
v1 <- graph$from_id [e1]
v2 <- graph$to_id [e1]
# insert new vertex in the middle of that randomly-selected edge:
graph2 <- dodgr_insert_vertex (graph, v1, v2)
nrow (graph)
#> [1] 6813
nrow (graph2) # new edges added to graph2
#> [1] 6815
```
