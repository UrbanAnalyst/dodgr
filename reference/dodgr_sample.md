# Sample a random but connected sub-component of a graph

Sample a random but connected sub-component of a graph

## Usage

``` r
dodgr_sample(graph, nverts = 1000)
```

## Arguments

- graph:

  A flat table of graph edges. Must contain columns labelled `from` and
  `to`, or `start` and `stop`. May also contain similarly labelled
  columns of spatial coordinates (for example `from_x`) or `stop_lon`).

- nverts:

  Number of vertices to sample

## Value

A connected sub-component of `graph`

## Note

Graphs may occasionally have `nverts + 1` vertices, rather than the
requested `nverts`.

## See also

Other misc:
[`compare_heaps()`](https://UrbanAnalyst.github.io/dodgr/reference/compare_heaps.md),
[`dodgr_flowmap()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_flowmap.md),
[`dodgr_full_cycles()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_full_cycles.md),
[`dodgr_fundamental_cycles()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_fundamental_cycles.md),
[`dodgr_insert_vertex()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_insert_vertex.md),
[`dodgr_sflines_to_poly()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_sflines_to_poly.md),
[`dodgr_vertices()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_vertices.md),
[`merge_directed_graph()`](https://UrbanAnalyst.github.io/dodgr/reference/merge_directed_graph.md),
[`summary.dodgr_dists_categorical()`](https://UrbanAnalyst.github.io/dodgr/reference/summary.dodgr_dists_categorical.md),
[`write_dodgr_wt_profile()`](https://UrbanAnalyst.github.io/dodgr/reference/write_dodgr_wt_profile.md)

## Examples

``` r
graph <- weight_streetnet (hampi)
nrow (graph) # 5,742
#> [1] 6813
graph <- dodgr_sample (graph, nverts = 200)
nrow (graph) # generally around 400 edges
#> [1] 399
nrow (dodgr_vertices (graph)) # 200
#> [1] 200
```
