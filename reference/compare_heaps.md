# Compare timings of different sort heaps for a given input graph.

Perform timing comparison between different kinds of heaps as well as
with equivalent routines from the igraph package. To do this, a random
sub-graph containing a defined number of vertices is first selected.
Alternatively, this random sub-graph can be pre-generated with the
`dodgr_sample` function and passed directly.

## Usage

``` r
compare_heaps(graph, nverts = 100, replications = 2)
```

## Arguments

- graph:

  `data.frame` object representing the network graph (or a sub-sample
  selected with `dodgr_sample`)

- nverts:

  Number of vertices used to generate random sub-graph. If a non-numeric
  value is given, the whole graph will be used.

- replications:

  Number of replications to be used in comparison

## Value

Result of [`bench::mark`](https://bench.r-lib.org/reference/mark.html)
comparison.

## See also

Other misc:
[`dodgr_flowmap()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_flowmap.md),
[`dodgr_full_cycles()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_full_cycles.md),
[`dodgr_fundamental_cycles()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_fundamental_cycles.md),
[`dodgr_insert_vertex()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_insert_vertex.md),
[`dodgr_sample()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_sample.md),
[`dodgr_sflines_to_poly()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_sflines_to_poly.md),
[`dodgr_vertices()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_vertices.md),
[`merge_directed_graph()`](https://UrbanAnalyst.github.io/dodgr/reference/merge_directed_graph.md),
[`summary.dodgr_dists_categorical()`](https://UrbanAnalyst.github.io/dodgr/reference/summary.dodgr_dists_categorical.md),
[`write_dodgr_wt_profile()`](https://UrbanAnalyst.github.io/dodgr/reference/write_dodgr_wt_profile.md)

## Examples

``` r
graph <- weight_streetnet (hampi)
if (FALSE) { # \dontrun{
compare_heaps (graph, nverts = 1000, replications = 1)
} # }
```
