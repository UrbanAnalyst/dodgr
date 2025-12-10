# Calculate fundamental cycles in a graph.

Calculate fundamental cycles in a graph.

## Usage

``` r
dodgr_fundamental_cycles(
  graph,
  vertices = NULL,
  graph_max_size = 10000,
  expand = 0.05
)
```

## Arguments

- graph:

  `data.frame` or equivalent object representing the contracted network
  graph (see Details).

- vertices:

  `data.frame` returned from
  [dodgr_vertices](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_vertices.md)`(graph)`.
  Will be calculated if not provided, but it's quicker to pass this if
  it has already been calculated.

- graph_max_size:

  Maximum size submitted to the internal C++ routines as a single chunk.
  Warning: Increasing this may lead to computer meltdown!

- expand:

  For large graphs which must be broken into chunks, this factor
  determines the relative overlap between chunks to ensure all cycles
  are captured. (This value should only need to be modified in special
  cases.)

## Value

List of cycle paths, in terms of vertex IDs in `graph` and, for spatial
graphs, the corresponding coordinates.

## Note

Calculation of fundamental cycles is VERY computationally demanding, and
this function should only be executed on CONTRACTED graphs (that is,
graphs returned from
[dodgr_contract_graph](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_contract_graph.md)),
and even than may take a long time to execute. Results for full graphs
can be obtained with the function
[dodgr_full_cycles](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_full_cycles.md).
The computational complexity can also not be calculated in advance, and
so the parameter `graph_max_size` will lead to graphs larger than that
(measured in numbers of edges) being cut into smaller parts. (Note that
that is only possible for spatial graphs, meaning that it is not at all
possible to apply this function to large, non-spatial graphs.) Each of
these smaller parts will be expanded by the specified amount (`expand`),
and cycles found within. The final result is obtained by aggregating all
of these cycles and removing any repeated ones arising due to overlap in
the expanded portions. Finally, note that this procedure of cutting
graphs into smaller, computationally manageable sub-graphs provides only
an approximation and may not yield all fundamental cycles.

## See also

Other misc:
[`compare_heaps()`](https://UrbanAnalyst.github.io/dodgr/reference/compare_heaps.md),
[`dodgr_flowmap()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_flowmap.md),
[`dodgr_full_cycles()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_full_cycles.md),
[`dodgr_insert_vertex()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_insert_vertex.md),
[`dodgr_sample()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_sample.md),
[`dodgr_sflines_to_poly()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_sflines_to_poly.md),
[`dodgr_vertices()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_vertices.md),
[`merge_directed_graph()`](https://UrbanAnalyst.github.io/dodgr/reference/merge_directed_graph.md),
[`summary.dodgr_dists_categorical()`](https://UrbanAnalyst.github.io/dodgr/reference/summary.dodgr_dists_categorical.md),
[`write_dodgr_wt_profile()`](https://UrbanAnalyst.github.io/dodgr/reference/write_dodgr_wt_profile.md)

## Examples

``` r
net <- weight_streetnet (hampi)
graph <- dodgr_contract_graph (net)
verts <- dodgr_vertices (graph)
cyc <- dodgr_fundamental_cycles (graph, verts)
```
