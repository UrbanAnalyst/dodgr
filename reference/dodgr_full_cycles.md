# Calculate fundamental cycles on a FULL (that is, non-contracted) graph.

Calculate fundamental cycles on a FULL (that is, non-contracted) graph.

## Usage

``` r
dodgr_full_cycles(graph, graph_max_size = 10000, expand = 0.05)
```

## Arguments

- graph:

  `data.frame` or equivalent object representing the contracted network
  graph (see Details).

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

This function converts the `graph` to its contracted form, calculates
the fundamental cycles on that version, and then expands these cycles
back onto the original graph. This is far more computationally efficient
than calculating fundamental cycles on a full (non-contracted) graph.

## See also

Other misc:
[`compare_heaps()`](https://UrbanAnalyst.github.io/dodgr/reference/compare_heaps.md),
[`dodgr_flowmap()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_flowmap.md),
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
if (FALSE) { # \dontrun{
net <- weight_streetnet (hampi)
graph <- dodgr_contract_graph (net)
cyc1 <- dodgr_fundamental_cycles (graph)
cyc2 <- dodgr_full_cycles (net)
} # }
# cyc2 has same number of cycles, but each one is generally longer, through
# including all points intermediate to junctions; cyc1 has cycles composed of
# junction points only.
```
