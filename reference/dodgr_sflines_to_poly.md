# Convert sf `LINESTRING` objects to `POLYGON` objects representing all fundamental cycles within the `LINESTRING` objects.

Convert sf `LINESTRING` objects to `POLYGON` objects representing all
fundamental cycles within the `LINESTRING` objects.

## Usage

``` r
dodgr_sflines_to_poly(sflines, graph_max_size = 10000, expand = 0.05)
```

## Arguments

- sflines:

  An sf `LINESTRING` object representing a network.

- graph_max_size:

  Maximum size submitted to the internal C++ routines as a single chunk.
  Warning: Increasing this may lead to computer meltdown!

- expand:

  For large graphs which must be broken into chunks, this factor
  determines the relative overlap between chunks to ensure all cycles
  are captured. (This value should only need to be modified in special
  cases.)

## Value

An [`sf::sfc`](https://r-spatial.github.io/sf/reference/sfc.html)
collection of `POLYGON` objects.

## See also

Other misc:
[`compare_heaps()`](https://UrbanAnalyst.github.io/dodgr/reference/compare_heaps.md),
[`dodgr_flowmap()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_flowmap.md),
[`dodgr_full_cycles()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_full_cycles.md),
[`dodgr_fundamental_cycles()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_fundamental_cycles.md),
[`dodgr_insert_vertex()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_insert_vertex.md),
[`dodgr_sample()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_sample.md),
[`dodgr_vertices()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_vertices.md),
[`merge_directed_graph()`](https://UrbanAnalyst.github.io/dodgr/reference/merge_directed_graph.md),
[`summary.dodgr_dists_categorical()`](https://UrbanAnalyst.github.io/dodgr/reference/summary.dodgr_dists_categorical.md),
[`write_dodgr_wt_profile()`](https://UrbanAnalyst.github.io/dodgr/reference/write_dodgr_wt_profile.md)

## Examples

``` r
if (FALSE) { # \dontrun{
p <- dodgr_sflines_to_poly (hampi)
} # }
```
