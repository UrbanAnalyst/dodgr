# Write `dodgr` weighting profiles to local file.

Write the `dodgr` street network weighting profiles to a local
`.json`-formatted file for manual editing and subsequent re-reading.

## Usage

``` r
write_dodgr_wt_profile(file = NULL)
```

## Arguments

- file:

  Full name (including path) of file to which to write. The `.json`
  suffix will be automatically appended.

## Value

TRUE if writing successful.

## See also

[weight_streetnet](https://UrbanAnalyst.github.io/dodgr/reference/weight_streetnet.md)

Other misc:
[`compare_heaps()`](https://UrbanAnalyst.github.io/dodgr/reference/compare_heaps.md),
[`dodgr_flowmap()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_flowmap.md),
[`dodgr_full_cycles()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_full_cycles.md),
[`dodgr_fundamental_cycles()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_fundamental_cycles.md),
[`dodgr_insert_vertex()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_insert_vertex.md),
[`dodgr_sample()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_sample.md),
[`dodgr_sflines_to_poly()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_sflines_to_poly.md),
[`dodgr_vertices()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_vertices.md),
[`merge_directed_graph()`](https://UrbanAnalyst.github.io/dodgr/reference/merge_directed_graph.md),
[`summary.dodgr_dists_categorical()`](https://UrbanAnalyst.github.io/dodgr/reference/summary.dodgr_dists_categorical.md)

## Examples

``` r
f <- tempfile (fileext = ".json")
write_dodgr_wt_profile (file = f)
wt_profiles <- jsonlite::read_json (f, simplify = TRUE)
```
