# Transform a result from [dodgr_dists_categorical](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_dists_categorical.md) to summary statistics

Transform a result from
[dodgr_dists_categorical](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_dists_categorical.md)
to summary statistics

## Usage

``` r
# S3 method for class 'dodgr_dists_categorical'
summary(object, ...)
```

## Arguments

- object:

  A 'dodgr_dists_categorical' object

- ...:

  Extra parameters currently not used

## Value

The summary statistics (invisibly)

## See also

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
[`write_dodgr_wt_profile()`](https://UrbanAnalyst.github.io/dodgr/reference/write_dodgr_wt_profile.md)

## Examples

``` r
# Prepare a graph for categorical routing by including an "edge_type" column
graph <- weight_streetnet (hampi, wt_profile = "foot")
graph <- graph [graph$component == 1, ]
graph$edge_type <- graph$highway
# Define start and end points for categorical distances; using all vertices
# here.
length (unique (graph$edge_type)) # Number of categories
#> [1] 8
v <- dodgr_vertices (graph)
from <- to <- v$id [1:100]
d <- dodgr_dists_categorical (graph, from, to)
# Internal 'summary' method to summarise results:
summary (d)
#> Proportional distances along each kind of edge:
#>   path: 0.8941
#>   primary: 0
#>   residential: 0
#>   secondary: 0.0159
#>   service: 0.0572
#>   steps: 0
#>   track: 0
#>   unclassified: 0.0328
```
