# Merge directed edges into equivalent undirected edges.

Merge directed edges into equivalent undirected values by aggregating
across directions. This function is primarily intended to aid
visualisation of directed graphs, particularly visualising the results
of the
[dodgr_flows_aggregate](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_flows_aggregate.md)
and
[dodgr_flows_disperse](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_flows_disperse.md)
functions, which return columns of aggregated flows directed along each
edge of a graph.

## Usage

``` r
merge_directed_graph(graph, col_names = c("flow"))
```

## Arguments

- graph:

  A undirected graph in which directed edges of the input graph have
  been merged through aggregation to yield a single, undirected edge
  between each pair of vertices.

- col_names:

  Names of columns to be merged through aggregation. Values for these
  columns in resultant undirected graph will be aggregated from directed
  values.

## Value

An equivalent graph in which all directed edges have been reduced to
single, undirected edges, and all values of the specified column(s) have
been aggregated across directions to undirected values.

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
[`summary.dodgr_dists_categorical()`](https://UrbanAnalyst.github.io/dodgr/reference/summary.dodgr_dists_categorical.md),
[`write_dodgr_wt_profile()`](https://UrbanAnalyst.github.io/dodgr/reference/write_dodgr_wt_profile.md)

## Examples

``` r
graph <- weight_streetnet (hampi)
from <- sample (graph$from_id, size = 10)
to <- sample (graph$to_id, size = 5)
to <- to [!to %in% from]
flows <- matrix (10 * runif (length (from) * length (to)),
    nrow = length (from)
)
graph <- dodgr_flows_aggregate (graph, from = from, to = to, flows = flows)
# graph then has an additonal 'flows` column of aggregate flows along all
# edges. These flows are directed, and can be aggregated to equivalent
# undirected flows on an equivalent undirected graph with:
graph_undir <- merge_directed_graph (graph)
# This graph will only include those edges having non-zero flows, and so:
nrow (graph)
#> [1] 6813
nrow (graph_undir) # the latter is much smaller
#> [1] 754
```
