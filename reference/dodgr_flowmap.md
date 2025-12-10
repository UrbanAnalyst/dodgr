# Create a map of `dodgr` flows.

Create a map of the output of
[dodgr_flows_aggregate](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_flows_aggregate.md)
or
[dodgr_flows_disperse](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_flows_disperse.md)

## Usage

``` r
dodgr_flowmap(net, bbox = NULL, linescale = 1)
```

## Arguments

- net:

  A street network with a `flow` column obtained from
  [dodgr_flows_aggregate](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_flows_aggregate.md)
  or
  [dodgr_flows_disperse](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_flows_disperse.md)

- bbox:

  If given, scale the map to this bbox, otherwise use entire extend of
  `net`

- linescale:

  Maximal thickness of plotted lines

## Value

Nothing; called for side-effect of producing plot.

## Note

`net` should be first passed through `merge_directed_graph` prior to
plotting, otherwise lines for different directions will be overlaid.

## See also

Other misc:
[`compare_heaps()`](https://UrbanAnalyst.github.io/dodgr/reference/compare_heaps.md),
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
from <- sample (graph$from_id, size = 10)
to <- sample (graph$to_id, size = 5)
to <- to [!to %in% from]
flows <- matrix (
    10 * runif (length (from) * length (to)),
    nrow = length (from)
)
graph <- dodgr_flows_aggregate (graph, from = from, to = to, flows = flows)
# graph then has an additonal 'flows` column of aggregate flows along all
# edges. These flows are directed, and can be aggregated to equivalent
# undirected flows on an equivalent undirected graph with:
graph_undir <- merge_directed_graph (graph)
if (FALSE) { # \dontrun{
dodgr_flowmap (graph_undir)
} # }
```
