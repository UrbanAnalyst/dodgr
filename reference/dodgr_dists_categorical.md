# Cumulative distances along different edge categories

The main
[dodgr_distances](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_distances.md)
function calculates distances between input sets of origin and
destination points, and returns a single matrix of numeric distances.
This function aggregates distances along categories of edges or
segments, and returns an overall distance matrix (identical to the
result of
[dodgr_distances](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_distances.md)),
along with one additional matrix for each edge category.

Edges types must be specified in a column of the input graph named
"edge_type". If this has two types of values (for example, "a" and "b"),
then the function will return two additional distance matrices, one of
total lengths of distances between all pairs of points traversed along
edges of the first type, "a", and one of aggregated distances along
edges of type "b".

See the description of
[dodgr_distances](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_distances.md)
for details on the expected format of input graphs.

## Usage

``` r
dodgr_dists_categorical(
  graph,
  from = NULL,
  to = NULL,
  proportions_only = FALSE,
  pairwise = FALSE,
  dlimit = NULL,
  heap = "BHeap",
  quiet = TRUE
)
```

## Arguments

- graph:

  `data.frame` or equivalent object representing the network graph which
  must have a column named "edge_type" which labels categories of edge
  types along which categorical distances are to be aggregated (see
  Note).

- from:

  Vector or matrix of points **from** which route distances are to be
  calculated, specified as one of the following:

  - Single character vector precisely matching node numbers or names
    given in `graph$from` or `graph$to`.

  - Single vector of integer-ish values, in which case these will be
    presumed to specify indices into
    [dodgr_vertices](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_vertices.md),
    and NOT to correspond to values in the 'from' or 'to' columns of the
    graph. See the example below for a demonstration.

  - Matrix or equivalent of longitude and latitude coordinates, in which
    case these will be matched on to the nearest coordinates of 'from'
    and 'to' points in the graph.

- to:

  Vector or matrix of points **to** which route distances are to be
  calculated. If `to` is `NULL`, pairwise distances will be calculated
  from all `from` points to all other nodes in `graph`. If both `from`
  and `to` are `NULL`, pairwise distances are calculated between all
  nodes in `graph`.

- proportions_only:

  If `FALSE`, return distance matrices for full distances and for each
  edge category; if `TRUE`, return single vector of proportional
  distances, like the `summary` function applied to full results. See
  Note.

- pairwise:

  If `TRUE`, calculate distances only between the ordered pairs of
  `from` and `to`. In this case, neither the `proportions_only` nor
  `dlimit` parameters have any effect, and the result is a single matrix
  with one row for each pair of `from`-`to` points, and one column for
  each category.

- dlimit:

  If no value to `to` is given, distances are aggregated from each
  `from` point out to the specified distance limit (in the same units as
  the edge distances of the input graph). `dlimit` only has any effect
  if `to` is not specified, in which case the `proportions_only`
  argument has no effect.

- heap:

  Type of heap to use in priority queue. Options include Fibonacci Heap
  (default; `FHeap`), Binary Heap (`BHeap`),
  `Trinomial Heap (`TriHeap`), Extended Trinomial Heap (`TriHeapExt`, and 2-3 Heap (`Heap23\`).

- quiet:

  If `FALSE`, display progress messages on screen.

## Value

If `to` is specified, a list of distance matrices of equal dimensions
(length(from), length(to)), the first of which ("distance") holds the
final distances, while the rest are one matrix for each unique value of
"edge_type", holding the distances traversed along those types of edges
only. Otherwise, a single matrix of total distances along all ways from
each point out to the specified value of `dlimit`, along with distances
along each of the different kinds of ways specified in the "edge_type"
column of the input graph.

## Note

The "edge_type" column in the graph can contain any kind of discrete or
categorical values, although integer values of 0 are not permissible.
`NA` values are ignored. The function requires one full distance matrix
to be stored for each category of "edge_type" (unless
`proportions_only = TRUE`). It is wise to keep numbers of discrete types
as low as possible, especially for large distance matrices.

Setting the `proportions_only` flag to `TRUE` may be advantageous for
large jobs, because this avoids construction of the full matrices. This
may speed up calculations, but perhaps more importantly it may make
possible calculations which would otherwise require distance matrices
too large to be directly stored.

Calculations are not able to be interrupted (for example, by `Ctrl-C`),
and can only be stopped by killing the R process.

## See also

Other distances:
[`dodgr_distances()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_distances.md),
[`dodgr_dists()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_dists.md),
[`dodgr_dists_nearest()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_dists_nearest.md),
[`dodgr_paths()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_paths.md),
[`dodgr_times()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_times.md)

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
#>   service: 0.0571
#>   steps: 0
#>   track: 0
#>   unclassified: 0.0329

class (d)
#> [1] "list"                    "dodgr_dists_categorical"
length (d)
#> [1] 9
sapply (d, dim)
#>      distances path primary residential secondary service steps track
#> [1,]       100  100     100         100       100     100   100   100
#> [2,]       100  100     100         100       100     100   100   100
#>      unclassified
#> [1,]          100
#> [2,]          100
# 9 distance matrices, all of same dimensions, first of which is standard
# distance matrix
s <- summary (d) # print summary as proportions along each "edge_type"
#> Proportional distances along each kind of edge:
#>   path: 0.8941
#>   primary: 0
#>   residential: 0
#>   secondary: 0.0159
#>   service: 0.0571
#>   steps: 0
#>   track: 0
#>   unclassified: 0.0329
# or directly calculate proportions only
dodgr_dists_categorical (graph, from, to,
    proportions_only = TRUE
)
#>         path      primary  residential    secondary      service        steps 
#>   0.89410618   0.00000000   0.00000000   0.01594045   0.05706168   0.00000000 
#>        track unclassified 
#>   0.00000000   0.03289169 

# Pairwise distances return single matrix with number of rows equal to 'from'
# / 'to', and number of columns equal to number of edge types plus one for
# total distances.
d <- dodgr_dists_categorical (graph, from, to, pairwise = TRUE)
class (d)
#> [1] "matrix" "array" 
dim (d)
#> [1] 100   9

# The 'dlimit' parameter can be used to calculate total distances along each
# category of edges from a set of points out to specified threshold:
dlimit <- 2000 # in metres
d <- dodgr_dists_categorical (graph, from, dlimit = dlimit)
dim (d) # length(from), length(unique(edge_type)) + 1
#> [1] 100   9
```
