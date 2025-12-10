# Distances On Directed GRaphs ("dodgr")

Distances on dual-weighted directed graphs using priority-queue shortest
paths. Weighted directed graphs have weights from A to B which may
differ from those from B to A. Dual-weighted directed graphs have two
sets of such weights. A canonical example is a street network to be used
for routing in which routes are calculated by weighting distances
according to the type of way and mode of transport, yet lengths of
routes must be calculated from direct distances.

## The Main Function

- [`dodgr_dists()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_dists.md):
  Calculate pair-wise distances between specified pairs of points in a
  graph.

## Functions to Obtain Graphs

- [`dodgr_streetnet()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_streetnet.md):
  Extract a street network in Simple Features (`sf`) form.

- [`weight_streetnet()`](https://UrbanAnalyst.github.io/dodgr/reference/weight_streetnet.md):
  Convert an `sf`-formatted street network to a `dodgr` graph through
  applying specified weights to all edges.

## Functions to Modify Graphs

- [`dodgr_components()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_components.md):
  Number all graph edges according to their presence in distinct
  connected components.

- [`dodgr_contract_graph()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_contract_graph.md):
  Contract a graph by removing redundant edges.

## Miscellaneous Functions

- [`dodgr_sample()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_sample.md):
  Randomly sample a graph, returning a single connected component of a
  defined number of vertices.

- [`dodgr_vertices()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_vertices.md):
  Extract all vertices of a graph.

- [`compare_heaps()`](https://UrbanAnalyst.github.io/dodgr/reference/compare_heaps.md):
  Compare the performance of different priority queue heap structures
  for a given type of graph.

## See also

Useful links:

- <https://UrbanAnalyst.github.io/dodgr/>

- <https://github.com/UrbanAnalyst/dodgr>

- Report bugs at <https://github.com/UrbanAnalyst/dodgr/issues>

## Author

**Maintainer**: Mark Padgham <mark.padgham@email.com>

Authors:

- Andreas Petutschnig

- David Cooley

Other contributors:

- Robin Lovelace \[contributor\]

- Andrew Smith \[contributor\]

- Malcolm Morgan \[contributor\]

- Andrea Gilardi ([ORCID](https://orcid.org/0000-0002-9424-7439))
  \[contributor\]

- Eduardo Leoni ([ORCID](https://orcid.org/0000-0003-0955-5232))
  \[contributor\]

- Shane Saunders (Original author of included code for priority heaps)
  \[copyright holder\]

- Stanislaw Adaszewski (author of include concaveman-cpp code)
  \[copyright holder\]
