<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Build
Status](https://travis-ci.org/ATFutures/dodgr.svg)](https://travis-ci.org/ATFutures/dodgr)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/ATFutures/dodgr?branch=master&svg=true)](https://ci.appveyor.com/project/ATFutures/dodgr)
[![codecov](https://codecov.io/gh/ATFutures/dodgr/branch/master/graph/badge.svg)](https://codecov.io/gh/ATFutures/dodgr)
[![Project Status:
Active](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/dodgr)](https://cran.r-project.org/package=dodgr)
[![CRAN
Downloads](http://cranlogs.r-pkg.org/badges/grand-total/dodgr?color=orange)](https://cran.r-project.org/package=dodgr)
[![CII Best
Practices](https://bestpractices.coreinfrastructure.org/projects/1396/badge)](https://bestpractices.coreinfrastructure.org/projects/1396)

# dodgr: Distances on Directed Graphs in R

`dodgr` is an R package for efficient calculation of many-to-many
pairwise distances on dual-weighted directed graphs, for aggregation of
flows throughout networks, and for highly realistic routing through
street networks (time-based routing considering incline, turn-angles,
surface quality, everything).

Note that most `dodgr` algorithms implement parallel computation with
the [`RcppParallel` library](https://rcppcore.github.io/RcppParallel/),
and by default use the maximal number of available cores or threads. If
you do not wish `dodgr`to use all available threads, please reduce the
number manually by first specifying a value via

``` r
RcppParallel::setThreadOptions (numThreads = <desired_number>)
```

## What’s so special?

Four aspects. First, while other packages exist for calculating
distances on directed graphs, notably [`igraph`](https://igraph.org/r),
even that otherwise fabulous package does not (readily) permit analysis
of *dual-weighted* graphs. Dual-weighted graphs have two sets of weights
for each edge, so routing can be evaluated with one set of weights,
while distances can be calculated with the other. A canonical example is
a street network, where *weighted distances* are assigned depending on
mode of transport (for example, weighted distances for pedestrians on
multi-lane vehicular roads are longer than equivalent distances along
isolated walking paths), yet the desired output remains direct,
unweighted distances. Accurate calculation of distances on street
networks requires a dual-weighted representation. In **R**, `dodgr` is
currently the only package that offers this functionality (without
excessive data wrangling).

Second, while [`igraph`](https://igraph.org/r) and almost all other
routing packages are primarily designed for one-to-one routing, `dodgr`
is specifically designed for many-to-many routing, and will generally
outperform equivalent packages in large routing tasks.

Third, `dodgr` goes beyond the functionality of comparable packages
through including routines to aggregate flows throughout a network,
through specifying origins, destinations, and flow densities between
each pair of points. Alternatively, flows can be aggregated according to
a network dispersal model from a set of origin points and associated
densities, and a user-specified dispersal model.

Fourth and finally, `dodgr` implements highly realistic and
fully-customisable profiles for routing through street networks with
various modes of transport, and using either distance- or time-based
routing. Routing can include such factors as waiting times at traffic
lights, delays for turning across oncoming traffic, and the effects of
elevation on both cyclists and pedestrians.

## Installation

You can install latest stable version of `dodgr` from CRAN with:

``` r
install.packages("dodgr") # current CRAN version
```

Alternatively, current development versions can be installed using any
of the following options:

``` r
# install.packages("remotes")
remotes::install_git("https://git.sr.ht/~mpadge/dodgr")
remotes::install_bitbucket("atfutures/dodgr")
remotes::install_gitlab("atfutures1/dodgr")
remotes::install_github("ATFutures/dodgr")
```

Then load with

``` r
library (dodgr)
packageVersion ("dodgr")
#> [1] '0.2.6'
```

## Usage: Sample Data and `dodgr` networks

To illustrate functionality, the package includes an example data set
containing the Open Street Map network for [Hampi,
India](https://www.openstreetmap.org/#map=15/15.3368/76.4601) (a
primarily pedestrian village in the middle of a large World Heritage
zone). These data are in [Simple Features
(`sf`)](https://cran.r-project.org/package=sf) format, as a collection
of `LINESTRING` objects. `dodgr` represents networks as a simple
rectangular graph, with each row representing an edge segment between
two points or vertices. `sf`-format objects can be converted to
equivalent `dodgr` representations with the `weight_streetnet()`
function:

``` r
class (hampi)
#> [1] "sf"         "data.frame"
dim (hampi)
#> [1] 203  15
graph <- weight_streetnet (hampi, wt_profile = "foot")
class (graph)
#> [1] "data.frame"      "dodgr_streetnet"
dim (graph)
#> [1] 5973   15
```

The `sf`-format network contained 203 `LINESTRING` objects, with the
`weight_streetnet()` function decomposing these into 5,973 distinct
edges, indicating that the `sf` representation had around 29 edges or
segments in each `LINESTRING` object. The `dodgr` network then looks
like this:

``` r
head (graph)
```

| geom\_num | edge\_id | from\_id   | from\_lon | from\_lat |   to\_id   |  to\_lon |  to\_lat |          d | d\_weighted | highway      | way\_id  | component |      time | time\_weighted |
| --------: | -------: | :--------- | --------: | --------: | :--------: | -------: | -------: | ---------: | ----------: | :----------- | :------- | --------: | --------: | -------------: |
|         1 |        1 | 339318500  |  76.47489 |  15.34169 | 339318502  | 76.47612 | 15.34173 | 132.442169 |   165.55271 | unclassified | 28565950 |         1 | 95.358362 |     119.197952 |
|         1 |        2 | 339318502  |  76.47612 |  15.34173 | 339318500  | 76.47489 | 15.34169 | 132.442169 |   165.55271 | unclassified | 28565950 |         1 | 95.358362 |     119.197952 |
|         1 |        3 | 339318502  |  76.47612 |  15.34173 | 2398958028 | 76.47621 | 15.34174 |   8.888670 |    11.11084 | unclassified | 28565950 |         1 |  6.399843 |       7.999803 |
|         1 |        4 | 2398958028 |  76.47621 |  15.34174 | 339318502  | 76.47612 | 15.34173 |   8.888670 |    11.11084 | unclassified | 28565950 |         1 |  6.399843 |       7.999803 |
|         1 |        5 | 2398958028 |  76.47621 |  15.34174 | 1427116077 | 76.47628 | 15.34179 |   9.326536 |    11.65817 | unclassified | 28565950 |         1 |  6.715106 |       8.393882 |
|         1 |        6 | 1427116077 |  76.47628 |  15.34179 | 2398958028 | 76.47621 | 15.34174 |   9.326536 |    11.65817 | unclassified | 28565950 |         1 |  6.715106 |       8.393882 |

The `geom_num` column maps directly onto the sequence of `LINESTRING`
objects within the `sf`-formatted data. The `highway` column is taken
directly from Open Street Map, and denotes the kind of “highway”
represented by each edge. The `component` column is an integer value
describing which of the connected components of the network each edge
belongs to (with `1` always being the largest component; `2` the second
largest; and so on).

Note that the `d_weighted` values are often greater than the geometric
distances, `d`. In the example shown, `service` highways are not ideal
for pedestrians, and so weighted distances are slightly greater than
actual distances. Compare this with:

``` r
head (graph [graph$highway == "path", ])
```

|    | geom\_num | edge\_id | from\_id   | from\_lon | from\_lat |   to\_id   |  to\_lon |  to\_lat |        d | d\_weighted | highway | way\_id  | component |     time | time\_weighted |
| -- | --------: | -------: | :--------- | --------: | --------: | :--------: | -------: | -------: | -------: | ----------: | :------ | :------- | --------: | -------: | -------------: |
| 47 |         2 |       47 | 338905220  |  76.47398 |  15.31224 | 338907543  | 76.47405 | 15.31241 | 19.70399 |    19.70399 | path    | 30643853 |         1 | 35.46718 |       35.46718 |
| 48 |         2 |       48 | 338907543  |  76.47405 |  15.31241 | 338905220  | 76.47398 | 15.31224 | 19.70399 |    19.70399 | path    | 30643853 |         1 | 35.46718 |       35.46718 |
| 49 |         2 |       49 | 338907543  |  76.47405 |  15.31241 | 2398957585 | 76.47410 | 15.31259 | 21.39172 |    21.39172 | path    | 30643853 |         1 | 38.50510 |       38.50510 |
| 50 |         2 |       50 | 2398957585 |  76.47410 |  15.31259 | 338907543  | 76.47405 | 15.31241 | 21.39172 |    21.39172 | path    | 30643853 |         1 | 38.50510 |       38.50510 |
| 51 |         2 |       51 | 2398957585 |  76.47410 |  15.31259 | 338907597  | 76.47413 | 15.31279 | 22.15205 |    22.15205 | path    | 30643853 |         1 | 39.87370 |       39.87370 |
| 52 |         2 |       52 | 338907597  |  76.47413 |  15.31279 | 2398957585 | 76.47410 | 15.31259 | 22.15205 |    22.15205 | path    | 30643853 |         1 | 39.87370 |       39.87370 |

A `"path"` offers ideal walking conditions, and so weighted distances
are equal to actual distances.

## Usage: Distances and Times

The many-to-many nature of `dodgr` means that the function to calculate
distances,
[`dodgr_distances()`](https://atfutures.github.io/dodgr/reference/dodgr_distances.html)
or, for street networks, times,
[`dodgr_times()`](https://atfutures.github.io/dodgr/reference/dodgr_times.html),
accepts two vectors or matrices of routing points as inputs (describing
origins and destinations), and returns a corresponding matrix of
pairwise distances. If an input graph has columns for both distances and
weighted distances, and/or times and weighted times, the weighted
versions are used to determine the effectively shortest or fastest
routes through a network, while actual distances or times are summed
along the routes to calculate final values. It is of course also
possible to calculate distances along fastest routes, times along
shortest routes, or any combination thereof, as detailed in the package
vignette on [street networks and time-based
routing](https://atfutures.github.io/dodgr/articles/times.html).

Routing points can, for example, be randomly selected from the vertices
of a graph. The vertices can in turn be extracted with the
`dodgr_vertices()` function:

``` r
v <- dodgr_vertices (graph)
head (v)
```

|    | id         |        x |        y | component | n |
| -- | :--------- | -------: | -------: | --------: | -: |
| 1  | 339318500  | 76.47489 | 15.34169 |         1 | 0 |
| 2  | 339318502  | 76.47612 | 15.34173 |         1 | 1 |
| 4  | 2398958028 | 76.47621 | 15.34174 |         1 | 2 |
| 6  | 1427116077 | 76.47628 | 15.34179 |         1 | 3 |
| 8  | 339318503  | 76.47641 | 15.34190 |         1 | 4 |
| 10 | 2398958034 | 76.47650 | 15.34199 |         1 | 5 |

For OSM data extracted with the `osmdata` package (or, equivalently, via
the `dodgr::dodgr_streetnet()` function), each object (vertices, ways,
and high-level relations between these objects) is assigned a unique
identifying number. These are retained both in `osmdata` and `dodgr`, as
the `way_id` column in the above `graph`, and as the `id` column in the
vertices. Random vertices may be generated in this case through
selecting `id` values:

``` r
from <- sample (v$id, size = 20)
to <- sample (v$id, size = 50)
d <- dodgr_dists (graph = graph, from = from, to = to)
dim (d)
#> [1] 20 50
```

Alternatively, the points may be specified as matrices of geographic
coordinates:

``` r
from_x <- min (graph$from_lon) + runif (20) * diff (range (graph$from_lon))
from_y <- min (graph$from_lat) + runif (20) * diff (range (graph$from_lat))
to_x <- min (graph$from_lon) + runif (50) * diff (range (graph$from_lon))
to_y <- min (graph$from_lat) + runif (50) * diff (range (graph$from_lat))
d <- dodgr_dists (graph = graph, from = cbind (from_x, from_y), to = cbind (to_x, to_y))
```

In this case, the random points will be mapped on to the nearest points
on the street network. This may, of course, map some points onto minor,
disconnected components of the graph. This can be controlled either by
reducing the graph to it’s largest connected component only:

``` r
graph <- graph [graph$component == 1, ]
nrow (graph)
```

or by explicitly using the `match_points_to_graph()` function with the
option `connected = TRUE`:

``` r
from <- match_points_to_graph (v, cbind (from_x, from_y), connected = TRUE)
to <- match_points_to_graph (v, cbind (to_x, to_y), connected = TRUE)
```

This function returns an index into the result of `dodgr_vertices`, and
so points to use for routing must then be extracted as follows:

``` r
from <- v$id [from] # or from <- v [from, c ("x", "y")]
to <- v$id [to]
d <- dodgr_dists (graph = graph, from = from, to = to)
```

## Usage: Flow Aggregation

Flow aggregation refers to the procedure of routing along multiple ways
according to specified densities of flow between defined origin and
destination points, and aggregating flows along each edge of the
network. The procedure is functionally similar to the above procedure
for distances, with the addition of a matrix specifying pairwise flow
densities between the input set of origin (`from`) and destination
(`to`) points. The following example illustrates use with a random “flow
matrix”:

``` r
flows <- array (runif (length (from) * length (to)), dim = c (length (from), length (to)))
length (from); length (to); dim (flows)
#> [1] 20
#> [1] 50
#> [1] 20 50
f <- dodgr_flows_aggregate (graph = graph, from = from, to = to, flows = flows)
```

The result is simply the input `graph` with an additional column
quantifying the aggregate flows along each edge:

``` r
head (f)
```

| geom\_num | edge\_id | from\_id   | from\_lon | from\_lat |   to\_id   |  to\_lon |  to\_lat |          d | d\_weighted | highway      | way\_id  | component |      time | time\_weighted |      flow |
| --------: | -------: | :--------- | --------: | --------: | :--------: | -------: | -------: | ---------: | ----------: | :----------- | :------- | --------: | --------: | -------------: | --------: |
|         1 |        1 | 339318500  |  76.47489 |  15.34169 | 339318502  | 76.47612 | 15.34173 | 132.442169 |   165.55271 | unclassified | 28565950 |         1 | 95.358362 |     119.197952 | 0.0822301 |
|         1 |        2 | 339318502  |  76.47612 |  15.34173 | 339318500  | 76.47489 | 15.34169 | 132.442169 |   165.55271 | unclassified | 28565950 |         1 | 95.358362 |     119.197952 | 0.0000000 |
|         1 |        3 | 339318502  |  76.47612 |  15.34173 | 2398958028 | 76.47621 | 15.34174 |   8.888670 |    11.11084 | unclassified | 28565950 |         1 |  6.399843 |       7.999803 | 0.0822301 |
|         1 |        4 | 2398958028 |  76.47621 |  15.34174 | 339318502  | 76.47612 | 15.34173 |   8.888670 |    11.11084 | unclassified | 28565950 |         1 |  6.399843 |       7.999803 | 0.0000000 |
|         1 |        5 | 2398958028 |  76.47621 |  15.34174 | 1427116077 | 76.47628 | 15.34179 |   9.326536 |    11.65817 | unclassified | 28565950 |         1 |  6.715106 |       8.393882 | 0.0822301 |
|         1 |        6 | 1427116077 |  76.47628 |  15.34179 | 2398958028 | 76.47621 | 15.34174 |   9.326536 |    11.65817 | unclassified | 28565950 |         1 |  6.715106 |       8.393882 | 0.0000000 |

An additional flow aggregation function can be applied in cases where
only densities at origin points are known, and movement throughout a
graph is dispersive:

``` r
f <- dodgr_flows_disperse (graph = graph, from = from, dens = runif (length (from)))
```

## Further detail

For more detail, see the [main package
vignette](https://atfutures.github.io/dodgr/articles/dodgr.html), and
the second vignette on [street networks and time-based
routing](https://atfutures.github.io/dodgr/articles/times.html)
