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

R package for calculating pairwise distances on dual-weighted directed
graphs using Priority Queue Shortest Paths. Dual-weighted directed
graphs are directed graphs with two sets of weights so that
`weight1(A->B) != weight1(B->A)`—the directed property—and
`weight2(A->B) != weight1(A->B)`—the dual property. `dodgr` calculates
shortest paths according to one weight, while distances along paths are
calculating using the other weight. A canonical example of a
dual-weighted directed graph is a street network to be used for routing.
Routes are usually calculated by weighting different kinds of streets or
ways according to a particular mode of transport, while the desired
output is a direct, unweighted distance.

## Installation

You can install `dodgr` from github with:

``` r
# install.packages("devtools")
devtools::install_github("ATFutures/dodgr")
```

## Usage

The primary function,

``` r
d <- dodgr_dists (graph = graph, from = pts, to = pts)
```

produces a square matrix of distances between all points listed in `pts`
and routed along the dual-weighted directed network given in `graph`. An
even simpler usage allows calculation of pair-wise distances between a
set of geographical coordinates (here, for a sizey chunk of New York
City):

``` r
xlim <- c (-74.12931, -73.99214)
ylim <- c (40.70347, 40.75354)
npts <- 1000
pts <- data.frame (x = xlim [1] + runif (npts) * diff (xlim),
                   y = ylim [1] + runif (npts) * diff (ylim))
st <- Sys.time ()
d <- dodgr_dists (from = pts)
Sys.time () - st
#> Time difference of 9.473597 secs
range (d, na.rm = TRUE)
#> [1]  0.00000 16.64285
```

This will automatically download the street network (using
[`osmdata`](https://cran.r-project.org/package=osmdata)), and even then
calculating distances between 1,000 points – that’s 1,000,000 pairwise
distances\! – can be done in around 10 seconds.

### Other Functions

The other main functions of `dodgr` are `dodgr_paths`, to return the
actual paths between pairs of vertices, and `dodgr_flows` to aggregate
flows as routed throughout a network according to sets of start and end
points (origins and destinations), and associated densities or numbers
travelling between each of these.

### The `dodgr` graph structure

A graph or network in `dodgr` is represented as a flat table
(`data.frame`, `tibble`, `data.table`, whatever) of minimally four
columns: `from`, `to`, `weight`, and `distance`. The first two can be of
arbitrary form (`numeric` or `character`); `weight` is used to evaluate
the shortest paths, and the desired distances are evaluated by summing
the values of `distance` along those paths. For a street network
example, `weight` will generally be the actual distance multiplied by a
priority weighting for a given mode of transport and type of way, while
`distance` will be the pysical distance.

`dodgr` includes the conversion functions:

1.  `dodgr_to_sfc` to convert spatial `dodgr` graphs into Simple
    Features format used by the [`sf`
    package](https://cran.r-project.org/package=sf).
2.  `dodgr_to_igraph` to convert (not necessarily spatial) `dodgr`
    graphs into [`igraph`](https://cran.r-project.org/package=igraph)
    format (not yet implemented); and
3.  `dodgr_to_tidygraph` to convert (not necessarily spatial) `dodgr`
    graphs into
    [`tidygraph`](https://cran.r-project.org/package=tidygraph) format
    (not yet implemented).

### Further detail

For more detail, see the [main package
vignette](https://atfutures.github.io/dodgr/articles/dodgr.html), along
with a second vignette detailing [benchmark
timings](https://atfutures.github.io/dodgr/articles/benchmark.html),
showing that under many circumstances, `dodgr` performs considerably
faster than equivalent routines from the `igraph` package.
