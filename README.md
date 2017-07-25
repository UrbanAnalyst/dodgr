<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Build Status](https://travis-ci.org/mpadge/pqspr.svg)](https://travis-ci.org/mpadge/pqspr) [![Project Status: Concept - Minimal or no implementation has been done yet.](http://www.repostatus.org/badges/0.1.0/concept.svg)](http://www.repostatus.org/#concept) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/pqspr)](http://cran.r-project.org/web/packages/pqspr) ![downloads](http://cranlogs.r-pkg.org/badges/grand-total/pqspr)

pqspr
=====

Priority Queue Shortest Paths in R. Simply because they are blindingly faster than anything else, yet do not seem to be implemented in any current package. Note that most of the source code which does the work here is only lightly adapted from original code of Shane Saunders obtained from <http://www.cosc.canterbury.ac.nz/tad.takaoka/alg/spalgs/spalgs.html>

Installation
------------

You can install `pqspr` from github with:

``` r
# install.packages("devtools")
devtools::install_github("mpadge/pqspr")
```

Timing Comparison
-----------------

Create a street network using `osmdata` and some code currently in the github repo [`osmprob`](https://github.com/osm-router/osmprob) which is installed first.

``` r
devtools::install_github ("osm-router/osmprob")
```

Use `osmprob` to get a street network graph:

``` r
start_pt <- c (-74.00150, 40.74178)
end_pt <- c (-74.07889, 40.71113)
graph <- download_graph (start_pt, end_pt)
```

The result has a `$original` graph containing all vertices and a `$compact` graph with redundant vertices removed. We'll perform routing on the latter which has this many edges:

``` r
nrow (graph$compact)
#> [1] 8177
```

Both of these graphs are simple data frames detailing all edges. Some light wrangling is now necessary to prepare both the `igraph` and the simple structure submitted to the `pqspr` routine:

``` r
edges <- cbind (graph$compact$from_id, graph$compact$to_id)
nodes <- unique (as.vector (edges)) # used below in test comparison
edges <- as.vector (t (edges))
igr <- igraph::make_directed_graph (edges)
igraph::E (igr)$weight <- graph$compact$d_weighted

graph <- graph$compact
indx <- which (names (graph) %in% c ("from_id", "to_id", "d", "d_weighted"))
graph <- graph [, indx]
graph$from_id <- paste0 (graph$from_id)
graph$to_id <- paste0 (graph$to_id)
```

Then the final timing comparison between `igraph::distances()`, which returns a matrix of distances between all vertices, and the equivalent and only function of this package, `test()`:

``` r
rbenchmark::benchmark (test (graph),
                   igraph::distances (igr, v = nodes, to = nodes, mode = "out"),
                       replications = 10)
#>                                                          test replications
#> 2 igraph::distances(igr, v = nodes, to = nodes, mode = "out")           10
#> 1                                                 test(graph)           10
#>   elapsed relative user.self sys.self user.child sys.child
#> 2  18.267   10.541    18.176    0.092          0         0
#> 1   1.733    1.000     1.664    0.068          0         0
```

And these priority queue routines are over ten times faster than the `igraph` equivalent.
