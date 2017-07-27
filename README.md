<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Build Status](https://travis-ci.org/mpadge/pqspr.svg)](https://travis-ci.org/mpadge/pqspr) [![Project Status: Concept - Minimal or no implementation has been done yet.](http://www.repostatus.org/badges/0.1.0/concept.svg)](http://www.repostatus.org/#concept) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/pqspr)](http://cran.r-project.org/web/packages/pqspr)

dodgr: Distances on Directed Graphs in R
========================================

R package for calculating pairwise distances on dual-weighted directed graphs using Priority Queue Shortest Paths. Dual-weighted directed graphs are directed graphs with two sets of weights so that `weight1(A->B) != weight1(B->A)`---the directed property---and `weight2(A->B) != weight1(A->B)`---the dual property. `dodgr` calculates shortest paths according to one weight, while distances along paths are calculating using the other weight. A canonical example of a dual-weighted directed graph is a street network to be used for routing. Routes are usually calculated by weighting different kinds of streets or ways according to a particular mode of transport, while the desired output is a direct distance which has to be calculated using a different set of weights.

The shortest path algorithm relies on priority queue code adapted from original code by Shane Saunders available here <http://www.cosc.canterbury.ac.nz/tad.takaoka/alg/spalgs/spalgs.html>

Installation
------------

You can install `dodgr` from github with:

``` r
# install.packages("devtools")
devtools::install_github("mpadge/pqspr")
```

Simple Usage
------------

The primary function,

``` r
d <- dodgr_dists (graph = graph, from = from)
```

produces a square matrix of distances between all points listed in `from` and routed along the dual-weighted directed network given in `graph`.

Detailed Usage
--------------

A graph or network in `dodgr` is represented as a flat table (`data.frame`, `tibble`, `data.table`, whatever) of minimally four columns: `from`, `to`, `weight`, and `distance`. The first two can be of arbitrary form (`numeric` or `character`); `weight` is used to evaluate the shortest paths, and the desired distances are evaluated by summing the values of `distance` along those paths. For a street network example, `weight` will generally be the actual distance multiplied by a priority weighting for a given mode of transport and type of way, while `distance` will be the pysical distance.

The primary function of `dodgr` is `dodgr_dists()`, which accepts the three main arguments of

1.  `graph` - the flat table described above;
2.  `from` - a list of origin points for which distances are to be calcualted; and
3.  `to` - a list of equivalent destination points. If missing, pairwise distances are calculated between all points in `from`.

For spatial graphs in which `graph` has additional longitude and latitude coordinates, `from` and `to` may also be given as lists of coordinates; otherwise they must match onto names or values given in `graph$from` and `graph$to`. Coordinates in `from` and `to` that do not exactly match on to values given in `graph` will be mapped on to the nearest `graph` nodes.
