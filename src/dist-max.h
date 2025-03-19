#pragma once

#include <vector>

#include <Rcpp.h>

struct OneEdge {
    double x, y;
};

Rcpp::NumericVector rcpp_dist_max (const Rcpp::DataFrame &graph);
