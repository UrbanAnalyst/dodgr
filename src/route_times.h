#pragma once

#include <vector>

#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

namespace routetimes {
bool isLess (double x0, double y0, double xa, double ya,
        double xb, double yb);
} // end namespace

Rcpp::NumericMatrix rcpp_route_times (const Rcpp::DataFrame graph);
