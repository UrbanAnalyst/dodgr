#pragma once

#include <vector>

#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

Rcpp::NumericMatrix rcpp_route_times (const Rcpp::DataFrame graph);
