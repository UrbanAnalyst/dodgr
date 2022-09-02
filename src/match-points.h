#include <string>
#include <cmath>

#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

constexpr float INFINITE_FLOAT =  std::numeric_limits<float>::max ();
constexpr double INFINITE_DOUBLE =  std::numeric_limits<double>::max ();
constexpr int INFINITE_INT =  std::numeric_limits<int>::max ();

int which_side_of_line (const double ax, const double ay,
        const double bx, const double by, const double x, const double y);

Rcpp::IntegerVector rcpp_points_index (const Rcpp::DataFrame &xy,
        Rcpp::DataFrame &pts);
Rcpp::NumericVector rcpp_points_to_edges_par (const Rcpp::DataFrame &graph,
        Rcpp::DataFrame &pts);
