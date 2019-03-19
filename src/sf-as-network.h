#include <string>
#include <cmath>

#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

const double earth = 6378.137; // WSG84 definition

constexpr float INFINITE_FLOAT =  std::numeric_limits<float>::max ();
constexpr double INFINITE_DOUBLE =  std::numeric_limits<double>::max ();
constexpr int INFINITE_INT =  std::numeric_limits<int>::max ();

double haversine (double x1, double y1, double x2, double y2);
Rcpp::List rcpp_sf_as_network (const Rcpp::List &sf_lines,
        const Rcpp::DataFrame &pr);
Rcpp::IntegerVector rcpp_points_index (const Rcpp::DataFrame &xy,
        Rcpp::DataFrame &pts);
Rcpp::IntegerVector rcpp_points_index_par (const Rcpp::DataFrame &xy,
        Rcpp::DataFrame &pts);
