#include <string>
#include <cmath>

#include <Rcpp.h>

const float INFINITE_FLOAT =  std::numeric_limits<float>::max ();
const double INFINITE_DOUBLE =  std::numeric_limits<double>::max ();
const int INFINITE_INT =  std::numeric_limits<int>::max ();

double haversine (double x1, double y1, double x2, double y2);
Rcpp::List rcpp_sf_as_network (const Rcpp::List &sf_lines,
        const Rcpp::DataFrame &pr);
Rcpp::NumericVector rcpp_points_index (const Rcpp::DataFrame &xy,
        Rcpp::DataFrame &pts);
