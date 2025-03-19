#include "dist-max.h"

//' rcpp_dist_max
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_dist_max (const Rcpp::DataFrame &graph)
{
    const size_t n = static_cast <size_t> (graph.nrow ());
    Rcpp::NumericVector result (4L);
    double dmax = 0;

    std::vector <double> x1 = graph ["from_lon"];
    std::vector <double> y1 = graph ["from_lat"];
    std::vector <double> x2 = graph ["to_lon"];
    std::vector <double> y2 = graph ["to_lat"];

    for (size_t i = 0; i < (n - 1L); i++) {
        for (size_t j = (i + 1L); j < n; j++) {

            const double d = fabs (x1 [i] - x1 [j]) + fabs (x2 [i] - x2 [j]);
            if (d > dmax) {
                dmax = d;
                result = {x1 [i], y1 [i], x2 [j], y2 [j]};
            }
        }
    }


    return result;
}
