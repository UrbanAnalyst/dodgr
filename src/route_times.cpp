
#include "route_times.h"


// https://stackoverflow.com/questions/6989100/sort-points-in-clockwise-order

bool routetimes::isLess (double x0, double y0, double xa, double ya,
        double xb, double yb)
{
    if (xa >= 0.0 && xb < 0.0)
        return true;
    else if (xa == 0.0 && xb == 0.0)
        return ya > yb;

    double det = (xa - x0) * (yb - y0) - (xb - x0) * (ya - y0);
    if (det < 0)
        return true;
    if (det > 0)
        return false;

    double d1 = (xa - x0) * (xa - x0) - (ya - x0) * (ya - y0);
    double d2 = (xb - x0) * (xb - x0) - (yb - x0) * (yb - y0);
    return d1 > d2;
}

//' rcpp_route_times
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_route_times (const Rcpp::DataFrame graph)
{
    Rcpp::NumericMatrix res;

    return res;
}
