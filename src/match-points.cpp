#include "match-points.h"

struct OnePointIndex : public RcppParallel::Worker
{
    const RcppParallel::RVector <double> xy_x, xy_y, pt_x, pt_y;
    const size_t nxy;
    RcppParallel::RVector <int> index;

    // constructor
    OnePointIndex (
            const RcppParallel::RVector <double> xy_x_in,
            const RcppParallel::RVector <double> xy_y_in,
            const RcppParallel::RVector <double> pt_x_in,
            const RcppParallel::RVector <double> pt_y_in,
            const size_t nxy_in,
            Rcpp::IntegerVector index_in) :
        xy_x (xy_x_in), xy_y (xy_y_in), pt_x (pt_x_in), pt_y (pt_y_in),
        nxy (nxy_in), index (index_in)
    {
    }

    // Parallel function operator
    void operator() (std::size_t begin, std::size_t end)
    {
        for (std::size_t i = begin; i < end; i++)
        {
            double dmin = INFINITE_DOUBLE;
            long int jmin = INFINITE_INT;
            for (size_t j = 0; j < nxy; j++)
            {
                double dij = (xy_x [j] - pt_x [i]) * (xy_x [j] - pt_x [i]) +
                    (xy_y [j] - pt_y [i]) * (xy_y [j] - pt_y [i]);
                if (dij < dmin)
                {
                    dmin = dij;
                    jmin = static_cast <long int> (j);
                }
            }
            index [i] = static_cast <int> (jmin);
        }
    }
                                   
};

//' rcpp_points_index_par
//'
//' Get index of nearest vertices to list of points
//'
//' @param xy Rcpp::DataFrame containing the vertex coordinates of the graph
//' @param pts Rcpp::DataFrame containing the points to be matched
//'
//' @return 0-indexed Rcpp::NumericVector index into graph of nearest points
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::IntegerVector rcpp_points_index_par (const Rcpp::DataFrame &xy,
        Rcpp::DataFrame &pts)
{
    Rcpp::NumericVector ptx = pts ["x"];
    Rcpp::NumericVector pty = pts ["y"];

    Rcpp::NumericVector vtx = xy ["x"];
    Rcpp::NumericVector vty = xy ["y"];

    size_t npts = static_cast <size_t> (pts.nrow ()),
           nxy = static_cast <size_t> (xy.nrow ());

    //Rcpp::IntegerVector index (n, Rcpp::IntegerVector::get_na ());
    Rcpp::IntegerVector index (npts);
    // Create parallel worker
    OnePointIndex one_pt_indx (RcppParallel::RVector <double> (vtx),
            RcppParallel::RVector <double> (vty),
            RcppParallel::RVector <double> (ptx),
            RcppParallel::RVector <double> (pty), nxy, index);

    RcppParallel::parallelFor (0, npts, one_pt_indx);

    return index;
}
