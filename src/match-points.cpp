#include "match-points.h"

//' Simple match of points to nearest vertices
//' @noRd
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

//' Match points to nearest edge of graph at which perpendicular from point
//' bisects edges. Uses psuedo-code from
//' https://stackoverflow.com/a/6853926
//' @noRd
struct OneEdgeIndex : public RcppParallel::Worker
{
    const RcppParallel::RVector <double> pt_x, pt_y,
          xfr, yfr, xto, yto;
    const size_t nxy;
    RcppParallel::RVector <int> index;

    // constructor
    OneEdgeIndex (
            const RcppParallel::RVector <double> pt_x_in,
            const RcppParallel::RVector <double> pt_y_in,
            const RcppParallel::RVector <double> xfr_in,
            const RcppParallel::RVector <double> yfr_in,
            const RcppParallel::RVector <double> xto_in,
            const RcppParallel::RVector <double> yto_in,
            const size_t nxy_in,
            Rcpp::IntegerVector index_in) :
        pt_x (pt_x_in), pt_y (pt_y_in),
        xfr (xfr_in), yfr (yfr_in), xto (xto_in), yto (yto_in),
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
                const double x1 = xfr [j], y1 = yfr [j];
                const double x2 = xto [j], y2 = yto [j];

                const double A = pt_x [i] - x1;
                const double B = pt_y [i] - y1;
                const double C = x2 - x1;
                const double D = y2 - y1;

                const double dot = A * C + B * D;
                const double len_sq = C * C + D * D;
                double param = -1.0;
                if (fabs (len_sq) < 1.0e-12)
                {
                    param = dot / len_sq;
                }

                double xx, yy;
                if (param < 0.0)
                {
                    xx = x1;
                    yy = y1;
                } else if (param > 1.0)
                {
                    xx = x2;
                    yy = y2;
                } else
                {
                    xx = x1 + param * C;
                    yy = y1 + param * D;
                }

                const double dx = pt_x [i] - xx;
                const double dy = pt_y [i] - yy;

                const double dij = sqrt (dx * dx + dy * dy);

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

//' rcpp_points_to_edges_par
//'
//' Get index of nearest edges to list of points
//'
//' @param graph Rcpp::DataFrame containing the full edge-based graph
//' @param pts Rcpp::DataFrame containing the points to be matched
//'
//' @return 0-indexed Rcpp::NumericVector index into graph of nearest points
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::IntegerVector rcpp_points_to_edges_par (const Rcpp::DataFrame &graph,
        Rcpp::DataFrame &pts)
{
    Rcpp::NumericVector ptx = pts ["x"];
    Rcpp::NumericVector pty = pts ["y"];

    Rcpp::NumericVector xfr = graph ["xfr"];
    Rcpp::NumericVector yfr = graph ["yfr"];
    Rcpp::NumericVector xto = graph ["xto"];
    Rcpp::NumericVector yto = graph ["yto"];

    const size_t npts = static_cast <size_t> (pts.nrow ()),
           nxy = static_cast <size_t> (graph.nrow ());

    Rcpp::IntegerVector index (npts);

    // Create parallel worker
    OneEdgeIndex one_edge_indx (RcppParallel::RVector <double> (ptx),
            RcppParallel::RVector <double> (pty),
            RcppParallel::RVector <double> (xfr),
            RcppParallel::RVector <double> (yfr),
            RcppParallel::RVector <double> (xto),
            RcppParallel::RVector <double> (yto),
            nxy, index);

    RcppParallel::parallelFor (0, npts, one_edge_indx);

    return index;
}
