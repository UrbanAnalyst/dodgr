#include "match-points.h"

//' Determine which side of intersecting line a point lies on.
//'
//' @param (ax, ay) Coordinates of one end of line
//' @param (bx, by) Coordinates of other end of line
//' @param (x, y) Coordinates of point
//' @return 0 if point on line, -1 if to the left; +1 if to the right.
//'
//' @noRd
int which_side_of_line (const double ax, const double ay,
        const double bx, const double by, const double x, const double y)
{
    double val = (bx - ax) * (y - ay) - (by - ay) * (x - ax);
    return (val > 0) ? 1 : ((val < 0) ? -1 : 0);
}

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
    const size_t nxy, npts;
    RcppParallel::RVector <double> index;

    // constructor
    OneEdgeIndex (
            const RcppParallel::RVector <double> pt_x_in,
            const RcppParallel::RVector <double> pt_y_in,
            const RcppParallel::RVector <double> xfr_in,
            const RcppParallel::RVector <double> yfr_in,
            const RcppParallel::RVector <double> xto_in,
            const RcppParallel::RVector <double> yto_in,
            const size_t nxy_in,
            const size_t npts_in,
            Rcpp::NumericVector index_in) :
        pt_x (pt_x_in), pt_y (pt_y_in),
        xfr (xfr_in), yfr (yfr_in), xto (xto_in), yto (yto_in),
        nxy (nxy_in), npts (npts_in), index (index_in)
    {
    }

    // Parallel function operator
    void operator() (std::size_t begin, std::size_t end)
    {
        for (std::size_t i = begin; i < end; i++)
        {
            double dmin = INFINITE_DOUBLE;
            double x_intersect = INFINITE_DOUBLE;
            double y_intersect = INFINITE_DOUBLE;
            long int jmin = INFINITE_INT;
            int which_side = INFINITE_INT;

            for (size_t j = 0; j < nxy; j++)
            {
                const double x1 = xfr [j], y1 = yfr [j];
                const double x2 = xto [j], y2 = yto [j];

                const double px = x2 - x1;
                const double py = y2 - y1;

                const double norm = px * px + py * py;

                double u = ((pt_x [i] - x1) * px + (pt_y [i] - y1) * py) / norm;
                u = (u > 1) ? 1 : ((u < 0) ? 0 : u);

                const double xx = x1 + u * px;
                const double yy = y1 + u * py;

                const double dx = xx - pt_x [i];
                const double dy = yy - pt_y [i];

                const double dij = sqrt (dx * dx + dy * dy);

                if (dij < dmin)
                {
                    dmin = dij;
                    jmin = static_cast <long int> (j);
                    which_side = which_side_of_line (xfr [j], yfr [j],
                            xto [j], yto [j], pt_x [i], pt_y [i]);
                    x_intersect = xx;
                    y_intersect = yy;
                }
            }
            index [i] = static_cast <double> (jmin);
            index [i + npts] = x_intersect;
            index [i + 2L * npts] = y_intersect;
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
Rcpp::NumericVector rcpp_points_to_edges_par (const Rcpp::DataFrame &graph,
        Rcpp::DataFrame &pts)
{
    Rcpp::NumericVector ptx = pts ["x"];
    Rcpp::NumericVector pty = pts ["y"];

    Rcpp::NumericVector xfr = graph ["xfr"];
    Rcpp::NumericVector yfr = graph ["yfr"];
    Rcpp::NumericVector xto = graph ["xto"];
    Rcpp::NumericVector yto = graph ["yto"];

    const size_t nxy = static_cast <size_t> (graph.nrow ()),
        npts = static_cast <size_t> (pts.nrow ());

    // index holds three vectors:
    // 1. the integer index
    // 2. the x coordinates of the intersection points
    // 3. the y coordinates of the intersection points
    Rcpp::NumericVector index (npts * 3L);

    // Create parallel worker
    OneEdgeIndex one_edge_indx (RcppParallel::RVector <double> (ptx),
            RcppParallel::RVector <double> (pty),
            RcppParallel::RVector <double> (xfr),
            RcppParallel::RVector <double> (yfr),
            RcppParallel::RVector <double> (xto),
            RcppParallel::RVector <double> (yto),
            nxy, npts, index);

    RcppParallel::parallelFor (0, npts, one_edge_indx);

    return index;
}
