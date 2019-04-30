#include "sf-as-network.h"

// Haversine great circle distance between two points
double sf::haversine (double x1, double y1, double x2, double y2)
{
    double xd = (x2 - x1) * M_PI / 180.0;
    double yd = (y2 - y1) * M_PI / 180.0;
    double d = sin (yd / 2.0) * sin (yd / 2.0) + cos (y2 * M_PI / 180.0) *
        cos (y1 * M_PI / 180.0) * sin (xd / 2.0) * sin (xd / 2.0);
    d = 2.0 * earth * asin (sqrt (d));
    return (d);
}

//' rcpp_sf_as_network
//'
//' Return OSM data from Simple Features format input
//'
//' @param sf_lines An sf collection of LINESTRING objects
//' @param pr Rcpp::DataFrame containing the weighting profile
//'
//' @return Rcpp::List objects of OSM data, one matrix of numeric and one of
//' character values. The former contain 7 columns:
//' 1. sf geom index
//' 2. from longitude
//' 3. from latitude
//' 4. to longitude
//' 5. to latitude
//' 6. distance
//' 7. weighted_distance
//' The character value matrix  has 4 columns of:
//' 1. from ID
//' 2. to ID
//' 3. highway type
//' 4. OSM way ID
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::List rcpp_sf_as_network (const Rcpp::List &sf_lines,
        const Rcpp::DataFrame &pr)
{
    std::unordered_map <std::string, double> profile;
    Rcpp::StringVector hw = pr ["way"];
    Rcpp::NumericVector val = pr ["value"];
    if (hw.size () > 1) // single NA -> 0-length RcppVector
    {
        for (int i = 0; i != hw.size (); i ++)
            profile.insert (std::make_pair (std::string (hw [i]), val [i]));
    } 


    Rcpp::CharacterVector nms = sf_lines.attr ("names");
    int one_way_index = -1;
    int highway_index = -1;
    int geom_index = -1;
    for (R_xlen_t i = 0; i < nms.size (); i++)
    {
        if (nms [i] == "oneway")
            one_way_index = static_cast <int> (i);
        if (nms [i] == "highway")
            highway_index = static_cast <int> (i);
        if (nms [i] == "geometry")
            geom_index = static_cast <int> (i);
    }
    if (geom_index < 0)
        throw std::runtime_error ("sf_lines have no geometry component");

    Rcpp::CharacterVector ow; // init length = 0
    Rcpp::CharacterVector highway;
    bool has_oneway = false;
    if (one_way_index >= 0)
    {
        ow = sf_lines [one_way_index];
        has_oneway = true;
    }
    if (highway_index >= 0)
        highway = sf_lines [highway_index];

    Rcpp::List geoms = sf_lines [geom_index];

    std::vector <std::string> att_names = geoms.attributeNames ();
    bool has_names = false;
    for (std::vector <std::string>::iterator it = att_names.begin ();
            it != att_names.end (); it++)
        if (*it == "names")
            has_names = true;
    std::vector <std::string> way_names;
    if (has_names)
        way_names = Rcpp::as <std::vector <std::string> > (geoms.attr ("names"));

    std::vector <bool> isOneWay (static_cast <size_t> (geoms.length ()));
    std::fill (isOneWay.begin (), isOneWay.end (), false);
    // Get dimension of matrix
    size_t nrows = 0;
    unsigned int ngeoms = 0;
    for (auto g = geoms.begin (); g != geoms.end (); ++g)
    {
        // Rcpp uses an internal proxy iterator here, NOT a direct copy
        Rcpp::NumericMatrix gi = (*g);
        size_t rows = static_cast <size_t> (gi.nrow () - 1);
        nrows += rows;
        if (has_oneway && (
                    ow [ngeoms] == "yes" || ow [ngeoms] == "true" ||
                    ow [ngeoms] == "Yes" || ow [ngeoms] == "True" ||
                    ow [ngeoms] == "YES" || ow [ngeoms] == "TRUE"))
            isOneWay [ngeoms] = true;
        else
            nrows += rows;
        ngeoms ++;
    }

    Rcpp::NumericMatrix nmat = Rcpp::NumericMatrix (Rcpp::Dimension (nrows, 7));
    Rcpp::CharacterMatrix idmat =
        Rcpp::CharacterMatrix (Rcpp::Dimension (nrows, 4));

    nrows = 0;
    ngeoms = 0;
    int fake_id = 0;
    for (auto g = geoms.begin (); g != geoms.end (); ++ g)
    {
        Rcpp::checkUserInterrupt ();
        Rcpp::NumericMatrix gi = (*g);
        std::string hway;
        double hw_factor = 1.0;
        if (profile.size () > 0)
        {
            hway = std::string (highway [ngeoms]);
            hw_factor = profile [hway];
            if (hw_factor > 0.0)
                hw_factor = 1.0 / hw_factor;
        }

        Rcpp::List ginames = gi.attr ("dimnames");
        Rcpp::CharacterVector rnms;
        if (ginames.length () > 0)
        {
            if (!Rf_isNull (ginames [0]))
                rnms = ginames [0];
        }
        if (rnms.size () == 0)
        {
            rnms = Rcpp::CharacterVector (gi.nrow ());
            for (int i = 0; i < gi.nrow (); i ++)
                rnms [i] = fake_id ++;
        }
        if (rnms.size () != gi.nrow ())
            throw std::runtime_error ("geom size differs from rownames");

        for (unsigned int i = 1;
                i < static_cast <unsigned int> (gi.nrow ()); i ++)
        {
            double d = sf::haversine (gi (i-1, 0), gi (i-1, 1), gi (i, 0),
                    gi (i, 1));
            sf::fill_one_row (ngeoms, gi, rnms, d, hw_factor, hway,
                    has_names, way_names, i, nrows, false, nmat, idmat);
            nrows++;

            if (!isOneWay [ngeoms])
            {
                sf::fill_one_row (ngeoms, gi, rnms, d, hw_factor, hway,
                        has_names, way_names, i, nrows, true, nmat, idmat);
                nrows++;
            }
        }
        ngeoms ++;
    }

    return Rcpp::List::create (
            Rcpp::Named ("numeric_values") = nmat,
            Rcpp::Named ("character_values") = idmat);
}

void sf::fill_one_row (const unsigned int ngeoms, const Rcpp::NumericMatrix &gi,
        const Rcpp::CharacterVector &rnms,
        const double &d, const double &hw_factor,
        const std::string &hway, const bool &has_names,
        const std::vector <std::string> &way_names,
        const size_t &grownum, const size_t &rownum, const bool &reverse,
        Rcpp::NumericMatrix &nmat, Rcpp::CharacterMatrix &idmat)
{
    size_t i_min_1 = grownum - 1, i = grownum;
    if (reverse)
    {
        i_min_1 = grownum;
        i = grownum - 1;
    }

    nmat (rownum, 0) = ngeoms;
    nmat (rownum, 1) = gi (i_min_1, 0);
    nmat (rownum, 2) = gi (i_min_1, 1);
    nmat (rownum, 3) = gi (i, 0);
    nmat (rownum, 4) = gi (i, 1);
    nmat (rownum, 5) = d;
    if (hw_factor > 0.0)
        nmat (rownum, 6) = d * hw_factor;
    else
        nmat (rownum, 6) = -1.0;

    idmat (rownum, 0) = rnms (i_min_1);
    idmat (rownum, 1) = rnms (i);
    idmat (rownum, 2) = hway;
    if (has_names)
        idmat (rownum, 3) = way_names [ngeoms];
}

struct OnePointIndex : public RcppParallel::Worker
{
    const Rcpp::NumericVector xy_x, xy_y, pt_x, pt_y;
    const size_t nxy;
    RcppParallel::RVector <int> index;

    // constructor
    OnePointIndex (
            const Rcpp::NumericVector xy_x_in,
            const Rcpp::NumericVector xy_y_in,
            const Rcpp::NumericVector pt_x_in,
            const Rcpp::NumericVector pt_y_in,
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
            long int li = static_cast <long int> (i);
            double dmin = INFINITE_DOUBLE;
            int jmin = INFINITE_INT;
            for (int j = 0; j < static_cast <int> (nxy); j++)
            {
                double dij = (xy_x [j] - pt_x [li]) * (xy_x [j] - pt_x [li]) +
                    (xy_y [j] - pt_y [li]) * (xy_y [j] - pt_y [li]);
                if (dij < dmin)
                {
                    dmin = dij;
                    jmin = j;
                }
            }
            index [i] = jmin;
        }
    }
                                   
};

//' rcpp_points_index_par
//'
//' Get index of nearest vertices to list of points
//'
//' @param graph Rcpp::DataFrame containing the graph
//' @param pts Rcpp::DataFrame containing the routing points
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
    OnePointIndex one_pt_indx (vtx, vty, ptx, pty, nxy, index);

    RcppParallel::parallelFor (0, npts, one_pt_indx);

    return index;
}
