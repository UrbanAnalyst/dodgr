#include "sf-as-network.h"

// Haversine great circle distance between two points
double haversine (double x1, double y1, double x2, double y2)
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
//' @return Rcpp::List objects of OSM data
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::List rcpp_sf_as_network (const Rcpp::List &sf_lines,
        const Rcpp::DataFrame &pr)
{
    std::map <std::string, double> profile;
    //Rcpp::StringVector hw = pr [1];
    //Rcpp::NumericVector val = pr [2];
    Rcpp::StringVector hw = pr ["way"];
    Rcpp::NumericVector val = pr ["value"];
    if (hw.size () > 1)
    {
        for (int i = 0; i != hw.size (); i ++)
            profile.insert (std::make_pair (std::string (hw [i]), val [i]));
    }

    Rcpp::CharacterVector nms = sf_lines.attr ("names");
    int one_way_index = -1;
    int one_way_bicycle_index = -1;
    int highway_index = -1;
    int geom_index = -1;
    for (R_xlen_t i = 0; i < nms.size (); i++)
    {
        if (nms [i] == "oneway")
            one_way_index = static_cast <int> (i);
        if (nms [i] == "oneway.bicycle")
            one_way_bicycle_index = static_cast <int> (i);
        if (nms [i] == "highway")
            highway_index = static_cast <int> (i);
        if (nms [i] == "geometry")
            geom_index = static_cast <int> (i);
    }
    if (geom_index < 0)
        throw std::runtime_error ("sf_lines have no geometry component");
    Rcpp::CharacterVector ow; // init length = 0
    Rcpp::CharacterVector owb;
    Rcpp::CharacterVector highway;
    if (one_way_index >= 0)
        ow = sf_lines [one_way_index];
    if (one_way_bicycle_index >= 0)
        owb = sf_lines [one_way_bicycle_index];
    if (highway_index >= 0)
        highway = sf_lines [highway_index];
    if (ow.size () > 0)
    {
        if (ow.size () == owb.size ())
        {
            for (R_xlen_t i = 0; i != ow.size (); ++ i)
                if (ow [i] == "NA" && owb [i] != "NA")
                    ow [i] = owb [i];
        } else if (owb.size () > ow.size ())
            ow = owb;
    }

    //Rcpp::List geoms = sf_lines [nms.size () - 1];
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
        if (ngeoms < static_cast <unsigned int> (ow.size ()))
        {
            //if (!(ow [ngeoms] == "yes" || ow [ngeoms] == "-1"))
            if (ow [ngeoms] == "yes" || ow [ngeoms] == "true")
                isOneWay [ngeoms] = true;
            else
                nrows += rows;
        } else if (ow.size () == 0)
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
            double d = haversine (gi (i-1, 0), gi (i-1, 1), gi (i, 0),
                    gi (i, 1));
            nmat (nrows, 0) = ngeoms;
            nmat (nrows, 1) = gi (i-1, 0);
            nmat (nrows, 2) = gi (i-1, 1);
            nmat (nrows, 3) = gi (i, 0);
            nmat (nrows, 4) = gi (i, 1);
            nmat (nrows, 5) = d;
            if (hw_factor > 0.0)
                nmat (nrows, 6) = d * hw_factor;
            else
                nmat (nrows, 6) = -1.0;
            idmat (nrows, 0) = rnms (i-1);
            idmat (nrows, 1) = rnms (i);
            idmat (nrows, 2) = hway;
            if (has_names)
                idmat (nrows, 3) = way_names [ngeoms];
            nrows ++;
            if (!isOneWay [ngeoms])
            {
                nmat (nrows, 0) = ngeoms;
                nmat (nrows, 1) = gi (i, 0);
                nmat (nrows, 2) = gi (i, 1);
                nmat (nrows, 3) = gi (i-1, 0);
                nmat (nrows, 4) = gi (i-1, 1);
                nmat (nrows, 5) = d;
                if (hw_factor > 0.0)
                    nmat (nrows, 6) = d * hw_factor;
                else
                    nmat (nrows, 6) = -1.0;
                idmat (nrows, 0) = rnms (i);
                idmat (nrows, 1) = rnms (i-1);
                idmat (nrows, 2) = hway;
                if (has_names)
                    idmat (nrows, 3) = way_names [ngeoms];
                nrows ++;
            }
        }
        ngeoms ++;
    }

    return Rcpp::List::create (
            Rcpp::Named ("numeric_values") = nmat,
            Rcpp::Named ("character_values") = idmat);
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

//' rcpp_points_index
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
Rcpp::IntegerVector rcpp_points_index (const Rcpp::DataFrame &xy,
        Rcpp::DataFrame &pts)
{
    Rcpp::NumericVector ptx = pts ["x"];
    Rcpp::NumericVector pty = pts ["y"];

    Rcpp::NumericVector vtx = xy ["x"];
    Rcpp::NumericVector vty = xy ["y"];

    Rcpp::IntegerVector index (pts.nrow ());

    for (long int i = 0; i < static_cast <long int> (pts.nrow ()); i++) // Rcpp::nrow is int!
    {
        double dmin = INFINITE_DOUBLE;
        int jmin = INFINITE_INT;
        for (int j = 0; j < xy.nrow (); j++)
        {
            double dij = (vtx [j] - ptx [i]) * (vtx [j] - ptx [i]) +
                (vty [j] - pty [i]) * (vty [j] - pty [i]);
            if (dij < dmin)
            {
                dmin = dij;
                jmin = j;
            }
        }
        if (jmin == INFINITE_INT)
            Rcpp::Rcout << "---ERROR---" << std::endl;
        index (static_cast <size_t> (i)) = jmin;
    }

    return index;
}
