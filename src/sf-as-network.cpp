#include "sf-as-network.h"

// Haversine great circle distance between two points
double haversine (double x1, double y1, double x2, double y2)
{
    double xd = (x2 - x1) * M_PI / 180.0;
    double yd = (y2 - y1) * M_PI / 180.0;
    double d = sin (yd / 2.0) * sin (yd / 2.0) + cos (y2 * M_PI / 180.0) *
        cos (y1 * M_PI / 180.0) * sin (xd / 2.0) * sin (xd / 2.0);
    d = 2.0 * 3671.0 * asin (sqrt (d));
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
    Rcpp::StringVector hw = pr [1];
    Rcpp::NumericVector val = pr [2];
    for (int i = 0; i != hw.size (); i ++)
        profile.insert (std::make_pair (std::string (hw [i]), val [i]));

    Rcpp::CharacterVector nms = sf_lines.attr ("names");
    if (nms [nms.size () - 1] != "geometry")
        throw std::runtime_error ("sf_lines have no geometry component");
    int one_way_index = -1;
    int one_way_bicycle_index = -1;
    int highway_index = -1;
    for (unsigned int i = 0; i < nms.size (); i++)
    {
        if (nms [i] == "oneway")
            one_way_index = static_cast <int> (i);
        if (nms [i] == "oneway.bicycle")
            one_way_bicycle_index = static_cast <int> (i);
        if (nms [i] == "highway")
            highway_index = static_cast <int> (i);
    }
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
            for (unsigned int i = 0; i != ow.size (); ++ i)
                if (ow [i] == "NA" && owb [i] != "NA")
                    ow [i] = owb [i];
        } else if (owb.size () > ow.size ())
            ow = owb;
    }

    Rcpp::List geoms = sf_lines [nms.size () - 1];
    std::vector <std::string> way_names = geoms.attr ("names");
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
        if (ngeoms < ow.size ())
        {
            if (!(ow [ngeoms] == "yes" || ow [ngeoms] == "-1"))
            {
                nrows += rows;
                isOneWay [ngeoms] = true;
            }
        }
        ngeoms ++;
    }

    Rcpp::NumericMatrix nmat = Rcpp::NumericMatrix (Rcpp::Dimension (nrows, 6));
    Rcpp::CharacterMatrix idmat =
        Rcpp::CharacterMatrix (Rcpp::Dimension (nrows, 4));

    nrows = 0;
    ngeoms = 0;
    int fake_id = 0;
    for (auto g = geoms.begin (); g != geoms.end (); ++ g)
    {
        Rcpp::checkUserInterrupt ();
        Rcpp::NumericMatrix gi = (*g);
        std::string hway = std::string (highway [ngeoms]);
        double hw_factor = profile [hway];
        if (hw_factor == 0.0) hw_factor = 1e-5;
        hw_factor = 1.0 / hw_factor;

        Rcpp::List ginames = gi.attr ("dimnames");
        Rcpp::CharacterVector rnms;
        if (ginames.length () > 0)
            rnms = ginames [0];
        else
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
            nmat (nrows, 0) = gi (i-1, 0);
            nmat (nrows, 1) = gi (i-1, 1);
            nmat (nrows, 2) = gi (i, 0);
            nmat (nrows, 3) = gi (i, 1);
            nmat (nrows, 4) = d;
            nmat (nrows, 5) = d * hw_factor;
            idmat (nrows, 0) = rnms (i-1);
            idmat (nrows, 1) = rnms (i);
            idmat (nrows, 2) = hway;
            idmat (nrows, 3) = way_names [ngeoms];
            nrows ++;
            if (isOneWay [ngeoms])
            {
                nmat (nrows, 0) = gi (i, 0);
                nmat (nrows, 1) = gi (i, 1);
                nmat (nrows, 2) = gi (i-1, 0);
                nmat (nrows, 3) = gi (i-1, 1);
                nmat (nrows, 4) = d;
                nmat (nrows, 5) = d * hw_factor;
                idmat (nrows, 0) = rnms (i);
                idmat (nrows, 1) = rnms (i-1);
                idmat (nrows, 2) = hway;
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
    Rcpp::NumericVector xy_x, xy_y, pt_x, pt_y;
    const size_t n;
    RcppParallel::RVector <int> index;

    // constructor
    OnePointIndex (
            const Rcpp::NumericVector xy_x,
            const Rcpp::NumericVector xy_y,
            const Rcpp::NumericVector pt_x,
            const Rcpp::NumericVector pt_y,
            const size_t n,
            Rcpp::IntegerVector index) :
        xy_x (xy_x), xy_y (xy_y), pt_x (pt_x), pt_y (pt_y), n (n),
        index (index)
    {
    }

    // Parallel function operator
    void operator() (std::size_t begin, std::size_t end)
    {
        for (std::size_t i = begin; i < end; i++)
        {
            double dmin = INFINITE_DOUBLE;
            int jmin = INFINITE_INT;
            for (int j = 0; j < n; j++)
            {
                double dij = (xy_x [j] - pt_x [i]) * (xy_x [j] - pt_x [i]) +
                    (xy_y [j] - pt_y [i]) * (xy_y [j] - pt_y [i]);
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

    for (unsigned int i = 0; i < static_cast <unsigned int> (pts.nrow ()); i++) // Rcpp::nrow is int!
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
        index (i) = jmin;
    }

    return index;
}

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

    size_t n = pts.nrow ();

    Rcpp::IntegerVector index (static_cast <int> (n),
            Rcpp::IntegerVector::get_na ());
    // Create parallel worker
    OnePointIndex one_pt_indx (vtx, vty, ptx, pty, n, index);

    RcppParallel::parallelFor (0, n, one_pt_indx);

    return index;
}

// from osmdata/src/get-bbox.cpp
Rcpp::NumericVector rcpp_get_bbox_sf (double xmin, double xmax, double ymin, double ymax)
{
    std::vector <std::string> names;
    names.push_back ("xmin");
    names.push_back ("ymin");
    names.push_back ("xmax");
    names.push_back ("ymax");

    Rcpp::NumericVector bbox (4, NA_REAL);
    bbox (0) = xmin;
    bbox (1) = xmax;
    bbox (2) = ymin;
    bbox (3) = ymax;

    bbox.attr ("names") = names;
    bbox.attr ("class") = "bbox";

    return bbox;
}

//' rcpp_aggregate_to_sf
//'
//' Aggregate a dodgr network data.frame to an sf LINESTRING data.frame
//'
//' @param graph_full Rcpp::DataFrame containing the **full** graph
//' @param graph_contr Rcpp::DataFrame containing the **contracted** graph
//' @param edge_map Rcpp::DataFrame containing the edge map returned from
//' \code{dodgr_contract_graph}
//'
//' @return Rcpp::List object of `sf::LINESTRING` geoms
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::List rcpp_aggregate_to_sf (const Rcpp::DataFrame &graph_full,
        const Rcpp::DataFrame &graph_contr, const Rcpp::DataFrame &edge_map,
        bool displ_progress)
{
    Rcpp::CharacterVector old_edges = edge_map ["edge_old"],
            new_edges = edge_map ["edge_new"],
            contr_edges = graph_contr ["edge_id"],
            idf_r = graph_full ["from_id"],
            idt_r = graph_full ["to_id"];

    Rcpp::NumericVector xf = graph_full ["from_lon"],
            yf = graph_full ["from_lat"],
            xt = graph_full ["to_lon"],
            yt = graph_full ["to_lat"];


    // store vector of names used to name the resultant sf-objects. A
    // corresponding set is also used below to insert non-contracted edges in
    // final result
    std::vector <std::string> new_edge_names;
    std::unordered_set <std::string> new_edge_set;
    new_edge_names.push_back (static_cast <std::string> (new_edges [0]));
    for (size_t i = 1; i < new_edges.size (); i++)
    {
        new_edge_set.emplace (new_edges [i]);
        if (static_cast <std::string> (new_edges [i]) != new_edge_names.back ())
            new_edge_names.push_back (static_cast <std::string> (new_edges [i]));
    }

    const size_t nedges = new_edge_names.size ();
    Rcpp::List edge_sequences (nedges);

    /* The first edge of a sequence won't be necessarily at the start of a
     * sequence, so two maps are made, one from each node to the next, and the
     * other holding the same values in reverse. The latter enables sequences to
     * be traced in reverse by matching end points.
     */
    Progress p (graph_contr.nrow (), displ_progress);
    for (size_t i = 0; i < graph_contr.nrow (); i++)
    {
        Rcpp::checkUserInterrupt ();
        std::map <std::string, std::string> idmap, idmap_rev;
        for (size_t j = 0; j < edge_map.nrow (); j++)
        {
            if (new_edges [j] == contr_edges [i])
            {
                // old_edges is 1-indexed!
                size_t the_edge = atoi (old_edges [j]) - 1;
                if (the_edge > idf_r.size ())
                    Rcpp::stop ("the_edge > size of fulll graph");

                std::string node_f = static_cast <std::string> (idf_r [the_edge]),
                    node_t = static_cast <std::string> (idt_r [the_edge]);
                idmap.emplace (node_f, node_t);
                idmap_rev.emplace (node_t, node_f);
            }
        }
        size_t nnodes = idmap.size ();

        /*
         * idmap then has a map between from and to IDs of original edges for
         * the contracted edge. These are not necessarily in any sequential
         * order, and so the following code constructs longest sequences of edge
         * IDs. If a "to" ID either already exists in idmap, or cannot be found,
         * then that segment is dumped as a list component of edge_ids.
         */

        if (idmap.size () > 0)
        {
            // Find the front of the sequence of which the first idmap pair is
            // part by stepping backwards through idmap_rev
            std::string front_node = idmap.begin ()->first;
            std::string back_node = idmap.begin ()->second;
            while (idmap_rev.find (front_node) != idmap_rev.end ())
                front_node = idmap_rev.find (front_node)->second;
            back_node = idmap.find (front_node)->second;

            std::vector <std::string> id (nnodes + 1);
            id [0] = front_node;
            id [1] = back_node;
            size_t count = 2;
            idmap.erase (front_node);

            while (idmap.size () > 0)
            {
                std::map <std::string, std::string>::iterator it =
                    idmap.find (back_node);
                back_node = it->second;
                id [count++] = back_node;
                idmap.erase (it);
            }

            edge_sequences [i] = id;
        }
        if (displ_progress)
            p.increment ();
    }

    // Then append the non-contracted edges that are in the contracted graph
    Rcpp::CharacterVector idf_r_c = graph_contr ["from_id"],
            idt_r_c = graph_contr ["to_id"];
    size_t edge_count = 0;
    for (size_t i = 0; i < graph_contr.nrow (); i++)
    {
        if (new_edge_set.find (static_cast <std::string> (contr_edges [i])) ==
                new_edge_set.end ())
        {
            edge_count++;
        }
    }
    Rcpp::List edge_sequences_new (edge_count);
    std::vector <std::string> old_edge_names (edge_count);
    edge_count = 0;
    for (size_t i = 0; i < graph_contr.nrow (); i++)
    {
        if (new_edge_set.find (static_cast <std::string> (contr_edges [i])) ==
                new_edge_set.end ())
        {
            old_edge_names [edge_count] = contr_edges [i];
            Rcpp::CharacterVector idvec (2);
            idvec [0] = idf_r_c [i];
            idvec [1] = idt_r_c [i];
            edge_sequences_new [edge_count++] = idvec;
        }
    }
    // Then just join the two edge_sequence Lists together, along with vectors
    // of edge names
    const size_t total_edges = edge_sequences.size () +
        edge_sequences_new.size ();
    Rcpp::List edge_sequences_all (total_edges);
    std::vector <std::string> all_edge_names (total_edges);
    for (size_t i = 0; i < edge_sequences.size (); i++)
    {
        all_edge_names [i] = new_edge_names [i];
        edge_sequences_all [i] = edge_sequences [i];
    }
    for (size_t i = 0; i < edge_sequences_new.size (); i++)
    {
        all_edge_names [edge_sequences.size () + i] = old_edge_names [i];
        edge_sequences_all [edge_sequences.size () + i] = edge_sequences_new [i];
    }

    // Finally the list of edge IDs just has to be converted to corresponding
    // list of coordinate matrices. These first require maps of from and to IDs
    // to integer indices into the graph_full columns
    std::unordered_map <std::string, size_t> edge_num_map;
    for (size_t i = 0; i < idf_r.size (); i++)
    {
        std::string idft = static_cast <std::string> (idf_r [i]) + "-" +
            static_cast <std::string> (idt_r [i]);
        edge_num_map.emplace (idft, i);
    }
    
    Rcpp::List res (total_edges);
    double xmin = INFINITE_DOUBLE, xmax = -INFINITE_DOUBLE,
           ymin = INFINITE_DOUBLE, ymax = -INFINITE_DOUBLE;
    for (size_t i = 0; i < total_edges; i++)
    {
        Rcpp::CharacterVector idvec = edge_sequences_all [i];
        Rcpp::NumericVector x (idvec.size ()), y (idvec.size ());
        // Fill first edge
        std::string id0 = static_cast <std::string> (idvec [0]),
            id1 = static_cast <std::string> (idvec [1]);
        std::string id01 = id0 + "-" + id1;
        size_t j = edge_num_map.find (id01)->second;
        x [0] = xf [j];
        y [0] = yf [j];
        x [1] = xt [j];
        y [1] = yt [j];
        // Then just remaining "to" values
        idvec.erase (0);
        idvec.erase (0);
        size_t count = 2;
        while (idvec.size () > 0)
        {
            id0 = id1;
            id1 = static_cast <std::string> (idvec [0]);
            id01 = id0 + "-" + id1;
            j = edge_num_map.find (id01)->second;
            x [count] = xt [j];
            y [count] = yt [j];
            count++;
            idvec.erase (0);
        }
        Rcpp::NumericMatrix mat (x.size (), 2);
        mat (Rcpp::_, 0) = x;
        mat (Rcpp::_, 1) = y;
        xmin = std::min (xmin, static_cast <double> (Rcpp::min (x)));
        xmax = std::max (xmax, static_cast <double> (Rcpp::max (x)));
        ymin = std::min (ymin, static_cast <double> (Rcpp::min (y)));
        ymax = std::max (ymax, static_cast <double> (Rcpp::max (y)));

        // Then some sf necessities:
        mat.attr ("class") = 
            Rcpp::CharacterVector::create ("XY", "LINESTRING", "sfg");
        res (i) = mat;
    }

    res.attr ("names") = all_edge_names;
    res.attr ("n_empty") = 0;
    res.attr ("class") = Rcpp::CharacterVector::create ("sfc_LINESTRING", "sfc");
    res.attr ("precision") = 0.0;
    res.attr ("bbox") = rcpp_get_bbox_sf (xmin, ymin, xmax, ymax);
    Rcpp::List crs = Rcpp::List::create ((int) 4326, osm_p4s);
    crs.attr ("names") = Rcpp::CharacterVector::create ("epsg", "proj4string");
    crs.attr ("class") = "crs";
    res.attr ("crs") = crs;

    return res;
}
