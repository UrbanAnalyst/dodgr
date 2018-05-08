#include "dodgr-to-sf.h"

//' Make unordered_set of all new edge names
//' @noRd
size_t dodgr_sf::make_edge_name_set (std::unordered_set <std::string> &new_edge_name_set,
        const Rcpp::CharacterVector &new_edges)
{
    new_edge_name_set.clear ();
    // Rcpp::CharacterVector requires long int index type
    for (long int i = 0; i < new_edges.size (); i++)
    {
        new_edge_name_set.emplace (static_cast <std::string> (new_edges [i]));
    }
    return new_edge_name_set.size ();
}

//' Make vector of all new edge names
//' @noRd
void dodgr_sf::make_edge_name_vec (const size_t n,
        const Rcpp::CharacterVector &new_edges,
        std::vector <std::string> &new_edge_name_vec)
{
    new_edge_name_vec.clear ();
    new_edge_name_vec.resize (n);
    new_edge_name_vec [0] = static_cast <std::string> (new_edges [0]);
    size_t count = 0;
    for (long int i = 1; i < new_edges.size (); i++)
    {
        std::string new_edge_i = static_cast <std::string> (new_edges [i]);
        if (new_edge_i != new_edge_name_vec [count])
            new_edge_name_vec [++count] = new_edge_i;
    }
}


//' Get numbers of old edges corresponding to each contracted edge
//' @noRd
size_t dodgr_sf::get_edgevec_sizes (const size_t nedges,
        const Rcpp::CharacterVector &new_edges,
        std::vector <size_t> &edgevec_sizes)
{
    edgevec_sizes.clear ();
    edgevec_sizes.resize (nedges);
    size_t count = 1, edgenum = 0;
    for (long int i = 1; i < new_edges.size (); i++)
    {
        if (new_edges [i] == new_edges [i - 1])
            count++;
        else
        {
            edgevec_sizes [edgenum++] = count;
            count = 1;
        }
    }
    // And last count too:
    edgevec_sizes [edgenum++] = count;
    return edgenum;
}

//' Collect sets of all from and to vertices for each set of edges corresponding
//' to each contracted edge. These sets aren't in any particular order, but the
//' two sets may be used to match the from and to vertices. These sets are then
//' arranged into sequences in the subsequent function,
//' \code{order_vert_sequences}.
//' @noRd
void dodgr_sf::get_edge_to_vert_maps (const std::vector <size_t> &edgevec_sizes,
        const Rcpp::DataFrame &graph_full,
        const Rcpp::CharacterVector &old_edges,
        const Rcpp::CharacterVector &new_edges,
        const std::vector <std::string> &new_edge_names,
        std::unordered_map <std::string,
                            std::vector <std::string> > &full_from_edge_map,
        std::unordered_map <std::string,
                            std::vector <std::string> > &full_to_edge_map)
{
    Rcpp::CharacterVector idf_r = graph_full ["from_id"],
            idt_r = graph_full ["to_id"];

    size_t count = 1, edgenum = 0;
    std::vector <std::string> from_node, to_node;
    from_node.resize (edgevec_sizes [edgenum]);
    to_node.resize (edgevec_sizes [edgenum]);
    // old_edges is 1-indexed!
    long int the_edge = static_cast <long int> (atoi (old_edges [0])) - 1;
    from_node [0] = static_cast <std::string> (idf_r [the_edge]);
    to_node [0] = static_cast <std::string> (idt_r [the_edge]);
    for (long int i = 1; i < new_edges.size (); i++)
    {
        if (new_edges [i] != new_edges [i - 1])
        {
            full_from_edge_map.emplace (new_edge_names [edgenum], from_node);
            full_to_edge_map.emplace (new_edge_names [edgenum], to_node);
            from_node.clear ();
            to_node.clear ();
            edgenum++;
            from_node.resize (edgevec_sizes [edgenum]);
            to_node.resize (edgevec_sizes [edgenum]);
            count = 0;
        }
        the_edge = atoi (old_edges [i]) - 1; // it's 1-indexed!
        from_node [count] = static_cast <std::string> (idf_r [the_edge]);
        to_node [count++] = static_cast <std::string> (idt_r [the_edge]);
    }
    full_from_edge_map.emplace (new_edge_names [edgenum], from_node);
    full_to_edge_map.emplace (new_edge_names [edgenum], to_node);

    from_node.clear ();
    to_node.clear ();
}

void dodgr_sf::order_vert_sequences (Rcpp::List &edge_sequences,
        std::vector <std::string> &new_edge_names,
        std::unordered_map <std::string,
                            std::vector <std::string> > &full_from_edge_map,
        std::unordered_map <std::string,
                            std::vector <std::string> > &full_to_edge_map)
{
    const size_t nedges = static_cast <size_t> (edge_sequences.size ());
    for (size_t i = 0; i < nedges; i++)
    {
        Rcpp::checkUserInterrupt ();
        std::map <std::string, std::string> idmap, idmap_rev;
        std::string this_edge = new_edge_names [i];
        std::vector <std::string> from_edges = full_from_edge_map [this_edge],
            to_edges = full_to_edge_map [this_edge];
        for (size_t j = 0; j < from_edges.size (); j++)
        {
            idmap.emplace (from_edges [j], to_edges [j]);
            idmap_rev.emplace (to_edges [j], from_edges [j]);
        }
        size_t nnodes = idmap.size ();

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
        idmap_rev.clear ();

        edge_sequences [static_cast <long int> (i)] = id;
    }
}

// Count number of non-contracted edges in the contracted graph (non-contracted
// because they can not be simplified).
size_t dodgr_sf::count_non_contracted_edges (const Rcpp::CharacterVector &contr_edges,
        std::unordered_set <std::string> &new_edge_name_set)
{
    size_t edge_count = 0;
    for (long int i = 0; i < contr_edges.size (); i++)
    {
        if (new_edge_name_set.find (static_cast <std::string>
                    (contr_edges [i])) == new_edge_name_set.end ())
        {
            edge_count++;
        }
    }
    return edge_count;
}

// Append non-contracted edges to contracted edges in all edge-to-vertex maps
void dodgr_sf::append_nc_edges (const size_t nc_edge_count, 
        const Rcpp::DataFrame &graph_contr,
        std::unordered_set <std::string> &new_edge_name_set,
        std::vector <std::string> &new_edge_name_vec,
        const Rcpp::List &edge_sequences_contr,
        std::vector <std::string> &all_edge_names,
        Rcpp::List &edge_sequences_all)
{
    Rcpp::List edge_sequences_new (nc_edge_count);
    std::vector <std::string> old_edge_names (nc_edge_count);
    size_t count = 0;
    Rcpp::CharacterVector idf_r_c = graph_contr ["from_id"],
            idt_r_c = graph_contr ["to_id"],
            contr_edges = graph_contr ["edge_id"];
    for (long int i = 0; i < graph_contr.nrow (); i++)
    {
        if (new_edge_name_set.find (static_cast <std::string>
                    (contr_edges [i])) == new_edge_name_set.end ())
        {
            old_edge_names [count] = contr_edges [i];
            Rcpp::CharacterVector idvec (2);
            idvec [0] = idf_r_c [i];
            idvec [1] = idt_r_c [i];
            edge_sequences_new [static_cast <long int> (count++)] = idvec;
        }
    }
    // Then just join the two edge_sequence Lists together, along with vectors
    // of edge names
    const size_t total_edges = static_cast <size_t> (edge_sequences_contr.size ()) +
        nc_edge_count;
    all_edge_names.resize (total_edges);
    // These two have different indexing types:
    for (size_t i = 0; i < static_cast <size_t> (edge_sequences_contr.size ()); i++)
        all_edge_names [i] = new_edge_name_vec [i];
    for (long int i = 0; i < edge_sequences_contr.size (); i++)
        edge_sequences_all [i] = edge_sequences_contr [i];
    for (size_t i = 0; i < static_cast <size_t> (edge_sequences_new.size ()); i++)
        all_edge_names [static_cast <size_t> (edge_sequences_contr.size ()) + i] =
            old_edge_names [i];
    for (long int i = 0; i < edge_sequences_new.size (); i++)
        edge_sequences_all [edge_sequences_contr.size () + i] = edge_sequences_new [i];
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

// Convert the Rcpp::List of vertex names to corresponding coordinates,
// construct matrices of sequential coordinates for each contracted edge, and
// convert these to sf geometry objects
void dodgr_sf::xy_to_sf (const Rcpp::DataFrame &graph_full,
        const Rcpp::List &edge_sequences, 
        const std::vector <std::string> &all_edge_names,
        Rcpp::List &res)
{
    const size_t total_edges = static_cast <size_t> (res.size ());

    Rcpp::CharacterVector idf_r = graph_full ["from_id"],
            idt_r = graph_full ["to_id"];

    Rcpp::NumericVector xf = graph_full ["from_lon"],
            yf = graph_full ["from_lat"],
            xt = graph_full ["to_lon"],
            yt = graph_full ["to_lat"];

    std::unordered_map <std::string, size_t> edge_num_map;
    for (long int i = 0; i < idf_r.size (); i++)
    {
        std::string idft = static_cast <std::string> (idf_r [i]) + "-" +
            static_cast <std::string> (idt_r [i]);
        edge_num_map.emplace (idft, i);
    }

    double xmin = INFINITE_DOUBLE, xmax = -INFINITE_DOUBLE,
           ymin = INFINITE_DOUBLE, ymax = -INFINITE_DOUBLE;
    for (size_t i = 0; i < total_edges; i++)
    {
        Rcpp::CharacterVector idvec = edge_sequences [static_cast <long int> (i)];
        Rcpp::NumericVector x (idvec.size ()), y (idvec.size ());
        // Fill first edge
        std::string id0 = static_cast <std::string> (idvec [0]),
            id1 = static_cast <std::string> (idvec [1]);
        std::string id01 = id0 + "-" + id1;
        long int j = static_cast <long int> (edge_num_map.find (id01)->second);
        x [0] = xf [j];
        y [0] = yf [j];
        x [1] = xt [j];
        y [1] = yt [j];
        // Then just remaining "to" values
        idvec.erase (0);
        idvec.erase (0);
        long int count = 2;
        while (idvec.size () > 0)
        {
            id0 = id1;
            id1 = static_cast <std::string> (idvec [0]);
            id01 = id0 + "-" + id1;
            j = static_cast <long int> (edge_num_map.find (id01)->second);
            x [count] = xt [j];
            y [count] = yt [j];
            count++;
            idvec.erase (0);
        }
        Rcpp::NumericMatrix mat (static_cast <const int> (x.size ()), 2);
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
    Rcpp::List crs = Rcpp::List::create (static_cast <int> (4326), osm_p4s);
    crs.attr ("names") = Rcpp::CharacterVector::create ("epsg", "proj4string");
    crs.attr ("class") = "crs";
    res.attr ("crs") = crs;
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
        const Rcpp::DataFrame &graph_contr, const Rcpp::DataFrame &edge_map)
{
    Rcpp::CharacterVector old_edges = edge_map ["edge_old"],
            new_edges = edge_map ["edge_new"],
            contr_edges = graph_contr ["edge_id"];

    // store vector of names used to name the resultant sf-objects. A
    // corresponding set is also used below to insert non-contracted edges in
    // final result
    std::vector <std::string> new_edge_names;
    std::unordered_set <std::string> new_edge_set;
    const size_t nedges = dodgr_sf::make_edge_name_set (new_edge_set, new_edges);
    dodgr_sf::make_edge_name_vec (nedges, new_edges, new_edge_names);

    // Then get numbers of original edges for each contracted edge
    std::vector <size_t> edgevec_sizes;
    size_t check = dodgr_sf::get_edgevec_sizes (nedges, new_edges, edgevec_sizes);
    if (check != nedges)
        Rcpp::stop ("number of new edges in contracted graph not the right size");

    /* Map contracted edge ids onto corresponding pairs of from and to vertices.
     * The first vertex of a sequence won't be necessarily at the start of a
     * sequence, so two maps are made, one from each node to the next, and the
     * other holding the same values in reverse. The latter enables sequences to
     * be traced in reverse by matching end points.
     */
    std::unordered_map <std::string,
        std::vector <std::string> > full_from_edge_map, full_to_edge_map;
    dodgr_sf::get_edge_to_vert_maps (edgevec_sizes, graph_full, old_edges, new_edges,
            new_edge_names, full_from_edge_map, full_to_edge_map);
    edgevec_sizes.clear ();

    Rcpp::List edge_sequences (nedges); // holds final sequences
    // The main work done in this whole file:
    dodgr_sf::order_vert_sequences (edge_sequences, new_edge_names, full_from_edge_map,
            full_to_edge_map);
    full_from_edge_map.clear ();
    full_to_edge_map.clear ();

    // Then append the non-contracted edges that are in the contracted graph
    size_t nc_edge_count = dodgr_sf::count_non_contracted_edges (contr_edges,
            new_edge_set);
    std::vector <std::string> all_edge_names;
    const size_t total_edges = nedges + nc_edge_count;
    Rcpp::List edge_sequences_all (total_edges);
    dodgr_sf::append_nc_edges (nc_edge_count, graph_contr,
            new_edge_set, new_edge_names, edge_sequences,
            all_edge_names, edge_sequences_all);
    new_edge_names.clear ();
    new_edge_set.clear ();

    Rcpp::List res (total_edges);
    dodgr_sf::xy_to_sf (graph_full, edge_sequences_all, all_edge_names, res);

    return res;
}
