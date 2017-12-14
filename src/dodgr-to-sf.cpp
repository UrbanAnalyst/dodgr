#include "dodgr-to-sf.h"

//' Make unordered_set of all new edge names
//' @noRd
size_t make_edge_name_set (std::unordered_set <std::string> &new_edge_name_set,
        const Rcpp::CharacterVector &new_edges)
{
    new_edge_name_set.clear ();
    for (size_t i = 0; i < new_edges.size (); i++)
    {
        new_edge_name_set.emplace (static_cast <std::string> (new_edges [i]));
    }
    return new_edge_name_set.size ();
}

//' Make vector of all new edge names
//' @noRd
void make_edge_name_vec (const size_t n,
        const Rcpp::CharacterVector &new_edges,
        std::vector <std::string> &new_edge_name_vec)
{
    new_edge_name_vec.clear ();
    new_edge_name_vec.resize (n);
    new_edge_name_vec [0] = static_cast <std::string> (new_edges [0]);
    size_t count = 0;
    for (size_t i = 1; i < new_edges.size (); i++)
    {
        std::string new_edge_i = static_cast <std::string> (new_edges [i]);
        if (new_edge_i != new_edge_name_vec [count])
            new_edge_name_vec [++count] = new_edge_i;
    }
}


//' Get numbers of old edges corresponding to each contracted edge
//' @noRd
size_t get_edgevec_sizes (const size_t nedges,
        const Rcpp::CharacterVector &new_edges,
        std::vector <size_t> &edgevec_sizes)
{
    edgevec_sizes.clear ();
    edgevec_sizes.resize (nedges);
    size_t count = 1, edgenum = 0;
    for (size_t i = 1; i < new_edges.size (); i++)
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

void get_edge_to_vert_maps (const std::vector <size_t> &edgevec_sizes,
        const Rcpp::CharacterVector &idf_r,
        const Rcpp::CharacterVector &idt_r,
        const Rcpp::CharacterVector &old_edges,
        const Rcpp::CharacterVector &new_edges,
        const std::vector <std::string> &new_edge_names,
        std::unordered_map <std::string,
                            std::vector <std::string> > &full_from_edge_map,
        std::unordered_map <std::string,
                            std::vector <std::string> > &full_to_edge_map)
{
    size_t count = 1, edgenum = 0;
    std::vector <std::string> from_node, to_node;
    from_node.resize (edgevec_sizes [edgenum]);
    to_node.resize (edgevec_sizes [edgenum]);
    size_t the_edge = atoi (old_edges [0]) - 1; // it's 1-indexed!
    from_node [0] = static_cast <std::string> (idf_r [the_edge]);
    to_node [0] = static_cast <std::string> (idt_r [the_edge]);
    for (size_t i = 1; i < new_edges.size (); i++)
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
        size_t the_edge = atoi (old_edges [i]) - 1; // it's 1-indexed!
        from_node [count] = static_cast <std::string> (idf_r [the_edge]);
        to_node [count++] = static_cast <std::string> (idt_r [the_edge]);
    }
    full_from_edge_map.emplace (new_edge_names [edgenum], from_node);
    full_to_edge_map.emplace (new_edge_names [edgenum], to_node);

    from_node.clear ();
    to_node.clear ();
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
        const Rcpp::DataFrame &graph_contr, const Rcpp::DataFrame &edge_map)
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
    const size_t nedges = make_edge_name_set (new_edge_set, new_edges);
    make_edge_name_vec (nedges, new_edges, new_edge_names);

    // Then get numbers of original edges for each contracted edge
    std::vector <size_t> edgevec_sizes;
    size_t check = get_edgevec_sizes (nedges, new_edges, edgevec_sizes);
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
    get_edge_to_vert_maps (edgevec_sizes, idf_r, idt_r, old_edges, new_edges,
            new_edge_names, full_from_edge_map, full_to_edge_map);

    Rcpp::List edge_sequences (nedges); // holds final sequences
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
            Rcpp::checkUserInterrupt ();
            std::map <std::string, std::string>::iterator it =
                idmap.find (back_node);
            back_node = it->second;
            id [count++] = back_node;
            idmap.erase (it);
        }
        idmap_rev.clear ();

        edge_sequences [i] = id;
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
