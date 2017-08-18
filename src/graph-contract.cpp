#include "graph.h"


// See docs/graph-contraction for explanation of the following code and
// associated vertex and edge maps.
void contract_graph (vertex_map_t &vertex_map, edge_map_t &edge_map,
        vert2edge_map_t &vert2edge_map)
{
    std::unordered_set <vertex_id_t> verts;
    for (auto v: vertex_map)
        verts.insert (v.first);

    int max_edge_id = 0;
    for (auto e: edge_map)
        if (e.second.getID () > max_edge_id)
            max_edge_id = e.second.getID ();
    max_edge_id++;

    std::vector<unsigned int> newe;

    while (verts.size () > 0)
    {
        std::unordered_set <vertex_id_t>::iterator vid = verts.begin ();
        vertex_id_t vtx_id = vertex_map.find (*vid)->first;
        vertex_t vtx = vertex_map.find (*vid)->second;
        std::set <unsigned int> edges = vert2edge_map [vtx_id];
        std::map <unsigned int, bool> edges_done;
        for (auto e: edges)
            edges_done.emplace (e, false);

        newe.clear ();
        newe.push_back (max_edge_id++);

        if ((vtx.is_intermediate_single () || vtx.is_intermediate_double ()) &&
                (edges.size () == 2 || edges.size () == 4))
        {
            if (edges.size () == 4) // is_intermediate_double as well!
                newe.push_back (max_edge_id++);

            // remove intervening vertex:
            auto nbs = vtx.get_all_neighbours (); // unordered_set <vertex_id_t>
            std::vector <vertex_id_t> two_nbs;
            for (vertex_id_t nb: nbs)
                two_nbs.push_back (nb); // size is always 2

            vertex_t vt0 = vertex_map [two_nbs [0]];
            vertex_t vt1 = vertex_map [two_nbs [1]];
            // Note that replace neighbour includes bi-directional replacement,
            // so this works for intermediate_double() too
            vt0.replace_neighbour (vtx_id, two_nbs [1]);
            vt1.replace_neighbour (vtx_id, two_nbs [0]);
            vertex_map [two_nbs [0]] = vt0;
            vertex_map [two_nbs [1]] = vt1;

            // construct new edge and remove old ones
            for (unsigned int ne: newe)
            {
                float d = 0.0, w = 0.0;
                vertex_id_t vt_from, vt_to;

                // insert "in" edge first
                std::set <unsigned int> old_edges;

                for (unsigned int e: edges)
                {
                    if (!edges_done [e])
                    {
                        edge_t ei = edge_map.find (e)->second;
                        vt_from = ei.get_from_vertex ();
                        if (vt_from == two_nbs [0] || vt_from == two_nbs [1])
                        {
                            d += ei.dist;
                            w += ei.weight;
                            if (ei.get_old_edges().size() > 0)
                                old_edges = ei.get_old_edges();
                            else
                                old_edges.insert (e);

                            erase_from_edge_map (vert2edge_map, vtx_id, e);
                            add_to_edge_map (vert2edge_map, vtx_id, ne);
                            erase_from_edge_map (vert2edge_map, vt_from, e);
                            add_to_edge_map (vert2edge_map, vt_from, ne);

                            vertex_t vt = vertex_map [vt_from];
                            vt.replace_neighbour (vtx_id, vt_to);
                            vertex_map [vt_from] = vt;

                            edges_done [e] = true;
                            break;
                        }
                    }
                }

                // then "out" edge
                for (unsigned int e: edges)
                {
                    if (!edges_done [e])
                    {
                        edge_t ei = edge_map.find (e)->second;
                        vt_to = ei.get_to_vertex ();
                        if (vt_to == two_nbs [0] || vt_to == two_nbs [1])
                        {
                            d += ei.dist;
                            w += ei.weight;
                            if (ei.get_old_edges().size() > 0)
                            {
                                std::set <unsigned int> tempe = ei.get_old_edges();
                                for (auto et: tempe)
                                    old_edges.insert (et);
                            } else
                                old_edges.insert (e);

                            erase_from_edge_map (vert2edge_map, vtx_id, e);
                            add_to_edge_map (vert2edge_map, vtx_id, ne);
                            erase_from_edge_map (vert2edge_map, vt_to, e);
                            add_to_edge_map (vert2edge_map, vt_to, ne);

                            vertex_t vt = vertex_map [vt_to];
                            vt.replace_neighbour (vtx_id, vt_from);
                            vertex_map [vt_to] = vt;

                            edges_done [e] = true;
                            break;
                        }
                    }
                }
                vert2edge_map[vtx_id].insert(ne);
                edge_t new_edge = edge_t (vt_from, vt_to, d, w, ne, old_edges);
                edge_map.emplace (ne, new_edge);
            }
            for (unsigned int e: edges)
            {
                if (edge_map.find (e) != edge_map.end ())
                    edge_map.erase (e);
            }
        }
        verts.erase (vtx_id);
    }
}

//' rcpp_contract_graph
//'
//' Removes nodes and edges from a graph that are not needed for routing
//'
//' @param graph graph to be processed
//' @param quiet If TRUE, display progress
//'
//' @return \code{Rcpp::List} containing one \code{data.frame} with the compact
//' graph, one \code{data.frame} with the original graph and one
//' \code{data.frame} containing information about the relating edge ids of the
//' original and compact graph.
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::List rcpp_contract_graph (Rcpp::DataFrame graph, bool quiet)
{
    vertex_map_t vertices;
    edge_map_t edge_map;
    std::unordered_map <vertex_id_t, int> components;
    vert2edge_map_t vert2edge_map;

    if (!quiet)
    {
        Rcpp::Rcout << "Constructing graph ... ";
        Rcpp::Rcout.flush ();
    }
    graph_from_df (graph, vertices, edge_map, vert2edge_map);
    if (!quiet)
    {
        Rcpp::Rcout << std::endl << "Determining connected components ... ";
        Rcpp::Rcout.flush ();
    }
    int largest_component = identify_graph_components (vertices, components);
    largest_component++; // suppress unused variable warning

    if (!quiet)
    {
        Rcpp::Rcout << std::endl << "Removing intermediate nodes ... ";
        Rcpp::Rcout.flush ();
    }
    vertex_map_t vertices2 = vertices;
    edge_map_t edge_map2 = edge_map;

    contract_graph (vertices2, edge_map2, vert2edge_map);

    if (!quiet)
    {
        Rcpp::Rcout << std::endl << "Mapping compact to original graph ... ";
        Rcpp::Rcout.flush ();
    }
    int nedges = edge_map2.size ();

    // These vectors are all for the contracted graph:
    Rcpp::StringVector from_vec (nedges), to_vec (nedges),
        highway_vec (nedges);
    Rcpp::NumericVector dist_vec (nedges), weight_vec (nedges),
        edgeid_vec (nedges);

    int map_size = 0; // size of edge map contracted -> original
    int edge_count = 0;
    for (auto e = edge_map2.begin (); e != edge_map2.end (); ++e)
    {
        vertex_id_t from = e->second.get_from_vertex ();
        vertex_id_t to = e->second.get_to_vertex ();
        vertex_t from_vtx = vertices2.at (from);
        vertex_t to_vtx = vertices2.at (to);

        from_vec (edge_count) = from;
        to_vec (edge_count) = to;
        dist_vec (edge_count) = e->second.dist;
        weight_vec (edge_count) = e->second.weight;
        edgeid_vec (edge_count) = e->second.getID ();

        edge_count++;

        map_size += e->second.get_old_edges ().size ();
    }

    Rcpp::NumericVector edge_id_orig (map_size), edge_id_comp (map_size);
    int pos = 0;
    for (auto e = edge_map2.begin (); e != edge_map2.end (); ++e)
    {
        int eid = e->second.getID ();
        std::set <unsigned int> edges = e->second.get_old_edges ();
        for (auto ei: edges)
        {
            edge_id_comp (pos) = eid;
            edge_id_orig (pos++) = ei;
        }
    }

    Rcpp::DataFrame compact = Rcpp::DataFrame::create (
            Rcpp::Named ("from_id") = from_vec,
            Rcpp::Named ("to_id") = to_vec,
            Rcpp::Named ("edge_id") = edgeid_vec,
            Rcpp::Named ("d") = dist_vec,
            Rcpp::Named ("d_weighted") = weight_vec);

    Rcpp::DataFrame rel = Rcpp::DataFrame::create (
            Rcpp::Named ("id_compact") = edge_id_comp,
            Rcpp::Named ("id_original") = edge_id_orig);

    if (!quiet)
        Rcpp::Rcout << std::endl;

    return Rcpp::List::create (
            Rcpp::Named ("compact") = compact,
            Rcpp::Named ("original") = graph,
            Rcpp::Named ("map") = rel);
}


//' rcpp_insert_vertices
//'
//' Insert routing vertices in compact graph
//'
//' @param fullgraph graph to be processed
//' @param compactgraph graph to be processed
//' @param pts_to_insert Index into graph of those points closest to desired
//' routing points. These are to be kept in the compact graph.
//' @return \code{Rcpp::List} containing one \code{data.frame} with the compact
//' graph, one \code{data.frame} with the original graph and one
//' \code{data.frame} containing information about the relating edge ids of the
//' original and compact graph.
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::List rcpp_insert_vertices (Rcpp::DataFrame fullgraph,
        Rcpp::DataFrame compactgraph, std::vector <int> pts_to_insert)
{
    vertex_map_t vertices;
    edge_map_t edge_map;
    std::unordered_map <vertex_id_t, int> components;
    vert2edge_map_t vert2edge_map;

    graph_from_df (fullgraph, vertices, edge_map, vert2edge_map);

    return Rcpp::List::create ();
}
