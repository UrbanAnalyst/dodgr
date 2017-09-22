#include "graph.h"

edge_id_t new_edge_id (edge_map_t &edge_map, std::mt19937 &rng)
{
    const int range = 1e8;
    std::uniform_int_distribution <unsigned int> unif (0, range);

    edge_id_t new_id = edge_map.begin()->first;
    while (edge_map.find (new_id) != edge_map.end ())
    {
        new_id = std::to_string (unif (rng));
    }
    return new_id;
}


// See docs/graph-contraction for explanation of the following code and
// associated vertex and edge maps.
void contract_graph (vertex_map_t &vertex_map, edge_map_t &edge_map,
        vert2edge_map_t &vert2edge_map)
{
    std::unordered_set <vertex_id_t> verts;
    for (auto v: vertex_map)
        verts.insert (v.first);

    std::vector<edge_id_t> newe;

    // Random generator for new_edge_id
    std::random_device rd;
    std::mt19937 rng (rd()); // mersenne twister

    while (verts.size () > 0)
    {
        std::unordered_set <vertex_id_t>::iterator vid = verts.begin ();
        vertex_id_t vtx_id = vertex_map.find (*vid)->first;
        vertex_t vtx = vertex_map.find (*vid)->second;
        std::unordered_set <edge_id_t> edges = vert2edge_map [vtx_id];
        std::unordered_map <edge_id_t, bool> edges_done;
        for (auto e: edges)
            edges_done.emplace (e, false);

        newe.clear ();
        newe.push_back (new_edge_id (edge_map, rng));

        if ((vtx.is_intermediate_single () || vtx.is_intermediate_double ()) &&
                (edges.size () == 2 || edges.size () == 4))
        {
            if (edges.size () == 4) // is_intermediate_double as well!
                newe.push_back (new_edge_id (edge_map, rng));

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
            for (edge_id_t ne: newe)
            {
                float d = 0.0, w = 0.0;
                vertex_id_t vt_from, vt_to;

                // insert "in" edge first
                std::set <edge_id_t> old_edges;

                for (edge_id_t e: edges)
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
                for (edge_id_t e: edges)
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
                                std::set <edge_id_t> tempe = ei.get_old_edges();
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
            for (edge_id_t e: edges)
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
//' @return \code{Rcpp::List} containing one \code{data.frame} with the
//' contracted graph, one \code{data.frame} with the original graph and one
//' \code{data.frame} containing information about the relating edge ids of the
//' original and contracted graph.
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::List rcpp_contract_graph (Rcpp::DataFrame graph, bool quiet)
{
    vertex_map_t vertices;
    edge_map_t edge_map;
    std::unordered_map <vertex_id_t, unsigned int> components;
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
        Rcpp::Rcout << std::endl << "Mapping contracted to original graph ... ";
        Rcpp::Rcout.flush ();
    }
    int nedges = edge_map2.size ();

    // These vectors are all for the contracted graph:
    Rcpp::StringVector from_vec (nedges), to_vec (nedges),
        edgeid_vec (nedges), highway_vec (nedges);
    Rcpp::NumericVector dist_vec (nedges), weight_vec (nedges);

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

    Rcpp::StringVector edge_id_orig (map_size), edge_id_comp (map_size);
    int pos = 0;
    for (auto e = edge_map2.begin (); e != edge_map2.end (); ++e)
    {
        edge_id_t eid = e->second.getID ();
        std::set <edge_id_t> edges = e->second.get_old_edges ();
        for (auto ei: edges)
        {
            edge_id_comp (pos) = eid;
            edge_id_orig (pos++) = ei;
        }
    }

    Rcpp::DataFrame contracted = Rcpp::DataFrame::create (
            Rcpp::Named ("edge_id") = edgeid_vec,
            Rcpp::Named ("from") = from_vec,
            Rcpp::Named ("to") = to_vec,
            Rcpp::Named ("d") = dist_vec,
            Rcpp::Named ("w") = weight_vec,
            Rcpp::_["stringsAsFactors"] = false);

    Rcpp::DataFrame map = Rcpp::DataFrame::create (
            Rcpp::Named ("id_contracted") = edge_id_comp,
            Rcpp::Named ("id_original") = edge_id_orig);

    if (!quiet)
        Rcpp::Rcout << std::endl;

    return Rcpp::List::create (
            Rcpp::Named ("contracted") = contracted,
            Rcpp::Named ("map") = map);
}


//' rcpp_insert_vertices
//'
//' Insert routing vertices in contracted graph
//'
//' @param fullgraph graph to be processed
//' @param contracted graph to be processed
//' @param map Map of old-to-new vertices returned from rcpp_contract_graph
//' @param pts_to_insert Index into graph of those points closest to desired
//' routing points. These are to be re-inserted in the contracted graph.
//' @return \code{Rcpp::List} containing one \code{data.frame} with the
//' contracted graph, one \code{data.frame} with the original graph and one
//' \code{data.frame} containing information about the relating edge ids of the
//' original and contracted graph.
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::List rcpp_insert_vertices (Rcpp::DataFrame fullgraph,
        Rcpp::DataFrame contracted, Rcpp::DataFrame map,
        std::vector <std::string> verts_to_insert)
{
    // verts_to_insert is actually vert_id_t, but that can't be exported to Rcpp
    vertex_map_t vertices;
    edge_map_t edge_map;
    std::unordered_map <vertex_id_t, int> components;
    vert2edge_map_t vert2edge_map;

    graph_from_df (fullgraph, vertices, edge_map, vert2edge_map);

    Rcpp::Rcout << "There are " << verts_to_insert.size () <<
        " vertices to insert" << std::endl;

    return Rcpp::List::create ();
}
