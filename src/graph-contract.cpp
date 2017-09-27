#include "graph.h"

edge_id_t get_new_edge_id (edge_map_t &edge_map, std::mt19937 &rng)
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

//' get_to_from
//'
//' Get one pair of two and from edges and vertices. Main task is to make sure
//' that bi-directed edges ("intermediate_double") correctly return the
//' **different** values of from and to vertices and edges.
//'
//' @noRd
void get_to_from (const edge_map_t &edge_map,
        const std::unordered_set <edge_id_t> &edges,
        const std::vector <vertex_id_t> &two_nbs,
        vertex_id_t &vt_from, vertex_id_t &vt_to,
        edge_id_t &edge_from_id, edge_id_t &edge_to_id)
{
    for (edge_id_t edge_id: edges)
    {
        edge_t edge = edge_map.find (edge_id)->second;
        if (edge_from_id == "")
        {
            if (edge.get_from_vertex () == two_nbs [0] &&
                    (edge_to_id == "" || vt_to != two_nbs [0]))
            {
                edge_from_id = edge_id;
                vt_from = two_nbs [0];
            } else if (edge.get_from_vertex () == two_nbs [1] &&
                    (edge_to_id == "" || vt_to != two_nbs [1]))
            {
                edge_from_id = edge_id;
                vt_from = two_nbs [1];
            }
        }
        if (edge_to_id == "")
        {
            if (edge.get_to_vertex () == two_nbs [0] &&
                    (edge_from_id == "" || vt_from != two_nbs [0]))
            {
                edge_to_id = edge_id;
                vt_to = two_nbs [0];
            } else if (edge.get_to_vertex () == two_nbs [1] &&
                    (edge_from_id == "" || vt_from != two_nbs [1]))
            {
                edge_to_id = edge_id;
                vt_to = two_nbs [1];
            }
        }
    }
}

void contract_one_edge (vert2edge_map_t &vert2edge_map,
        vertex_map_t &vertex_map, edge_map_t &edge_map,
        const std::unordered_set <edge_id_t> &edgelist,
        const vertex_id_t vtx_id, const vertex_id_t vt_from,
        const vertex_id_t vt_to,
        const edge_id_t edge_from_id, const edge_id_t edge_to_id,
        const edge_id_t new_edge_id)
{
    edge_t edge_from = edge_map.find (edge_from_id)->second,
           edge_to = edge_map.find (edge_to_id)->second;
    float d = edge_from.dist + edge_to.dist,
          w = edge_from.weight + edge_to.weight;

    std::set <edge_id_t> old_edges,
        old_edges_fr = edge_from.get_old_edges (),
        old_edges_to = edge_to.get_old_edges ();
    for (auto e: old_edges_fr)
        old_edges.insert (e);
    for (auto e: old_edges_to)
        old_edges.insert (e);
    // Only insert edge IDs that are in the original list
    if (edgelist.find (edge_from_id) != edgelist.end ())
        old_edges.insert (edge_from_id);
    if (edgelist.find (edge_to_id) != edgelist.end ())
        old_edges.insert (edge_to_id);

    erase_from_v2e_map (vert2edge_map, vtx_id, edge_from_id);
    erase_from_v2e_map (vert2edge_map, vtx_id, edge_to_id);
    erase_from_v2e_map (vert2edge_map, vt_from, edge_from_id);
    erase_from_v2e_map (vert2edge_map, vt_to, edge_to_id);
    add_to_v2e_map (vert2edge_map, vt_from, new_edge_id);
    add_to_v2e_map (vert2edge_map, vt_to, new_edge_id);

    vertex_t vt = vertex_map [vt_from];
    vt.replace_neighbour (vtx_id, vt_from);
    vertex_map [vt_from] = vt;
    vt = vertex_map [vt_to];
    vt.replace_neighbour (vtx_id, vt_to);
    vertex_map [vt_to] = vt;

    edge_map.erase (edge_from_id);
    edge_map.erase (edge_to_id);
    edge_t new_edge = edge_t (vt_from, vt_to, d, w,
            new_edge_id, old_edges);
    edge_map.emplace (new_edge_id, new_edge);
}


// See docs/graph-contraction for explanation of the following code and
// associated vertex and edge maps.
void contract_graph (vertex_map_t &vertex_map, edge_map_t &edge_map,
        vert2edge_map_t &vert2edge_map)
{
    std::unordered_set <vertex_id_t> verts;
    for (auto v: vertex_map)
        verts.insert (v.first);
    std::unordered_set <edge_id_t> edgelist;
    for (auto e: edge_map)
        edgelist.insert (e.first);

    std::vector<edge_id_t> new_edge_ids;

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

        new_edge_ids.clear ();
        new_edge_ids.push_back (get_new_edge_id (edge_map, rng));

        if ((vtx.is_intermediate_single () || vtx.is_intermediate_double ()) &&
                (edges.size () == 2 || edges.size () == 4))
        {
            if (edges.size () == 4) // is_intermediate_double as well!
                new_edge_ids.push_back (get_new_edge_id (edge_map, rng));

            // remove intervening vertex:
            std::unordered_set <vertex_id_t> nbs = vtx.get_all_neighbours ();
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

            // construct new edge(s) and remove old ones. There are 2
            // new_edge_ids only for intermediate double vertices
            // (that is, bi-directional).
            for (edge_id_t new_edge_id: new_edge_ids)
            {
                //float d = 0.0, w = 0.0;
                vertex_id_t vt_from = "", vt_to = "";
                edge_id_t edge_from_id = "", edge_to_id = "";

                // get the from and to edges and vertices
                get_to_from (edge_map, edges, two_nbs,
                        vt_from, vt_to, edge_from_id, edge_to_id);
                edges.erase (edge_from_id);
                edges.erase (edge_to_id);

                contract_one_edge (vert2edge_map, vertex_map, edge_map,
                        edgelist, vtx_id, vt_from, vt_to,
                        edge_from_id, edge_to_id, new_edge_id);
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
    vertex_map_t vertices_contracted = vertices;
    edge_map_t edge_map_contracted = edge_map;

    contract_graph (vertices_contracted, edge_map_contracted, vert2edge_map);

    if (!quiet)
    {
        Rcpp::Rcout << std::endl << "Mapping contracted to original graph ... ";
        Rcpp::Rcout.flush ();
    }
    int nedges = edge_map_contracted.size ();

    // These vectors are all for the contracted graph:
    Rcpp::StringVector from_vec (nedges), to_vec (nedges),
        edgeid_vec (nedges), highway_vec (nedges);
    Rcpp::NumericVector dist_vec (nedges), weight_vec (nedges);

    int map_size = 0; // size of edge map contracted -> original
    int edge_count = 0;
    for (auto e = edge_map_contracted.begin ();
            e != edge_map_contracted.end (); ++e)
    {
        vertex_id_t from = e->second.get_from_vertex ();
        vertex_id_t to = e->second.get_to_vertex ();
        vertex_t from_vtx = vertices_contracted.at (from);
        vertex_t to_vtx = vertices_contracted.at (to);

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
    for (auto e = edge_map_contracted.begin ();
            e != edge_map_contracted.end (); ++e)
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
            Rcpp::Named ("id_original") = edge_id_orig,
            Rcpp::_["stringsAsFactors"] = false);

    if (!quiet)
        Rcpp::Rcout << std::endl;

    return Rcpp::List::create (
            Rcpp::Named ("contracted") = contracted,
            Rcpp::Named ("map") = map);
}


//' rcpp_insert_vertices
//'
//' Re-insert routing vertices in contracted graph, actually done just by
//' returning lists of edge IDs to be re-inserted from full graph, and
//' corresponding contracted edge IDs to be removed from contracted graph.
//'
//' @param fullgraph graph to be processed
//' @param contracted graph to be processed
//' @param map Map of old-to-new vertices returned from rcpp_contract_graph
//' @param pts_to_insert Names of vertices to be re-inserted.
//' @return \code{Rcpp::List} of names of edges to re-insert and edges to remove
//' from contracted graph.
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
    vert2edge_map_t vert2edge_map;

    graph_from_df (fullgraph, vertices, edge_map, vert2edge_map);

    // Get maps of new-to-old and old-to-new edges
    Rcpp::StringVector edge_orig = map ["id_original"];
    Rcpp::StringVector edge_contr = map ["id_contracted"];

    std::unordered_map <edge_id_t, edge_id_t> edge_old2new_map;
    std::unordered_map <edge_id_t, std::set <edge_id_t> > edge_new2old_map;
    for (int i = 0; i < edge_orig.size (); i++)
    {
        edge_id_t eor = Rcpp::as <edge_id_t> (edge_orig (i)),
                  eco = Rcpp::as <edge_id_t> (edge_contr (i));
        edge_old2new_map.emplace (eor, eco);

        std::set <edge_id_t> these_edges;
        if (edge_new2old_map.find (eco) != edge_new2old_map.end ())
        {
            these_edges = edge_new2old_map.find (eco)->second;
            edge_new2old_map.erase (eco);
        }
        these_edges.emplace (eor);
        edge_new2old_map.emplace (eco, these_edges);
    }

    // Convert verts_to_insert to a list of edges_to_insert (according to edge
    // IDs from original graph), and edges_to_erase. Note that verts_to_insert
    // only includes vertices **not** present in contracted graph, meaning they
    // must map onto new (contracted) edge IDs. Note also that vert2edge_map is
    // made from the original graph. Original edge IDs thus have to be mapped
    // onto contracted edge IDs, and these in turn re-mapped onto the full lists
    // of **all** original edges which they have replaced.
    std::unordered_set <edge_id_t> edges_to_insert, edges_to_erase;
    for (auto i: verts_to_insert)
    {
        std::unordered_set <edge_id_t> edges =
            vert2edge_map.find (i)->second;
        for (auto e: edges)
        {
            // Not all edges will have been replaced, so this is necssary:
            if (edge_old2new_map.find (e) != edge_old2new_map.end ())
            {
                edge_id_t new_edge = edge_old2new_map.find (e)->second;
                std::set <edge_id_t> all_old_edges = edge_new2old_map.find (new_edge)->second;
                for (auto ei: all_old_edges)
                    edges_to_insert.emplace (ei);
                edges_to_erase.emplace (new_edge);
            }
        }
    }

    // Then just return lists of edges_to_insert and edges_to_erase
    Rcpp::StringVector edge_insert (edges_to_insert.size ()),
        edge_erase (edges_to_erase.size ());
    unsigned int i = 0;
    for (const auto &e: edges_to_insert)
        edge_insert (i++) = e;
    i = 0;
    for (const auto &e: edges_to_erase)
        edge_erase (i++) = e;

    return Rcpp::List::create (
            Rcpp::Named ("insert") = edge_insert,
            Rcpp::Named ("erase") = edge_erase);
}
