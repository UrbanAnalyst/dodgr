#include "graph.h"


void add_to_edge_map (vert2edge_map_t &vert2edge_map, osm_id_t vid,
        unsigned int eid)
{
    std::set <unsigned int> edge_ids;
    if (vert2edge_map.find (vid) == vert2edge_map.end ())
    {
        edge_ids.emplace (eid);
        vert2edge_map.emplace (vid, edge_ids);
    } else
    {
        edge_ids = vert2edge_map [vid];
        edge_ids.emplace (eid);
        vert2edge_map [vid] = edge_ids;
    }
}

void erase_from_edge_map (vert2edge_map_t &vert2edge_map, osm_id_t vid,
        unsigned int eid)
{
    std::set <unsigned int> edge_ids = vert2edge_map [vid];
    if (edge_ids.find (eid) != edge_ids.end ())
    {
        edge_ids.erase (eid);
        vert2edge_map [vid] = edge_ids;
    }
}

void graph_from_df (Rcpp::DataFrame gr, vertex_map_t &vm,
        edge_map_t &edge_map, vert2edge_map_t &vert2edge_map,
        bool is_spatial)
{
    Rcpp::StringVector from = gr ["from_id"];
    Rcpp::StringVector to = gr ["to_id"];
    Rcpp::StringVector hw;
    Rcpp::NumericVector from_lon, from_lat, to_lon, to_lat;
    if (is_spatial)
    {
        from_lon = gr ["from_lon"];
        from_lat = gr ["from_lat"];
        to_lon = gr ["to_lon"];
        to_lat = gr ["to_lat"];
        hw = gr ["highway"];
    } else
    {
        hw = Rcpp::StringVector (from.size (), "");
    }
    Rcpp::NumericVector edge_id = gr ["edge_id"];
    Rcpp::NumericVector dist = gr ["d"];
    Rcpp::NumericVector weight = gr ["d_weighted"];

    for (int i = 0; i < to.length (); i ++)
    {
        edge_id (i) -= 1; // now 0-indexed!

        osm_id_t from_id = std::string (from [i]);
        osm_id_t to_id = std::string (to [i]);

        if (vm.find (from_id) == vm.end ())
        {
            osm_vertex_t fromV = osm_vertex_t ();
            if (is_spatial)
            {
                fromV.set_lat (from_lat [i]);
                fromV.set_lon (from_lon [i]);
            }
            vm.emplace (from_id, fromV);
        }
        osm_vertex_t from_vtx = vm.at (from_id);
        from_vtx.add_neighbour_out (to_id);
        vm [from_id] = from_vtx;

        if (vm.find (to_id) == vm.end ())
        {
            osm_vertex_t toV = osm_vertex_t ();
            if (is_spatial)
            {
                toV.set_lat (to_lat [i]);
                toV.set_lon (to_lon [i]);
            }
            vm.emplace (to_id, toV);
        }
        osm_vertex_t to_vtx = vm.at (to_id);
        to_vtx.add_neighbour_in (from_id);
        vm [to_id] = to_vtx;

        std::set <unsigned int> replacement_edges;
        osm_edge_t edge = osm_edge_t (from_id, to_id, dist [i], weight [i],
                std::string (hw [i]), edge_id [i], replacement_edges);
        edge_map.emplace (edge_id [i], edge);
        add_to_edge_map (vert2edge_map, from_id, edge_id [i]);
        add_to_edge_map (vert2edge_map, to_id, edge_id [i]);
    }
}

int get_largest_graph_component (vertex_map_t &v,
        std::unordered_map <osm_id_t, int> &com)
{
    int largest_id = -1;

    // initialize components map
    for (auto it = v.begin (); it != v.end (); ++ it)
        com.insert (std::make_pair (it -> first, -1));

    std::unordered_set <osm_id_t> all_verts, component, nbs_todo, nbs_done;
    for (auto it = v.begin (); it != v.end (); ++ it)
        all_verts.insert (it -> first);
    osm_id_t vt = (*all_verts.begin ());
    nbs_todo.insert (vt);
    int compnum = 0;
    while (all_verts.size () > 0)
    {
        vt = (*nbs_todo.begin ());
        component.insert (vt);
        com.at (vt) = compnum;
        all_verts.erase (vt);

        osm_vertex_t vtx = v.find (vt)->second;
        std::unordered_set <osm_id_t> nbs = vtx.get_all_neighbours ();
        for (auto n: nbs)
        {
            component.insert (n);
            com.at (vt) = compnum;
            if (nbs_done.find (n) == nbs_done.end ())
                nbs_todo.insert (n);
        }
        nbs_todo.erase (vt);
        nbs_done.insert (vt);

        if (nbs_todo.size () == 0 && all_verts.size () > 0)
        {
            nbs_todo.insert (*all_verts.begin ());
            compnum++;
        }
    }

    std::vector <int> comp_sizes (compnum, 0);
    for (auto c: com)
        comp_sizes [c.second]++;
    auto maxi = std::max_element (comp_sizes.begin (), comp_sizes.end ());

    largest_id = std::distance (comp_sizes.begin (), maxi);

    return largest_id;
}


// See docs/graph-contraction for explanation of the following code and
// associated vertex and edge maps.
void contract_graph (vertex_map_t &vertex_map, edge_map_t &edge_map,
        vert2edge_map_t &vert2edge_map)
{
    std::unordered_set <osm_id_t> verts;
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
        std::unordered_set <osm_id_t>::iterator vid = verts.begin ();
        osm_id_t vtx_id = vertex_map.find (*vid)->first;
        osm_vertex_t vtx = vertex_map.find (*vid)->second;
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
            auto nbs = vtx.get_all_neighbours (); // unordered_set <osm_id_t>
            std::vector <osm_id_t> two_nbs;
            for (osm_id_t nb: nbs)
                two_nbs.push_back (nb); // size is always 2

            osm_vertex_t vt0 = vertex_map [two_nbs [0]];
            osm_vertex_t vt1 = vertex_map [two_nbs [1]];
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
                // NOTE: There is no check that types of highways are consistent!
                std::string hw;
                osm_id_t vt_from, vt_to;

                // insert "in" edge first
                std::set <unsigned int> old_edges;

                for (unsigned int e: edges)
                {
                    if (!edges_done [e])
                    {
                        osm_edge_t ei = edge_map.find (e)->second;
                        vt_from = ei.get_from_vertex ();
                        if (vt_from == two_nbs [0] || vt_from == two_nbs [1])
                        {
                            d += ei.dist;
                            w += ei.weight;
                            hw = ei.highway;
                            if (ei.get_old_edges().size() > 0)
                                old_edges = ei.get_old_edges();
                            else
                                old_edges.insert (e);

                            erase_from_edge_map (vert2edge_map, vtx_id, e);
                            add_to_edge_map (vert2edge_map, vtx_id, ne);
                            erase_from_edge_map (vert2edge_map, vt_from, e);
                            add_to_edge_map (vert2edge_map, vt_from, ne);

                            osm_vertex_t vt = vertex_map [vt_from];
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
                        osm_edge_t ei = edge_map.find (e)->second;
                        vt_to = ei.get_to_vertex ();
                        if (vt_to == two_nbs [0] || vt_to == two_nbs [1])
                        {
                            d += ei.dist;
                            w += ei.weight;
                            hw = ei.highway;
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

                            osm_vertex_t vt = vertex_map [vt_to];
                            vt.replace_neighbour (vtx_id, vt_from);
                            vertex_map [vt_to] = vt;

                            edges_done [e] = true;
                            break;
                        }
                    }
                }
                vert2edge_map[vtx_id].insert(ne);
                osm_edge_t new_edge = osm_edge_t (vt_from, vt_to,
                        d, w, hw, ne, old_edges);
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


//' rcpp_sample_graph
//'
//' Randomly sample one connected componnent of a graph
//'
//' @param graph graph to be processed
//' @param nverts_to_sample Number of vertices to sample
//' @param e0 Random edge of graph from which to get first vertex to include in
//' sample
//' @param is_spatial Is the graph spatial or not?
//' @param quiet If TRUE, display progress
//'
//' @return Smaller sub-set of \code{graph}
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_sample_graph (Rcpp::DataFrame graph,
        unsigned int nverts_to_sample, unsigned int e0, bool is_spatial)
{
    std::random_device rd;
    std::mt19937 rng (rd()); // mersenne twister

    vertex_map_t vertices;
    edge_map_t edge_map;
    std::unordered_map <osm_id_t, int> components;
    vert2edge_map_t vert2edge_map;

    graph_from_df (graph, vertices, edge_map, vert2edge_map, is_spatial);
    int largest_component = get_largest_graph_component (vertices, components);

    if (edge_map.find (e0) == edge_map.end ())
        throw std::runtime_error ("edge number not in range of graph");

    Rcpp::NumericVector index;
    if (vertices.size () <= nverts_to_sample)
        return index;

    osm_id_t this_vert;
    bool in_largest = false;
    while (!in_largest)
    {
        osm_edge_t this_edge = edge_map.find (e0++)->second;
        this_vert = this_edge.get_from_vertex ();
        if (components [this_vert] == largest_component)
            in_largest = true;
        if (e0 >= edge_map.size ())
            e0 = 0;
    }

    // Samples are built by randomly tranwing a vertex list, and inspecting
    // edges that extend from it. The only effective way to randomly sample a
    // C++ container is with a std::vector, even though that requires a
    // std::find each time prior to insertion. It's also useful to know the
    // size, so this vector is **NOT** reserved, even though it easily could be.
    // Maybe not the best solution?
    std::vector <osm_id_t> vertlist;
    std::unordered_set <unsigned int> edgelist;
    vertlist.push_back (this_vert);

    while (vertlist.size () < nverts_to_sample)
    {
        // initialise random int generator:
        // TODO: Is this quicker to use a single unif and round each time?
        std::uniform_int_distribution <int> uni (0, vertlist.size () - 1);
        unsigned int randv = uni (rng);
        this_vert = vertlist [randv];

        std::set <unsigned int> edges = vert2edge_map [this_vert];
        for (auto e: edges)
        {
            edgelist.insert (e);
            osm_edge_t this_edge = edge_map.find (e)->second;
            osm_id_t vt = this_edge.get_from_vertex ();
            if (std::find (vertlist.begin(), vertlist.end(), vt) ==
                    vertlist.end())
                vertlist.push_back (vt);
            vt = this_edge.get_to_vertex ();
            if (std::find (vertlist.begin(), vertlist.end(), vt) ==
                    vertlist.end())
                vertlist.push_back (vt);
        }
    }

    int nedges = edgelist.size ();

    // edgelist is an unordered set, so has to be iteratively inserted
    index = Rcpp::NumericVector (nedges);
    unsigned int i = 0;
    for (auto e: edgelist)
        index (i++) = e;

    return index;
}

//' rcpp_get_components
//'
//' Get component numbers for each edge of graph
//'
//' @param graph graph to be processed
//' @param is_spatial Is the graph spatial or not?
//'
//' @return Vector of component numbers, one for each edge.
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_get_components (Rcpp::DataFrame graph)
{
    vertex_map_t vertices;
    edge_map_t edge_map;
    std::unordered_map <osm_id_t, int> components;
    vert2edge_map_t vert2edge_map;

    bool is_spatial = true;

    graph_from_df (graph, vertices, edge_map, vert2edge_map, is_spatial);
    int largest_component = get_largest_graph_component (vertices, components);

    // Then map component numbers of vertices onto edges
    Rcpp::NumericVector ret (edge_map.size (), -1);
    for (auto v: vertices)
    {
        osm_id_t vi = v.first;
        std::set <unsigned int> edges = vert2edge_map [vi];
        for (unsigned int e: edges)
        {
            //if (e >= edge_map.size ())
            //    throw std::runtime_error ("edge number exceeds graph size");
            //else
            if (e < edge_map.size ())
                ret (e) = components [vi] + 1; // 1-indexed!
        }
    }

    return ret;
}

//' rcpp_make_compact_graph
//'
//' Removes nodes and edges from a graph that are not needed for routing
//'
//' @param graph graph to be processed
//' @param is_spatial Is the graph spatial or not?
//' @param quiet If TRUE, display progress
//'
//' @return \code{Rcpp::List} containing one \code{data.frame} with the compact
//' graph, one \code{data.frame} with the original graph and one
//' \code{data.frame} containing information about the relating edge ids of the
//' original and compact graph.
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::List rcpp_make_compact_graph (Rcpp::DataFrame graph, 
        bool is_spatial, bool quiet)
{
    vertex_map_t vertices;
    edge_map_t edge_map;
    std::unordered_map <osm_id_t, int> components;
    vert2edge_map_t vert2edge_map;

    if (!quiet)
    {
        Rcpp::Rcout << "Constructing graph ... ";
        Rcpp::Rcout.flush ();
    }
    graph_from_df (graph, vertices, edge_map, vert2edge_map, is_spatial);
    if (!quiet)
    {
        Rcpp::Rcout << std::endl << "Determining connected components ... ";
        Rcpp::Rcout.flush ();
    }
    int largest_component = get_largest_graph_component (vertices, components);

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
    Rcpp::NumericVector from_lat_vec (nedges), from_lon_vec (nedges),
        to_lat_vec (nedges), to_lon_vec (nedges), dist_vec (nedges),
        weight_vec (nedges), edgeid_vec (nedges);

    int map_size = 0; // size of edge map contracted -> original
    int edge_count = 0;
    for (auto e = edge_map2.begin (); e != edge_map2.end (); ++e)
    {
        osm_id_t from = e->second.get_from_vertex ();
        osm_id_t to = e->second.get_to_vertex ();
        osm_vertex_t from_vtx = vertices2.at (from);
        osm_vertex_t to_vtx = vertices2.at (to);

        from_vec (edge_count) = from;
        to_vec (edge_count) = to;
        highway_vec (edge_count) = e->second.highway;
        dist_vec (edge_count) = e->second.dist;
        weight_vec (edge_count) = e->second.weight;
        from_lat_vec (edge_count) = from_vtx.getLat ();
        from_lon_vec (edge_count) = from_vtx.getLon ();
        to_lat_vec (edge_count) = to_vtx.getLat ();
        to_lon_vec (edge_count) = to_vtx.getLon ();
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
            Rcpp::Named ("d_weighted") = weight_vec,
            Rcpp::Named ("from_lat") = from_lat_vec,
            Rcpp::Named ("from_lon") = from_lon_vec,
            Rcpp::Named ("to_lat") = to_lat_vec,
            Rcpp::Named ("to_lon") = to_lon_vec,
            Rcpp::Named ("highway") = highway_vec);

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
        Rcpp::DataFrame compactgraph, std::vector <int> pts_to_insert,
        bool is_spatial)
{
    vertex_map_t vertices;
    edge_map_t edge_map;
    std::unordered_map <osm_id_t, int> components;
    int largest_component;
    vert2edge_map_t vert2edge_map;

    graph_from_df (fullgraph, vertices, edge_map, vert2edge_map, is_spatial);

    return Rcpp::List::create ();
}
