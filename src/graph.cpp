#include "graph.h"


void add_to_edge_map (vert2edge_map_t &vert2edge_map, vertex_id_t vid,
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

void erase_from_edge_map (vert2edge_map_t &vert2edge_map, vertex_id_t vid,
        unsigned int eid)
{
    std::set <unsigned int> edge_ids = vert2edge_map [vid];
    if (edge_ids.find (eid) != edge_ids.end ())
    {
        edge_ids.erase (eid);
        vert2edge_map [vid] = edge_ids;
    }
}
//' graph_from_df
//'
//' Convert a standard graph data.frame into an object of class graph. Graphs
//' are standardised with the function \code{convert_graph()$graph}, and contain
//' only the four columns [from, to, d, w]
//' @noRd
void graph_from_df (Rcpp::DataFrame gr, vertex_map_t &vm,
        edge_map_t &edge_map, vert2edge_map_t &vert2edge_map)
{
    if (gr.ncol () != 5)
        throw std::runtime_error ("graph must have 5 columns: run convert_graph() first");

    Rcpp::NumericVector edge_id = gr ["edge_id"];
    Rcpp::StringVector from = gr ["from"];
    Rcpp::StringVector to = gr ["to"];
    Rcpp::NumericVector dist = gr ["d"];
    Rcpp::NumericVector weight = gr ["w"];

    for (int i = 0; i < to.length (); i ++)
    {
        edge_id (i) -= 1; // now 0-indexed!

        vertex_id_t from_id = std::string (from [i]);
        vertex_id_t to_id = std::string (to [i]);

        if (vm.find (from_id) == vm.end ())
        {
            vertex_t fromV = vertex_t ();
            vm.emplace (from_id, fromV);
        }
        vertex_t from_vtx = vm.at (from_id);
        from_vtx.add_neighbour_out (to_id);
        vm [from_id] = from_vtx;

        if (vm.find (to_id) == vm.end ())
        {
            vertex_t toV = vertex_t ();
            vm.emplace (to_id, toV);
        }
        vertex_t to_vtx = vm.at (to_id);
        to_vtx.add_neighbour_in (from_id);
        vm [to_id] = to_vtx;

        std::set <unsigned int> replacement_edges;
        edge_t edge = edge_t (from_id, to_id, dist [i], weight [i],
                edge_id [i], replacement_edges);
        edge_map.emplace (edge_id [i], edge);
        add_to_edge_map (vert2edge_map, from_id, edge_id [i]);
        add_to_edge_map (vert2edge_map, to_id, edge_id [i]);
    }
}

int identify_graph_components (vertex_map_t &v,
        std::unordered_map <vertex_id_t, int> &com)
{
    int largest_id = -1;

    // initialize components map
    for (auto it = v.begin (); it != v.end (); ++ it)
        com.insert (std::make_pair (it -> first, -1));

    std::unordered_set <vertex_id_t> all_verts, component, nbs_todo, nbs_done;
    for (auto it = v.begin (); it != v.end (); ++ it)
        all_verts.insert (it -> first);
    vertex_id_t vt = (*all_verts.begin ());
    nbs_todo.insert (vt);
    int compnum = 0;
    while (all_verts.size () > 0)
    {
        vt = (*nbs_todo.begin ());
        component.insert (vt);
        com.at (vt) = compnum;
        all_verts.erase (vt);

        vertex_t vtx = v.find (vt)->second;
        std::unordered_set <vertex_id_t> nbs = vtx.get_all_neighbours ();
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


//' rcpp_get_component_vector
//'
//' Get component numbers for each edge of graph
//'
//' @param graph graph to be processed
//'
//' @return Vector of component numbers, one for each edge.
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_get_component_vector (Rcpp::DataFrame graph)
{
    vertex_map_t vertices;
    edge_map_t edge_map;
    std::unordered_map <vertex_id_t, int> components;
    vert2edge_map_t vert2edge_map;

    graph_from_df (graph, vertices, edge_map, vert2edge_map);
    int largest_component = identify_graph_components (vertices, components);
    largest_component++; // suppress unused variable warning

    // Then map component numbers of vertices onto edges
    Rcpp::NumericVector ret (edge_map.size (), -1);
    for (auto v: vertices)
    {
        vertex_id_t vi = v.first;
        std::set <unsigned int> edges = vert2edge_map [vi];
        for (unsigned int e: edges)
        {
            if (e >= edge_map.size ())
                throw std::runtime_error ("edge number exceeds graph size");
            else
                ret (e) = components [vi] + 1; // 1-indexed!
        }
    }

    return ret;
}
