#include "graph.h"


//' sample_one_edge_no_comps
//'
//' Sample one edge for graph that has no pre-calculated components. Only used
//' in \code{sample_one_vertex}
//'
//' @param edge_map edge_map
//' @return std::vector of 2 elements: [0] with value of largest connected 
//' component; [1] with random index to one edge that is part of that component.
//' @noRd
edge_component sample_one_edge_no_comps (vertex_map_t &vertices,
        edge_map_t &edge_map)
{
    // TDOD: FIX edge_id_t type defs here!
    std::unordered_map <vertex_id_t, unsigned int> components;
    std::random_device rd;
    std::mt19937 rng (rd()); // mersenne twister

    unsigned int largest_component = identify_graph_components (vertices, components);

    bool in_largest = false;
    std::uniform_int_distribution <unsigned int> uni0 (0, edge_map.size () - 1);
    unsigned int e0 = uni0 (rng);
    while (!in_largest)
    {
        // TODO: The following is an O(N) lookup; maybe just use
        // edge_map.find (std::to_string (e0)) and handle not found; that'd be
        // constant time.
        edge_t this_edge = std::next (edge_map.begin (), e0++)->second;
        vertex_id_t this_vert = this_edge.get_from_vertex ();
        if (components [this_vert] == largest_component)
            in_largest = true;
        if (e0 >= edge_map.size ())
            e0 = 0;
    }

    edge_component edge_comp;
    edge_comp.component = largest_component;
    edge_comp.edge = std::next (edge_map.begin (), e0)->first;

    return edge_comp;
}

//' sample_one_edge_with_comps
//'
//' Sample one edge for graph that has pre-calculated components. Only used in
//' \code{sample_one_vertex}
//'
//' @param edge_map edge_map
//' @return Random index to one edge that is part of the largest connected
//' component.
//' @noRd
edge_id_t sample_one_edge_with_comps (Rcpp::DataFrame graph,
        edge_map_t &edge_map)
{
    std::random_device rd;
    std::mt19937 rng (rd()); // mersenne twister

    Rcpp::NumericVector component = graph ["component"];
    std::uniform_int_distribution <unsigned int> uni (0, graph.nrow () - 1);
    unsigned int e0 = uni (rng);
    while (component (e0) > 1)
        e0 = uni (rng);

    return std::next (edge_map.begin (), e0)->first;
}

vertex_id_t select_random_vert (Rcpp::DataFrame graph,
        edge_map_t &edge_map, vertex_map_t &vertices)
{
    vertex_id_t this_vert;
    if (graph_has_components (graph))
    {
        edge_id_t e0 = sample_one_edge_with_comps (graph, edge_map);
        edge_t this_edge = edge_map.find (e0)->second;
        this_vert = this_edge.get_from_vertex ();
    } else
    {
        edge_component edge_comp = sample_one_edge_no_comps (vertices, edge_map);
        edge_t this_edge = edge_map.find (edge_comp.edge)->second; // random edge
        this_vert = this_edge.get_from_vertex ();
    }

    return this_vert;
}


//' rcpp_sample_graph
//'
//' Randomly sample one connected componnent of a graph
//'
//' @param graph graph to be processed
//' @param nverts_to_sample Number of vertices to sample
//' @param quiet If TRUE, display progress
//'
//' @return Smaller sub-set of \code{graph}
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::StringVector rcpp_sample_graph (Rcpp::DataFrame graph,
        unsigned int nverts_to_sample)
{
    std::random_device rd;
    std::default_random_engine rng (rd()); // safest to use here

    vertex_map_t vertices;
    edge_map_t edge_map;
    vert2edge_map_t vert2edge_map;

    graph_from_df (graph, vertices, edge_map, vert2edge_map);

    Rcpp::StringVector edges_out;
    if (vertices.size () <= nverts_to_sample)
        return edges_out; // return empty vector

    std::unordered_map <vertex_id_t, unsigned int> components;
    unsigned int largest_component = identify_graph_components (vertices, components);
    // simple vert_ids set for quicker random selection:
    std::unordered_set <vertex_id_t> vert_ids;
    for (auto v: vertices)
        if (components [v.first] == largest_component)
            vert_ids.emplace (v.first);
    if (vert_ids.size () < nverts_to_sample)
    {
        Rcpp::Rcout << "Largest connected component only has " <<
            vert_ids.size () << " vertices" << std::endl;
        nverts_to_sample = vert_ids.size ();
    }

    vertex_id_t this_vert = select_random_vert (graph, edge_map, vertices);

    // Samples are built by randomly tranwing a vertex list, and inspecting
    // edges that extend from it. The only effective way to randomly sample a
    // C++ container is with a std::vector, even though that requires a
    // std::find each time prior to insertion. It's also useful to know the
    // size, so this vector is **NOT** reserved, even though it easily could be.
    // Maybe not the best solution?
    std::vector <vertex_id_t> vertlist;
    std::unordered_set <edge_id_t> edgelist;
    vertlist.push_back (this_vert);


    unsigned int count = 0;
    while (vertlist.size () < nverts_to_sample)
    {
        // initialise random int generator:
        // TODO: Is this quicker to use a single unif and round each time?
        std::uniform_int_distribution <unsigned int> uni (0, vertlist.size () - 1);
        unsigned int randv = uni (rng);
        this_vert = vertlist [randv];

        std::unordered_set <edge_id_t> edges = vert2edge_map [this_vert];
        for (auto e: edges)
        {
            edgelist.insert (e);
            edge_t this_edge = edge_map.find (e)->second;
            vertex_id_t vt = this_edge.get_from_vertex ();
            if (std::find (vertlist.begin(), vertlist.end(), vt) ==
                    vertlist.end())
                vertlist.push_back (vt);
            if (vertlist.size () >= nverts_to_sample)
            {
                break;
            } else
            {
                vt = this_edge.get_to_vertex ();
                if (std::find (vertlist.begin(), vertlist.end(), vt) ==
                        vertlist.end())
                    vertlist.push_back (vt);
            }
            if (vertlist.size () >= nverts_to_sample)
                break;
        }
        count++;
        if (count > (100 * nverts_to_sample))
        {
            // likely stuck in some one-way part of graph that can't connect, so
            // reset to another start node
            edgelist.clear ();
            this_vert = select_random_vert (graph, edge_map, vertices);
            vertlist.clear ();
            vertlist.push_back (this_vert);
            count = 0;
        }
    }

    unsigned int nedges = edgelist.size ();

    // edgelist is an unordered set, so has to be iteratively inserted
    edges_out = Rcpp::StringVector (nedges);
    unsigned int i = 0;
    for (auto e: edgelist)
        edges_out (i++) = e;

    return edges_out;
}
