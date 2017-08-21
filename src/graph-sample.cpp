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
std::vector <unsigned int>  sample_one_edge_no_comps (vertex_map_t &vertices,
        edge_map_t &edge_map)
{
    // TDOD: FIX edge_id_t type defs here!
    std::unordered_map <vertex_id_t, int> components;
    std::random_device rd;
    std::mt19937 rng (rd()); // mersenne twister

    int largest_component = identify_graph_components (vertices, components);

    bool in_largest = false;
    std::uniform_int_distribution <unsigned int> uni0 (0, edge_map.size () - 1);
    unsigned int e0 = uni0 (rng);
    while (!in_largest)
    {
        edge_t this_edge = edge_map.find (e0++)->second;
        vertex_id_t this_vert = this_edge.get_from_vertex ();
        if (components [this_vert] == largest_component)
            in_largest = true;
        if (e0 >= edge_map.size ())
            e0 = 0;
    }

    std::vector <unsigned int> res;
    res.reserve (2);
    res [0] = (unsigned int) largest_component;
    res [1] = e0;

    return res;
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
edge_id_t sample_one_edge_with_comps (Rcpp::DataFrame graph)
{
    std::random_device rd;
    std::mt19937 rng (rd()); // mersenne twister

    Rcpp::NumericVector component = graph ["component"];
    std::uniform_int_distribution <unsigned int> uni (0, graph.nrow () - 1);
    unsigned int e0 = uni (rng);
    while (component (e0) > 1)
        e0 = uni (rng);

    return e0;
}

//' graph_has_components
//'
//' Does a graph have a vector of connected component IDs? Only used in
//' \code{sample_one_vertex}
//' @noRd
bool graph_has_components (Rcpp::DataFrame graph)
{
    Rcpp::CharacterVector graph_names = graph.attr ("names");
    bool has_comps = false;
    for (auto n: graph_names)
        if (n == "component")
            has_comps = true;

    return has_comps;
}

//' sample_one_vertex
//' 
//' Get a random vertex in graph that is part of the largest connected component
//' @noRd
vertex_id_t sample_one_vertex (Rcpp::DataFrame graph, vertex_map_t &vertices,
        edge_map_t &edge_map)
{
    vertex_id_t this_vert;

    if (graph_has_components (graph))
    {
        unsigned int e0 = sample_one_edge_with_comps (graph);
        edge_t this_edge = edge_map.find (e0)->second;
        this_vert = this_edge.get_from_vertex ();
    } else
    {
        std::vector <unsigned int> es;
        es = sample_one_edge_no_comps (vertices, edge_map);
        edge_t this_edge = edge_map.find (es [1])->second; // random edge
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
Rcpp::NumericVector rcpp_sample_graph (Rcpp::DataFrame graph,
        unsigned int nverts_to_sample)
{
    std::random_device rd;
    std::mt19937 rng (rd()); // mersenne twister

    vertex_map_t vertices;
    edge_map_t edge_map;
    std::unordered_map <vertex_id_t, int> components;
    vert2edge_map_t vert2edge_map;

    graph_from_df (graph, vertices, edge_map, vert2edge_map);

    Rcpp::NumericVector index;
    if (vertices.size () <= nverts_to_sample)
        return index;

    vertex_id_t this_vert;
    if (graph_has_components (graph))
    {
        unsigned int e0 = sample_one_edge_with_comps (graph);
        edge_t this_edge = edge_map.find (e0)->second;
        this_vert = this_edge.get_from_vertex ();
    } else
    {
        std::vector <unsigned int> es = sample_one_edge_no_comps (vertices, edge_map);
        edge_t this_edge = edge_map.find (es [1])->second; // random edge
        this_vert = this_edge.get_from_vertex ();
    }

    // Samples are built by randomly tranwing a vertex list, and inspecting
    // edges that extend from it. The only effective way to randomly sample a
    // C++ container is with a std::vector, even though that requires a
    // std::find each time prior to insertion. It's also useful to know the
    // size, so this vector is **NOT** reserved, even though it easily could be.
    // Maybe not the best solution?
    std::vector <vertex_id_t> vertlist;
    std::unordered_set <edge_id_t> edgelist;
    vertlist.push_back (this_vert);

    while (vertlist.size () < nverts_to_sample)
    {
        // initialise random int generator:
        // TODO: Is this quicker to use a single unif and round each time?
        std::uniform_int_distribution <unsigned int> uni (0, vertlist.size () - 1);
        unsigned int randv = uni (rng);
        this_vert = vertlist [randv];

        std::set <edge_id_t> edges = vert2edge_map [this_vert];
        for (auto e: edges)
        {
            edgelist.insert (e);
            edge_t this_edge = edge_map.find (e)->second;
            vertex_id_t vt = this_edge.get_from_vertex ();
            if (std::find (vertlist.begin(), vertlist.end(), vt) ==
                    vertlist.end())
                vertlist.push_back (vt);
            if (vertlist.size () < nverts_to_sample)
            {
                vt = this_edge.get_to_vertex ();
                if (std::find (vertlist.begin(), vertlist.end(), vt) ==
                        vertlist.end())
                    vertlist.push_back (vt);
            }
        }
    }

    unsigned int nedges = edgelist.size ();

    // edgelist is an unordered set, so has to be iteratively inserted
    index = Rcpp::NumericVector (nedges);
    unsigned int i = 0;
    for (auto e: edgelist)
        index (i++) = e;

    return index;
}
