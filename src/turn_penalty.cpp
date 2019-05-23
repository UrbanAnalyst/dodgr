#include "turn_penalty.h"

/* This code adds time penalties for turning across oncoming traffic in
 * insersections. Works by sorting all outgoing edges in a clockwise direction,
 * and applying penalties to edges with sorted indices >= 2. (Right-side travel
 * simply reverses indices to effectively be anti-clockwise. Incoming edges are
 * also sorted, but sort order is not used.) Implements the following steps:
 *
 * 1. The `fill_edges` function creates unordered_map between centre vertices of
 *    junctions and a `std::pair <RTEdgeSet, RTEdgeSet>` of (incoming, outgoing)
 *    vertices, where `RTEdgeSet = std::set <OneEdge, clockwise_sort>`. This
 *    function also fills an unordered_set of `junction_vertices`.
 * 2. The `replace_junctions` function scans through that map and creates new
 *    edges connecting the incoming vertex to the outgoing vertex in each
 *    direction, and applies penalties to all edges with sorted order > 2.
 */


void routetimes::fill_edges (const Rcpp::DataFrame &graph,
        std::unordered_map <std::string,
                            std::pair <RTEdgeSet, RTEdgeSet> > &the_edges,
        std::unordered_set <std::string> &junction_vertices)
{
    std::vector <std::string> vx0 = graph [".vx0"],
                              vx1 = graph [".vx1"],
                              edge_ = graph ["edge_"];
    std::vector <double> vx0_x = graph [".vx0_x"],
                         vx0_y = graph [".vx0_y"],
                         vx1_x = graph [".vx1_x"],
                         vx1_y = graph [".vx1_y"];

    const size_t n = static_cast <size_t> (graph.nrow ());

    for (size_t i = 0; i < n; i++)
    {
        OneEdge edge;
        edge.v0 = vx0 [i];
        edge.v1 = vx1 [i];
        edge.edge = edge_ [i];
        edge.x = vx1_x [i] - vx0_x [i];
        edge.y = vx1_y [i] - vx0_y [i];

        // the incoming edge: vx0->vx1 with key = vx1
        routetimes::replace_one_map_edge (the_edges, vx1 [i], edge, true);
        // the outgoing edge: vx0->vx1 with key = vx0
        routetimes::replace_one_map_edge (the_edges, vx0 [i], edge, false);
    }

    erase_non_junctions (the_edges, junction_vertices);
}

void routetimes::replace_one_map_edge (
        std::unordered_map <std::string,
                            std::pair <RTEdgeSet, RTEdgeSet> > &the_edges,
        std::string key, OneEdge edge, bool incoming)
{
    std::pair <RTEdgeSet, RTEdgeSet> edge_set;
    if (the_edges.find (key) != the_edges.end ())
    {
        edge_set = the_edges.at (key);
        the_edges.erase (key);
    }
    if (incoming)
        edge_set.first.emplace (edge);
    else
        edge_set.second.emplace (edge);
    
    the_edges.emplace (key, edge_set);
}

/* Remove all edge items with < 3 outgoing edges. A junction has 2 or more
 * out_edges / in_edges, but presume that flow in and out of 1->2 junctions is
 * regulated by lights, and that only proper cross-junctions (1->3) incur
 * additional waiting times. The result: min_edges = 4
 */
void routetimes::erase_non_junctions (
        std::unordered_map <std::string,
                            std::pair <RTEdgeSet, RTEdgeSet> > &the_edges,
        std::unordered_set <std::string> &junction_vertices)
{
    const int min_edges = 4;

    std::unordered_set <std::string> removes;
    for (auto e: the_edges)
    {
        if (e.second.second.size () < min_edges) // set of outgoing edges
            removes.emplace (e.first);
        else
        {
            for (auto r: e.second.second) // RTEdgeSet
                junction_vertices.emplace (r.v0);
        }
    }
    for (auto r: removes)
        the_edges.erase (r);
}

/* Replace junction mid-points with direct connections from in -> 0 -> out with
 * (in, out) pairs and a binary flag for turning penalty.
 */
void routetimes::replace_junctions (
        const std::unordered_map <std::string,
                                  std::pair <RTEdgeSet, RTEdgeSet> > &the_edges,
        std::vector <OneCompoundEdge> &junctions,
        bool left_side)
{
    // Estimate size of resultant vector - this is not simply n * (n - 1),
    // because one way streeets mean in-edges aren't always the same as out, so
    // numbers have to be explicitly counted.
    size_t n_junctions = 0;
    for (auto e: the_edges)
        for (auto ei: e.second.first)  // incoming edge
        {
            std::unordered_set <std::string> out_edges;
            for (auto ej: e.second.second) // outgoing edge
                if (ej.edge != ei.edge)
                    out_edges.emplace (ej.edge);

            for (auto ej: e.second.second)
                if (out_edges.find (ej.edge) != out_edges.end ())
                    n_junctions++;
        }
    junctions.resize (n_junctions);

    size_t count = 0;
    for (auto e: the_edges)
    {
        // iterate over each incoming edge, and establish penalties for all
        // other ongoing.
        for (auto ei: e.second.first)  // incoming edge
        {
            std::unordered_map <std::string, size_t> out_edges;
            size_t i = 0;
            for (auto ej: e.second.second) // outgoing edge
                if (ej.edge != ei.edge)
                {
                    out_edges.emplace (ej.edge, i++);
                }
            // For right-side travel, simply reverse the clockwise sequence of
            // values, to give anti-clockwise sequential counts:
            if (!left_side)
                for (auto oe: out_edges)
                {
                    out_edges [oe.first] = out_edges.size () - oe.second;
                }

            // Then use those out_edges to construct new edge lists:
            for (auto ej: e.second.second) // outgoing edge
                if (out_edges.find (ej.edge) != out_edges.end ())
                {
                    OneCompoundEdge edge;
                    edge.v0 = ei.v0;
                    edge.v1 = ej.v1;
                    edge.edge0 = ei.edge;
                    edge.edge1 = ej.edge;
                    edge.penalty = false;
                    if (out_edges.at (ej.edge) > 2)
                        edge.penalty = true;
                    junctions [count++] = edge;
                }
        }
    }
}

// This is hard-coded for weight_streetnet.sc structure
Rcpp::DataFrame routetimes::expand_edges (const Rcpp::DataFrame &graph, 
        std::vector <OneCompoundEdge> &junctions, int turn_penalty)
{
    const size_t hash_len = 10; // for new edge IDs

    std::unordered_map <std::string, double> map_x, map_y, map_d, map_dw,
        map_time, map_timew;
    std::unordered_map <std::string, std::string> map_object, map_highway;

    std::vector <std::string> vx0 = graph [".vx0"],
                              vx1 = graph [".vx1"],
                              edge_ = graph ["edge_"],
                              object = graph ["object_"],
                              highway = graph ["highway"];
    std::vector <double> vx0_x = graph [".vx0_x"],
                         vx0_y = graph [".vx0_y"],
                         vx1_x = graph [".vx1_x"],
                         vx1_y = graph [".vx1_y"],
                         d = graph ["d"],
                         dw = graph ["d_weighted"],
                         time = graph ["time"],
                         timew = graph ["time_weighted"];

    std::unordered_set <std::string> edge_set;
    for (auto j: junctions)
    {
        edge_set.emplace (j.edge0);
        edge_set.emplace (j.edge1);
    }

    for (size_t i = 0; i < static_cast <size_t> (graph.nrow ()); i++)
    {
        if (edge_set.find (edge_ [i]) != edge_set.end ())
        {
            map_x.emplace (vx0 [i], vx0_x [i]);
            map_y.emplace (vx0 [i], vx0_y [i]);
            map_x.emplace (vx1 [i], vx1_x [i]);
            map_y.emplace (vx1 [i], vx1_y [i]);

            map_d.emplace (edge_ [i], d [i]);
            map_dw.emplace (edge_ [i], dw [i]);
            map_time.emplace (edge_ [i], time [i]);
            map_timew.emplace (edge_ [i], timew [i]);

            map_object.emplace (edge_ [i], object [i]);
            map_highway.emplace (edge_ [i], highway [i]);
        }
    }

    const size_t n = junctions.size ();
    std::vector <std::string> vx0_out (n),
                              vx1_out (n),
                              edge_out (n),
                              object_out (n),
                              highway_out (n),
                              old_edge_in (n),
                              old_edge_out (n);
    std::vector <double> vx0_x_out (n),
                         vx0_y_out (n),
                         vx1_x_out (n),
                         vx1_y_out (n),
                         d_out (n),
                         dw_out (n),
                         time_out (n),
                         timew_out (n);

    for (size_t i = 0; i < n; i++)
    {
        vx0_out [i] = junctions [i].v0;
        vx0_x_out [i] = map_x.at (junctions [i].v0);
        vx0_y_out [i] = map_y.at (junctions [i].v0);
        vx1_out [i] = junctions [i].v1;
        vx1_x_out [i] = map_x.at (junctions [i].v1);
        vx1_y_out [i] = map_y.at (junctions [i].v1);
        edge_out [i] = "j_" + sc::random_id (hash_len);
        old_edge_in [i] = junctions [i].edge0;
        old_edge_out [i] = junctions [i].edge1;

        // Map all others to properties of out edge:
        object_out [i] = map_object.at (junctions [i].edge1);
        highway_out [i] = map_highway.at (junctions [i].edge1);

        d_out [i] = map_d.at (junctions [i].edge0) + 
                    map_d.at (junctions [i].edge1);
        dw_out [i] = map_dw.at (junctions [i].edge0) + 
                     map_dw.at (junctions [i].edge1);
        time_out [i] = map_time.at (junctions [i].edge0) +
                       map_time.at (junctions [i].edge1);
        timew_out [i] = map_timew.at (junctions [i].edge0) +
                       map_timew.at (junctions [i].edge1);
        if (junctions [i].penalty)
        {
            time_out [i] += turn_penalty;
            timew_out [i] += turn_penalty;
        }
    }

    Rcpp::DataFrame res = Rcpp::DataFrame::create (
            Rcpp::Named (".vx0") = vx0_out,
            Rcpp::Named (".vx1") = vx1_out,
            Rcpp::Named ("edge_") = edge_out, // empty here
            Rcpp::Named (".vx0_x") = vx0_x_out,
            Rcpp::Named (".vx0_y") = vx0_y_out,
            Rcpp::Named (".vx1_x") = vx1_x_out,
            Rcpp::Named (".vx1_y") = vx1_y_out,
            Rcpp::Named ("d") = d_out,
            Rcpp::Named ("object_") = object_out,
            Rcpp::Named ("highway") = highway_out,
            Rcpp::Named ("d_weighted") = dw_out,
            Rcpp::Named ("time") = time_out,
            Rcpp::Named ("time_weighted") = timew_out,
            Rcpp::Named ("old_edge_in") = old_edge_in,
            Rcpp::Named ("old_edge_out") = old_edge_out,
            Rcpp::_["stringsAsFactors"] = false);

    return res;
}

//' rcpp_route_times
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::List rcpp_route_times (const Rcpp::DataFrame graph,
        bool left_side, int turn_penalty)
{
    std::unordered_map <std::string, std::pair <RTEdgeSet, RTEdgeSet> > the_edges;
    std::unordered_set <std::string> junction_vertices;

    routetimes::fill_edges (graph, the_edges, junction_vertices);
    std::vector <OneCompoundEdge> junctions;
    routetimes::replace_junctions (the_edges, junctions, left_side);

    Rcpp::DataFrame expanded_graph = routetimes::expand_edges (graph,
            junctions, turn_penalty);

    Rcpp::CharacterVector junction_vec (junction_vertices.size ());
    size_t i = 0;
    for (auto j: junction_vertices)
        junction_vec (i++) = j;

    return Rcpp::List::create (
            Rcpp::Named ("graph") = expanded_graph,
            Rcpp::Named ("junction_vertices") = junction_vec);
}
