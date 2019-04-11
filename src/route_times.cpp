#include "route_times.h"

/* This code add time penalties for turning across oncoming traffic in
 * insersections. Works by sorting all outgoing edges in a clockwise direction,
 * and applying penalties to edges with sorted indices >= 2. (Right-side travel
 * simply reverses indices to effectively be anti-clockwise. Incoming edges are
 * also sorted, but sort order is not used.) Implements the following steps:
 *
 * 1. The `fill_edges` function creates unordered_map between centre vertices of
 *    junctions and a `std::pair <RTEdgeSet, RTEdgeSet>` of (incoming, outgoing)
 *    vertices, where `RTEdgeSet = std::set <OneEdge, clockwise_sort>`. This
 *    function also fills an unordered_set of `edges_to_remove`, which are all
 *    original outgoing edges from junction vertices.
 * 2. The `replace_junctions` function scans through that map and creates new
 *    edges connecting the incoming vertex to the outgoing vertex in each
 *    direction, and applies penalties to all edges with sorted order > 2.
 */


void routetimes::fill_edges (const Rcpp::DataFrame &graph,
        std::unordered_map <std::string,
                            std::pair <RTEdgeSet, RTEdgeSet> > &the_edges,
        std::unordered_set <std::string> &edges_to_remove,
        bool ignore_oneway)
{
    std::vector <std::string> vx0 = graph [".vx0"],
                              vx1 = graph [".vx1"],
                              edge_ = graph ["edge_"];
    std::vector <double> vx0_x = graph [".vx0_x"],
                         vx0_y = graph [".vx0_y"],
                         vx1_x = graph [".vx1_x"],
                         vx1_y = graph [".vx1_y"];
    std::vector <bool> oneway = graph ["oneway"];

    const int n = graph.nrow ();

    for (int i = 0; i < n; i++)
    {
        OneEdge edge;
        edge.v0 = vx0 [i];
        edge.v1 = vx1 [i];
        edge.edge = edge_ [i];
        edge.x = vx1_x [i] - vx0_x [i];
        edge.y = vx1_y [i] - vx0_y [i];

        // the incoming edge:
        routetimes::replace_one_map_edge (the_edges, vx1 [i], edge, true);
        // the outgoing edge:
        routetimes::replace_one_map_edge (the_edges, vx0 [i], edge, false);
        // TODO: oneway_bicycle tag, although this is hardly used
        if (ignore_oneway || !oneway [i])
        {
            routetimes::replace_one_map_edge (the_edges, vx1 [i], edge, false);
            routetimes::replace_one_map_edge (the_edges, vx0 [i], edge, true);
        }
    }

    erase_non_junctions (the_edges, edges_to_remove);
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
        std::unordered_set <std::string> &edges_to_remove)
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
                edges_to_remove.emplace (r.edge);
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
    // Estimate size of resultant vector - this will usually be an over-estimate
    // because any one-way streets will yield less than this number of
    // combinations in in/out edges.
    size_t count = 0;
    for (auto e: the_edges)
        count += e.second.second.size () * (e.second.second.size () - 1);
    junctions.resize (count);
    count = 0;

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
    junctions.resize (count - 1);
}

void routetimes::new_graph (graph, junctions)
{
}

//' rcpp_route_times
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::DataFrame rcpp_route_times (const Rcpp::DataFrame graph,
        bool ignore_oneway, bool left_side)
{
    std::unordered_map <std::string, std::pair <RTEdgeSet, RTEdgeSet> > the_edges;
    std::unordered_set <std::string> edges_to_remove;

    routetimes::fill_edges (graph, the_edges, edges_to_remove, ignore_oneway);
    std::vector <OneCompoundEdge> junctions;
    routetimes::replace_junctions (the_edges, junctions, left_side);

    routetimes::new_graph (graph, junctions);

    // dummy values to demonstrate that clockwise sorting works.
    RTEdgeSet x;
    OneEdge e;

    e.x = 1.0;  e.y = -3.0;     e.v0 = "0"; e.v1 = "1";    x.emplace (e);
    e.x = 0.0;  e.y = 1.0;      e.v0 = "0"; e.v1 = "2";    x.emplace (e);
    e.x = -1.0; e.y = 2.0;      e.v0 = "0"; e.v1 = "3";    x.emplace (e);
    e.x = 2.0;  e.y = -0.5;     e.v0 = "0"; e.v1 = "4";    x.emplace (e);
    e.x = 3.0;  e.y = 1.5;      e.v0 = "0"; e.v1 = "5";    x.emplace (e);
    e.x = -2.5; e.y = -0.7;     e.v0 = "0"; e.v1 = "6";    x.emplace (e);

    const size_t len = x.size ();
    Rcpp::NumericVector xout (len), yout (len);
    Rcpp::CharacterVector v0 (len), v1 (len);
    int count = 0;
    for (auto i: x)
    {
        xout [count] = i.x;
        yout [count] = i.y;
        v0 [count] = i.v0;
        v1 [count++] = i.v1;
    }
    Rcpp::DataFrame res = Rcpp::DataFrame::create (
            Rcpp::Named ("x") = xout,
            Rcpp::Named ("y") = yout,
            Rcpp::Named ("v0") = v0,
            Rcpp::Named ("v1") = v1,
            Rcpp::_["stringsAsFactors"] = false);

    return res;
}
