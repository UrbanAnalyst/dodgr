
#include "route_times.h"


// https://stackoverflow.com/questions/6989100/sort-points-in-clockwise-order

// x and y values pre-converted by subtracting the centre values
bool routetimes::isLess (OneNode a, OneNode b)
{
    if (a.x >= 0.0 && b.x < 0.0)
        return true;
    if (a.x < 0.0 && b.x >= 0)
        return false;
    if (a.x == 0.0 && b.x == 0.0)
    {
        if (a.y >= 0.0 || b.y >= 0.0)
            return a.y > b.y;
        return b.y > a.y;
    }

    double det = a.x * b.y - a.y * b.x;
    if (det < 0)
        return true;
    if (det > 0)
        return false;

    double d1 = a.x * a.x + a.y * a.y;
    double d2 = b.x * b.x + b.y * b.y;
    return d1 > d2;
}

void routetimes::replace_one_map_edge (
        std::unordered_map <std::string, std::vector <std::string> > &the_edges,
        std::string key, std::string value)
{
    std::vector <std::string> edges;
    if (the_edges.find (key) != the_edges.end ())
    {
        edges = the_edges.at (key);
        the_edges.erase (key);
    }
    edges.push_back (value);
    the_edges.emplace (key, edges);
}

void routetimes::erase_non_junctions (
        std::unordered_map <std::string, std::vector <std::string> > &the_edges)
{
    // A junction has 2 or more out_edges / in_edges, but presume that flow in
    // and out of 1->2 junctions is regulated by lights, and that only proper
    // cross-junctions (1->3) incur additional waiting times:
    const int min_edges = 3;

    std::unordered_set <std::string> removes;
    for (auto e: the_edges)
        if (e.second.size () < min_edges)
            removes.emplace (e.first);
    for (auto r: removes)
        the_edges.erase (r);
}

void routetimes::fill_edges (const Rcpp::DataFrame &graph,
        std::unordered_map <std::string, double> &x0,
        std::unordered_map <std::string, double> &y0,
        std::unordered_map <std::string, std::vector <std::string> > &out_edges)
{
    std::vector <std::string> vx0 = graph [".vx0"],
                              vx1 = graph [".vx1"];
    std::vector <double> vx0_x = graph [".vx0_x"],
                         vx0_y = graph [".vx0_y"],
                         vx1_x = graph [".vx1_x"],
                         vx1_y = graph [".vx1_y"];

    const int n = graph.nrow ();

    for (int i = 0; i < n; i++)
    {
        x0.emplace (vx0 [i], vx0_x [i]);
        x0.emplace (vx1 [i], vx1_x [i]);
        y0.emplace (vx0 [i], vx0_y [i]);
        y0.emplace (vx1 [i], vx1_y [i]);

        routetimes::replace_one_map_edge (out_edges, vx0 [i], vx1 [i]);
    }

    erase_non_junctions (out_edges);
}

// clockwise sorting of edges
void routetimes::sort_edges (
        const std::unordered_map <std::string, std::vector <std::string> > &edges_in,
        std::unordered_map <std::string, std::vector <std::string> > &edges_sorted,
        const std::unordered_map <std::string, double> &x0,
        const std::unordered_map <std::string, double> &y0)
{
    for (auto e: edges_in)
    {
        std::vector <std::string> edge_set = e.second;
        std::vector <OneNode> nodes (edge_set.size ());
        size_t count = 0;
        bool chk = true;
        for (auto ei: edge_set)
        {
            OneNode n;
            if (x0.find (ei) == x0.end () || y0.find (ei) == y0.end ())
                chk = false;
            else
            {
                n.x = x0.at (ei);
                n.y = y0.at (ei);
                n.id = ei;
                nodes [count++] = n;
            }
        }
        if (!chk)
            continue;

        std::sort (nodes.begin (), nodes.end (), routetimes::isLess);
        std::vector <std::string> sorted_edges (nodes.size ());
        count = 0;
        for (auto n: nodes)
            sorted_edges [count++] = n.id;
        edges_sorted.emplace (e.first, sorted_edges);
    }
}

/* Replace junction mid-points with direct connections from in -> 0 -> out with
 * (in, out) pairs and a binary flag for turning penalty. Each juction in input
 * edges has key = centre, values = clockwise-sorted out vertices. These are
 * replaced with a bunch of standard pairs subsituting single (in->centre) with
 * all combinations of (in->out) and binary flag for turning penalty.  edges_in
 * are sorted version created in sort_edges
 */
void routetimes::replace_junctions (
        const std::unordered_map <std::string, std::vector <std::string> > &edges,
        std::vector <Junction> junctions,
        bool left_side)
{
    size_t count = 0;
    for (auto e: edges)
        count += e.second.size () * (e.second.size () - 1);
    junctions.resize (count);
    count = 0;

    for (auto e: edges)
    {
        // iterate over each outgoing edge, treating that as incoming and
        // establishing penalties for all other ongoing.
        int ei_pos = 0;
        for (auto ei: e.second) 
        {
            std::unordered_map <std::string, int> edge_map;
            int ej_pos = 0;
            for (auto ej: e.second)
                if (ej != ei)
                {
                    int dir = ej_pos++ - ei_pos;
                    if (dir < 0)
                        dir += e.second.size () - 1;
                    /* dir then counts clockwise from 0 = first left turn up to
                     * e.second.size () for last right turn. Penalties can then
                     * be applied to any (dir > 1). For right-side driving,
                     * these dir values are simply reversed: */
                    if (!left_side)
                        dir = e.second.size () - 1 - dir;
                    Junction j;
                    j.in = ei;
                    j.centre = e.first;
                    j.out = ej;
                    j.penalty = false;
                    if (dir > 1)
                        j.penalty = true;
                    count++;
                }
        }
    }
}

//' rcpp_route_times
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::DataFrame rcpp_route_times (const Rcpp::DataFrame graph, bool left_side)
{
    std::unordered_map <std::string, double> x0, y0;
    // these must all be vectors because the IDs are arbitrarily sorted in
    // clockwise order around the junction point.
    std::unordered_map <std::string, std::vector <std::string> >
        out_edges, out_edges_sorted;

    routetimes::fill_edges (graph, x0, y0, out_edges);
    routetimes::sort_edges (out_edges, out_edges_sorted, x0, y0);
    std::vector <Junction> junctions;
    routetimes::replace_junctions (out_edges_sorted, junctions, left_side);

    // dummy values to demonstrate that routetimes::isLess works:
    const size_t len = 6;
    std::vector <OneNode> x (len);
    OneNode n;

    n.x = 1.0;  n.y = -3.0;     n.id = "1";     x [0] = n;
    n.x = 0.0;  n.y = 1.0;      n.id = "2";     x [1] = n;
    n.x = -1.0; n.y = 2.0;      n.id = "3";     x [2] = n;
    n.x = 2.0;  n.y = -0.5;     n.id = "4";     x [3] = n;
    n.x = 3.0;  n.y = 1.5;      n.id = "5";     x [4] = n;
    n.x = -2.5; n.y = -0.7;     n.id = "6";     x [5] = n;

    std::sort (x.begin (), x.end (), routetimes::isLess);

    Rcpp::NumericVector xout (len), yout (len);
    Rcpp::CharacterVector idout (len);
    int count = 0;
    for (auto i: x)
    {
        xout [count] = i.x;
        yout [count] = i.y;
        idout [count++] = i.id;
    }
    Rcpp::DataFrame res = Rcpp::DataFrame::create (
            Rcpp::Named ("x") = xout,
            Rcpp::Named ("y") = yout,
            Rcpp::Named ("id") = idout,
            Rcpp::_["stringsAsFactors"] = false);

    return res;
}
