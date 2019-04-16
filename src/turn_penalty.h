#pragma once

#include "sc-as-network.h"

#include <vector>

#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

struct OneEdge {
    std::string v0, v1, edge;
    double x, y;
};

struct OneCompoundEdge {
    std::string v0, v1, edge0, edge1;
    bool penalty;
};

// ordering function to sort OneEdge structs in clockwise order, for which x and
// y values pre-converted by subtracting the centre values.
// https://stackoverflow.com/questions/6989100/sort-points-in-clockwise-order
struct clockwise_sort
{
    bool operator () (const OneEdge &a, const OneEdge &b)
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
};

typedef std::set <OneEdge, clockwise_sort> RTEdgeSet;

namespace routetimes {

void fill_edges (const Rcpp::DataFrame &graph,
        std::unordered_map <std::string,
                            std::pair <RTEdgeSet, RTEdgeSet> > &the_edges,
        std::unordered_set <std::string> &junction_vertices);

void replace_one_map_edge (
        std::unordered_map <std::string,
                            std::pair <RTEdgeSet, RTEdgeSet> > &the_edges,
        std::string key, OneEdge edge, bool incoming);

void erase_non_junctions (
        std::unordered_map <std::string,
                            std::pair <RTEdgeSet, RTEdgeSet> > &the_edges,
        std::unordered_set <std::string> &junction_vertices);

void replace_junctions (
        const std::unordered_map <std::string,
                                  std::pair <RTEdgeSet, RTEdgeSet> > &the_edges,
        std::vector <OneCompoundEdge> &junctions,
        bool left_side);

Rcpp::DataFrame new_graph (const Rcpp::DataFrame &graph, 
        std::vector <OneCompoundEdge> &junctions, int turn_penalty);

} // end namespace

Rcpp::List rcpp_route_times (const Rcpp::DataFrame graph,
        bool left_side, int turn_penalty);
