#pragma once

#include <vector>

#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

struct OneNode {
    double x, y;
    std::string id;
};

struct OneEdge {
    std::string start, end, id;
    double x0, y0, x1, y1;
    bool penalty;
};

struct Junction {
    std::string in, centre, out;
    bool penalty;
};

namespace routetimes {

bool isLess (OneNode a, OneNode b);

void replace_one_map_edge (
        std::unordered_map <std::string, std::vector <std::string> > &the_edges,
        std::string key, std::string value);

void erase_non_junctions (
        std::unordered_map <std::string, std::vector <std::string> > &the_edges);

void fill_edges (const Rcpp::DataFrame &graph,
        std::unordered_map <std::string, double> &x0,
        std::unordered_map <std::string, double> &y0,
        std::unordered_map <std::string, std::vector <std::string> > &out_edges);

void sort_edges (
        const std::unordered_map <std::string, std::vector <std::string> > &edges_in,
        std::unordered_map <std::string, std::vector <std::string> > &edges_sorted,
        const std::unordered_map <std::string, double> &x0,
        const std::unordered_map <std::string, double> &y0);

void replace_junctions (
        const std::unordered_map <std::string, std::vector <std::string> > &edges,
        std::vector <Junction> junctions,
        bool left_side);

} // end namespace

Rcpp::DataFrame rcpp_route_times (const Rcpp::DataFrame graph);
