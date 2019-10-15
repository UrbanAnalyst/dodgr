#pragma once

#include <memory>
#include <vector>
#include <algorithm> // std::fill, std::reverse
#include <iostream>
#include <fstream>
#include <deque> // used in centrality

#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

#include "pathfinders.h"
#include "centrality.h"

class DGraph;
class PathFinder;

struct by_w
{
    template <class T1, class T2>
    bool operator () (const std::pair <T1, T2>& lhs,
            const std::pair <T1, T2>& rhs)
    {
        if (lhs.second == rhs.second)
            return (lhs.first > rhs.first);
        else
            return lhs.second > rhs.second;
    }
}; 

Rcpp::NumericVector rcpp_centrality (const Rcpp::DataFrame graph,
        const Rcpp::DataFrame vert_map_in,
        const std::string& heap_type);

