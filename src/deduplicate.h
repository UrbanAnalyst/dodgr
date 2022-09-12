#pragma once

#include <vector>
#include <unordered_set>
#include <unordered_map>

#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel,RcppThread)]]
#include <RcppThread.h>
#include <RcppParallel.h>

namespace deduplicate {

typedef std::pair <std::string, std::string> str_pair;
     
struct str_pair_hash
{
    std::size_t operator() (const std::pair <std::string, std::string> &pair) const {
        return std::hash <std::string> () (pair.first) ^ std::hash <std::string> () (pair.second);
    }
};

typedef std::unordered_map <deduplicate::str_pair, double, deduplicate::str_pair_hash> EdgeMapType;

void update_dupl_edge_map (deduplicate::EdgeMapType &edge_map,
        const str_pair &this_pair, const double &val);

} // end namespace deduplicate

Rcpp::DataFrame rcpp_deduplicate (const Rcpp::DataFrame &graph, const std::string fr_col, const std::string to_col,
        const std::string d_col, const std::string t_col);
