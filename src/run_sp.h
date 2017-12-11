#pragma once

#include <memory>
#include <vector>
#include <algorithm> // std::fill, std::reverse

#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

#include "dijkstra.h"

class DGraph;
class Dijkstra;

const float INFINITE_FLOAT =  std::numeric_limits<float>::max ();
const double INFINITE_DOUBLE =  std::numeric_limits<double>::max ();
const int INFINITE_INT =  std::numeric_limits<int>::max ();


//----------------------------
//----- functions in run_sp.cpp
//----------------------------

// ancilliary functions
std::shared_ptr <HeapDesc> getHeapImpl(const std::string& heap_type);

size_t make_vert_map (Rcpp::DataFrame &vert_map_in,
        std::vector <std::string> &vert_map_id,
        std::vector <unsigned int> &vert_map_n,
        std::map <std::string, unsigned int> &vert_map);

// the main functions
Rcpp::NumericMatrix rcpp_get_sp_dists (Rcpp::DataFrame graph,
        Rcpp::DataFrame vert_map_in,
        std::vector <int> fromi,
        std::vector <int> toi,
        const std::string& heap_type);

Rcpp::NumericMatrix rcpp_get_sp_dists_par (Rcpp::DataFrame graph,
        Rcpp::DataFrame vert_map_in,
        Rcpp::NumericVector fromi,
        Rcpp::NumericVector toi,
        std::string heap_type);

Rcpp::NumericMatrix rcpp_get_sp_dists (Rcpp::DataFrame graph,
        Rcpp::DataFrame vert_map_in,
        std::vector <int> fromi,
        std::vector <int> toi,
        std::string heap_type);

Rcpp::List rcpp_get_paths (Rcpp::DataFrame graph,
        Rcpp::DataFrame vert_map_in,
        std::vector <int> fromi,
        std::vector <int> toi,
        std::string heap_type);

Rcpp::NumericVector rcpp_aggregate_flows (Rcpp::DataFrame graph,
        Rcpp::DataFrame vert_map_in,
        std::vector <int> fromi,
        std::vector <int> toi,
        Rcpp::NumericMatrix flows,
        std::string heap_type);

Rcpp::NumericVector rcpp_aggregate_all_flows (Rcpp::DataFrame graph,
        Rcpp::DataFrame vert_map_in,
        std::vector <int> fromi,
        double k,
        Rcpp::NumericMatrix flows,
        std::string heap_type);

