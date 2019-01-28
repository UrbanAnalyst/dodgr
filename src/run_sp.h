#pragma once

#include <memory>
#include <vector>
#include <algorithm> // std::fill, std::reverse
#include <iostream>
#include <fstream>

#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

#include "dijkstra.h"

class DGraph;
class Dijkstra;

constexpr float INFINITE_FLOAT =  std::numeric_limits<float>::max ();
constexpr double INFINITE_DOUBLE =  std::numeric_limits<double>::max ();
constexpr int INFINITE_INT =  std::numeric_limits<int>::max ();

//----------------------------
//----- functions in run_sp.cpp
//----------------------------

namespace run_sp {

std::shared_ptr <HeapDesc> getHeapImpl(const std::string& heap_type);

size_t make_vert_map (const Rcpp::DataFrame &vert_map_in,
        const std::vector <std::string> &vert_map_id,
        const std::vector <unsigned int> &vert_map_n,
        std::map <std::string, unsigned int> &vert_map);

size_t get_fromi_toi (const Rcpp::DataFrame &vert_map_in,
        Rcpp::IntegerVector &fromi, Rcpp::IntegerVector &toi,
        Rcpp::NumericVector &id_vec);

size_t get_fromi (const Rcpp::DataFrame &vert_map_in,
        Rcpp::IntegerVector &fromi, Rcpp::NumericVector &id_vec);

void make_vert_to_edge_maps (const std::vector <std::string> &from,
        const std::vector <std::string> &to, const std::vector <double> &wt,
        std::unordered_map <std::string, unsigned int> &verts_to_edge_map,
        std::unordered_map <std::string, double> &verts_to_dist_map);

} // end namespace run_sp

Rcpp::NumericMatrix rcpp_get_sp_dists (const Rcpp::DataFrame graph,
        const Rcpp::DataFrame vert_map_in,
        Rcpp::IntegerVector fromi,
        Rcpp::IntegerVector toi,
        const std::string& heap_type);

Rcpp::NumericMatrix rcpp_get_sp_dists_par (const Rcpp::DataFrame graph,
        const Rcpp::DataFrame vert_map_in,
        Rcpp::IntegerVector fromi,
        Rcpp::IntegerVector toi,
        const std::string& heap_type);

Rcpp::List rcpp_get_paths (const Rcpp::DataFrame graph,
        const Rcpp::DataFrame vert_map_in,
        Rcpp::IntegerVector fromi,
        Rcpp::IntegerVector toi,
        const std::string& heap_type);

void rcpp_flows_aggregate_par (const Rcpp::DataFrame graph,
        const Rcpp::DataFrame vert_map_in,
        Rcpp::IntegerVector fromi,
        Rcpp::IntegerVector toi,
        Rcpp::NumericMatrix flows,
        const std::string dirtxt,
        const std::string heap_type);

Rcpp::NumericVector rcpp_aggregate_files (const Rcpp::CharacterVector file_names,
        const int len);

Rcpp::NumericVector rcpp_flows_disperse (const Rcpp::DataFrame graph,
        const Rcpp::DataFrame vert_map_in,
        Rcpp::IntegerVector fromi,
        double k,
        Rcpp::NumericMatrix flows,
        std::string heap_type);
