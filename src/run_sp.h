#pragma once

#include <memory>
#include <vector>
#include <algorithm> // std::fill, std::reverse
#include <iostream>
#include <fstream>

#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel,RcppThread)]]
#include <RcppThread.h>
#include <RcppParallel.h>

#include "pathfinders.h"

class DGraph;
class PathFinder;

//----------------------------
//----- functions in run_sp.cpp
//----------------------------

namespace run_sp {

std::shared_ptr <HeapDesc> getHeapImpl(const std::string& heap_type);

size_t make_vert_map (const Rcpp::DataFrame &vert_map_in,
        const std::vector <std::string> &vert_map_id,
        const std::vector <size_t> &vert_map_n,
        std::map <std::string, size_t> &vert_map);

void make_vert_to_edge_maps (const std::vector <std::string> &from,
        const std::vector <std::string> &to, const std::vector <double> &wt,
        std::unordered_map <std::string, size_t> &verts_to_edge_map,
        std::unordered_map <std::string, double> &verts_to_dist_map);

size_t get_chunk_size (const size_t nfrom);
} // end namespace run_sp

namespace categorical {

size_t num_edge_types (const std::vector <size_t> &edge_type);

} // end namespace categorical

namespace centrality {

    struct edge_pair_hash {
        inline std::size_t operator()(const std::pair <size_t, size_t> & v) const {
            return v.first * 31 + v.second;
        }
    };

} // end namespace centrality


Rcpp::NumericMatrix rcpp_get_sp_dists (const Rcpp::DataFrame graph,
        const Rcpp::DataFrame vert_map_in,
        Rcpp::IntegerVector fromi,
        Rcpp::IntegerVector toi_in,
        const std::string& heap_type);

Rcpp::NumericMatrix rcpp_get_sp_dists_par (const Rcpp::DataFrame graph,
        const Rcpp::DataFrame vert_map_in,
        Rcpp::IntegerVector fromi,
        Rcpp::IntegerVector toi_in,
        const std::string& heap_type,
        const bool is_spatial);

Rcpp::NumericMatrix rcpp_get_sp_dists_paired_par (const Rcpp::DataFrame graph,
        const Rcpp::DataFrame vert_map_in,
        Rcpp::IntegerVector fromi,
        Rcpp::IntegerVector toi,
        const std::string& heap_type,
        const bool is_spatial);

Rcpp::NumericMatrix rcpp_get_iso (const Rcpp::DataFrame graph,
        const Rcpp::DataFrame vert_map_in,
        Rcpp::IntegerVector fromi,
        Rcpp::NumericVector dlim,
        const std::string& heap_type);

Rcpp::List rcpp_get_paths (const Rcpp::DataFrame graph,
        const Rcpp::DataFrame vert_map_in,
        Rcpp::IntegerVector fromi,
        Rcpp::IntegerVector toi_in,
        const std::string& heap_type);

Rcpp::List rcpp_get_paths_pairwise (const Rcpp::DataFrame graph,
                           const Rcpp::DataFrame vert_map_in,
                           Rcpp::IntegerVector fromi,
                           Rcpp::IntegerVector toi_in,
                           const std::string& heap_type);

// in run_sp_categorical:
Rcpp::NumericMatrix rcpp_get_sp_dists_categorical (const Rcpp::DataFrame graph,
        const Rcpp::DataFrame vert_map_in,
        Rcpp::IntegerVector fromi,
        Rcpp::IntegerVector toi_in,
        const std::string& heap_type,
        const bool proportions_only);

Rcpp::NumericMatrix rcpp_get_sp_dists_cat_threshold (const Rcpp::DataFrame graph,
        const Rcpp::DataFrame vert_map_in,
        Rcpp::IntegerVector fromi,
        const double dlimit,
        const std::string& heap_type);

//----------------------------
//----- functions in centrality.cpp
//----------------------------

Rcpp::NumericVector rcpp_centrality (const Rcpp::DataFrame graph,
        const Rcpp::DataFrame vert_map_in,
        const std::string& heap_type,
        const double dist_threshold,
        const bool edge_centrality,
        const int sample);
