#pragma once

#include <Rcpp.h>
#include <algorithm> // std::fill, std::reverse

const float INFINITE_FLOAT =  std::numeric_limits<float>::max ();
const double INFINITE_DOUBLE =  std::numeric_limits<double>::max ();
const int INFINITE_INT =  std::numeric_limits<int>::max ();

template <typename T>
void inst_graph (DGraph *g, unsigned int nedges,
        std::map <std::string, unsigned int> &vert_map,
        std::vector <std::string> &from,
        std::vector <std::string> &to,
        std::vector <T> &dist,
        std::vector <T> &wt);
template void inst_graph <double> (DGraph *g, unsigned int nedges,
        std::map <std::string, unsigned int> &vert_map,
        std::vector <std::string> &from,
        std::vector <std::string> &to,
        std::vector <double> &dist,
        std::vector <double> &wt);


//----------------------------
//----- functions in run_sp.cpp
//----------------------------
Dijkstra * dijkstra_bheap (unsigned int nverts);
Dijkstra * dijkstra_fheap (unsigned int nverts);
Dijkstra * dijkstra_heap23 (unsigned int nverts);
Dijkstra * dijkstra_triheap (unsigned int nverts);
Dijkstra * dijkstra_triheapext (unsigned int nverts);
Dijkstra * dijkstra_radix (unsigned int nverts);

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


//----------------------------
//----- functions in sf_as_network.cpp
//----------------------------
double haversine (double x1, double y1, double x2, double y2);
