#include <string>
#include <cmath>

#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

const float INFINITE_FLOAT =  std::numeric_limits<float>::max ();
const double INFINITE_DOUBLE =  std::numeric_limits<double>::max ();
const int INFINITE_INT =  std::numeric_limits<int>::max ();

const std::string osm_p4s = "+proj=longlat +datum=WGS84 +no_defs";

size_t make_edge_name_set (std::unordered_set <std::string> &new_edge_name_set,
        const Rcpp::CharacterVector &new_edges);
void make_edge_name_vec (const size_t n,
        const Rcpp::CharacterVector &new_edges,
        std::vector <std::string> &new_edge_name_vec);
size_t get_edgevec_sizes (const size_t nedges,
        const Rcpp::CharacterVector &new_edges,
        std::vector <size_t> &edgevec_sizes);
Rcpp::NumericVector rcpp_get_bbox_sf (double xmin, double xmax,
        double ymin, double ymax);
Rcpp::List rcpp_aggregate_to_sf (const Rcpp::DataFrame &graph_full,
        const Rcpp::DataFrame &graph_contr, const Rcpp::DataFrame &edge_map);
