#include <string>
#include <cmath>

#include <Rcpp.h>

const float INFINITE_FLOAT =  std::numeric_limits<float>::max ();
const double INFINITE_DOUBLE =  std::numeric_limits<double>::max ();
const int INFINITE_INT =  std::numeric_limits<int>::max ();

const std::string osm_p4s = "+proj=longlat +datum=WGS84 +no_defs";

namespace dodgr_sf {

size_t make_edge_name_set (std::unordered_set <std::string> &new_edge_name_set,
        const Rcpp::CharacterVector &new_edges);
void make_edge_name_vec (const size_t n,
        const Rcpp::CharacterVector &new_edges,
        std::vector <std::string> &new_edge_name_vec);
size_t get_edgevec_sizes (const size_t nedges,
        const Rcpp::CharacterVector &new_edges,
        std::vector <size_t> &edgevec_sizes);
void get_edge_to_vert_maps (const std::vector <size_t> &edgevec_sizes,
        const Rcpp::DataFrame &graph_full,
        const Rcpp::CharacterVector &old_edges,
        const Rcpp::CharacterVector &new_edges,
        const std::vector <std::string> &new_edge_names,
        std::unordered_map <std::string,
                            std::vector <std::string> > &full_from_edge_map,
        std::unordered_map <std::string,
                            std::vector <std::string> > &full_to_edge_map);
void order_vert_sequences (Rcpp::List &edge_sequences,
        std::vector <std::string> &new_edge_names,
        std::unordered_map <std::string,
                            std::vector <std::string> > &full_from_edge_map,
        std::unordered_map <std::string,
                            std::vector <std::string> > &full_to_edge_map);
size_t count_non_contracted_edges (const Rcpp::CharacterVector &contr_edges,
        std::unordered_set <std::string> &new_edge_name_set);
void append_nc_edges (const size_t nc_edge_count, 
        const Rcpp::DataFrame &graph_contr,
        std::unordered_set <std::string> &new_edge_name_set,
        std::vector <std::string> &new_edge_name_vec,
        const Rcpp::List &edge_sequences_contr,
        std::vector <std::string> &all_edge_names,
        Rcpp::List &edge_sequences_all);
void xy_to_sf (const Rcpp::DataFrame &graph_full,
        const Rcpp::List &edge_sequences, 
        const std::vector <std::string> &all_edge_names,
        Rcpp::List &res);

} // end namespace dodgr_sf

Rcpp::NumericVector rcpp_get_bbox_sf (double xmin, double xmax,
        double ymin, double ymax);

Rcpp::List rcpp_aggregate_to_sf (const Rcpp::DataFrame &graph_full,
        const Rcpp::DataFrame &graph_contr, const Rcpp::DataFrame &edge_map);
