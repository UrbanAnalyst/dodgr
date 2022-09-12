
#include "deduplicate.h"

void deduplicate::update_dupl_edge_map (deduplicate::EdgeMapType &edge_map,
        const deduplicate::str_pair &this_pair, const double &val)
{
    if (edge_map.find (this_pair) == edge_map.end ())
    {
        edge_map.emplace (this_pair, val);
    } else
    {
        double val_min = edge_map.find (this_pair)->second;
        if (val < val_min)
        {
            edge_map.erase (this_pair);
            edge_map.emplace (this_pair, val);
        }
    }
}

//' De-duplicate edges by replacing with minimal weighted distances and times
//'
//' @param graph Full graph in any form
//' @param fr_col Name of column holding from edge labels
//' @param to_col Name of column holding to edge labels
//' @param d_col Name of column holding weighted distances
//' @param t_col Name of column holding weighted times (or "" if no weighted
//' times).
//' @return A `data.frame` of 3 columns: 'from', 'to', and 'd', where 'd' is the
//' minimal value taken from all duplicated edges. If 't_col' is specified, the
//' equivalent minimal times are in the lower half of the result.
//' @noRd
// [[Rcpp::export]]
Rcpp::DataFrame rcpp_deduplicate (const Rcpp::DataFrame &graph, const std::string fr_col, const std::string to_col,
        const std::string d_col, const std::string t_col)
{
    const bool has_time = t_col.length () > 0L;

    std::unordered_set < deduplicate::str_pair, deduplicate::str_pair_hash > edge_pair_set;
    std::unordered_set < deduplicate::str_pair, deduplicate::str_pair_hash > edge_pair_dupl;

    const std::vector <std::string> fr = graph [fr_col];
    const std::vector <std::string> to = graph [to_col];
    const std::vector <double> d = graph [d_col];
    const std::vector <double> t = graph [t_col];
    
    const size_t n = static_cast <size_t> (graph.nrow ());

    // First loop to collect duplicated edges
    for (size_t i = 0; i < n; i++)
    {
        deduplicate::str_pair this_pair {fr [i], to [i]};
        if (edge_pair_set.find (this_pair) != edge_pair_set.end ())
        {
            edge_pair_dupl.emplace (this_pair);
        }
        edge_pair_set.emplace (this_pair);
    }

    // Then aggregate values for each duplicate
    deduplicate::EdgeMapType edge_map_d, edge_map_t;

    for (size_t i = 0; i < n; i++)
    {
        deduplicate::str_pair this_pair {fr [i], to [i]};

        if (edge_pair_dupl.find (this_pair) == edge_pair_dupl.end ())
        {
            continue;
        }

        deduplicate::update_dupl_edge_map (edge_map_d, this_pair, d [i]);

        if (has_time)
        {
            deduplicate::update_dupl_edge_map (edge_map_t, this_pair, t [i]);
        }

    }

    const size_t ndupl = edge_map_d.size ();
    std::vector <std::string> fr_dupl (ndupl), to_dupl (ndupl);
    std::vector <double> d_dupl (ndupl);

    size_t i = 0;
    for (auto e: edge_map_d)
    {
        fr_dupl [i] = e.first.first;
        to_dupl [i] = e.first.second;
        d_dupl [i++] = e.second;
    }
    if (has_time)
    {
        fr_dupl.resize (ndupl * 2);
        to_dupl.resize (ndupl * 2);
        d_dupl.resize (ndupl * 2);
        for (auto e: edge_map_t)
        {
            fr_dupl [i] = e.first.first;
            to_dupl [i] = e.first.second;
            d_dupl [i++] = e.second;
        }
    }

    Rcpp::DataFrame res = Rcpp::DataFrame::create (
            Rcpp::Named ("from") = fr_dupl,
            Rcpp::Named ("to") = to_dupl,
            Rcpp::Named ("d") = d_dupl,
            Rcpp::_["stringsAsFactors"] = false);

    return res;
}
