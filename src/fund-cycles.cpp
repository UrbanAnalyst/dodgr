#include "fund-cycles.h"

//' @noRd
// [[Rcpp::export]]
Rcpp::List rcpp_fundamental_cycles (Rcpp::DataFrame graph,
        Rcpp::DataFrame verts)
{
    std::vector <std::string> vert_id = verts ["id"];
    std::vector <std::string> from = graph ["from"];
    std::vector <std::string> to = graph ["to"];

    std::unordered_map <std::string, size_t> vert_index;
    for (size_t i = 0; i < vert_id.size (); i++)
        vert_index.emplace (vert_id [i], i);

    std::vector <size_t> edge_array (from.size() * 2);
    for (size_t i = 0; i < from.size (); i++)
    {
        edge_array [i * 2] = vert_index.at (from [i]);
        edge_array [i * 2 + 1] = vert_index.at (to [i]);
    }
    graph::Graph <std::string> gr (vert_id, vert_id.size (),
            edge_array, from.size ());
	gr.computeFundamentalCycles();

    Rcpp::List ret (gr.m_fundamentalCycles.size ());
    std::vector<graph::AdjacencyMatrix>::iterator cycle_iter;
    for (cycle_iter = gr.m_fundamentalCycles.begin();
            cycle_iter != gr.m_fundamentalCycles.end (); cycle_iter++)
    {
        graph::Graph<std::string>::NodePath path =
            gr.cycleMatrix2nodePath (*cycle_iter);
        std::vector <std::string> pathi (path.size ());
        size_t i = 0;
        for (const std::string* obj: path)
        {
            pathi [i++] = *obj;
        }
        ret [std::distance (gr.m_fundamentalCycles.begin(), cycle_iter)] = pathi;
    }
    
    return ret;
}
