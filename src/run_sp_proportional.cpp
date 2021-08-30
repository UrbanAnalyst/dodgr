// Modified functions from run_sp.cpp to calculate proportional distances along
// defined kinds of ways (issue #144)

#include "run_sp.h"

#include "dgraph.h"
#include "heaps/heap_lib.h"
#include <unordered_set>

// # nocov start
template <typename T>
void inst_graph (std::shared_ptr<DGraph> g, size_t nedges,
        const std::map <std::string, size_t>& vert_map,
        const std::vector <std::string>& from,
        const std::vector <std::string>& to,
        const std::vector <T>& dist,
        const std::vector <T>& wt)
{
    for (size_t i = 0; i < nedges; ++i)
    {
        size_t fromi = vert_map.at(from [i]);
        size_t toi = vert_map.at(to [i]);
        g->addNewEdge (fromi, toi, dist [i], wt [i], i);
    }
}
// # nocov end

struct OneProportionalDist : public RcppParallel::Worker
{
    RcppParallel::RVector <int> dp_fromi;
    const std::vector <size_t> toi;
    const std::vector <int> edge_type;
    const size_t nverts;
    const std::vector <double> vx;
    const std::vector <double> vy;
    const std::shared_ptr <DGraph> g;
    const std::string heap_type;
    bool is_spatial;

    RcppParallel::RMatrix <double> dout;

    // constructor
    OneProportionalDist (
            const RcppParallel::RVector <int> fromi,
            const std::vector <size_t> toi_in,
            const std::vector <int> edge_type_in,
            const size_t nverts_in,
            const std::vector <double> vx_in,
            const std::vector <double> vy_in,
            const std::shared_ptr <DGraph> g_in,
            const std::string & heap_type_in,
            const bool & is_spatial_in,
            RcppParallel::RMatrix <double> dout_in) :
        dp_fromi (fromi), toi (toi_in), edge_type (edge_type_in),
        nverts (nverts_in), vx (vx_in), vy (vy_in),
        g (g_in), heap_type (heap_type_in), is_spatial (is_spatial_in),
        dout (dout_in)
    {
    }

    // Parallel function operator
    void operator() (std::size_t begin, std::size_t end)
    {
        std::shared_ptr<PF::PathFinder> pathfinder =
            std::make_shared <PF::PathFinder> (nverts,
                    *run_sp::getHeapImpl (heap_type), g);
        std::vector <double> w (nverts);
        std::vector <double> d (nverts);
        std::vector <long int> prev (nverts);

        std::vector <double> heuristic (nverts, 0.0);

        for (std::size_t i = begin; i < end; i++)
        {
            size_t from_i = static_cast <size_t> (dp_fromi [i]);

            if (is_spatial)
            {
                for (size_t j = 0; j < nverts; j++)
                {
                    const double dx = vx [j] - vx [from_i],
                        dy = vy [j] - vy [from_i];
                    heuristic [j] = sqrt (dx * dx + dy * dy);
                }
                pathfinder->AStar (d, w, prev, heuristic, from_i, toi);
            } else if (heap_type.find ("set") == std::string::npos)
                pathfinder->Dijkstra (d, w, prev, from_i, toi);
            else
                pathfinder->Dijkstra_set (d, w, prev, from_i);

            for (size_t j = 0; j < toi.size (); j++)
            {
                if (w [toi [j]] < INFINITE_DOUBLE)
                {
                    dout (i, j) = d [toi [j]];
                }
            }
        }
    }
                                   
};

size_t proportional::num_edge_types (const std::vector <int> &edge_type)
{
    std::unordered_set <int> type_set;
    for (auto e: edge_type)
        type_set.emplace (e);

    return type_set.size ();
}


//' rcpp_get_sp_dists_par
//'
//' Implemented in parallal form only; no single-threaded version
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_get_sp_dists_proportional (const Rcpp::DataFrame graph,
        const Rcpp::DataFrame vert_map_in,
        Rcpp::IntegerVector fromi,
        Rcpp::IntegerVector toi_in,
        const std::string& heap_type,
        const bool is_spatial)
{
    std::vector <size_t> toi =
        Rcpp::as <std::vector <size_t> > ( toi_in);

    size_t nfrom = static_cast <size_t> (fromi.size ());
    size_t nto = static_cast <size_t> (toi.size ());

    const std::vector <std::string> from = graph ["from"];
    const std::vector <std::string> to = graph ["to"];
    const std::vector <double> dist = graph ["d"];
    const std::vector <double> wt = graph ["d_weighted"];
    const std::vector <int> edge_type = graph ["edge_type"];

    const size_t num_types = proportional::num_edge_types (edge_type);

    const size_t nedges = static_cast <size_t> (graph.nrow ());
    std::map <std::string, size_t> vert_map;
    std::vector <std::string> vert_map_id = vert_map_in ["vert"];
    std::vector <size_t> vert_map_n = vert_map_in ["id"];
    const size_t nverts = run_sp::make_vert_map (vert_map_in, vert_map_id,
            vert_map_n, vert_map);

    std::vector <double> vx (nverts), vy (nverts);
    if (is_spatial)
    {
        vx = Rcpp::as <std::vector <double> > (vert_map_in ["x"]);
        vy = Rcpp::as <std::vector <double> > (vert_map_in ["y"]);
    }

    std::shared_ptr <DGraph> g = std::make_shared <DGraph> (nverts);
    inst_graph (g, nedges, vert_map, from, to, dist, wt);

    Rcpp::NumericVector na_vec = Rcpp::NumericVector (nfrom * nto,
            Rcpp::NumericVector::get_na ());
    Rcpp::NumericMatrix dout (static_cast <int> (nfrom),
            static_cast <int> (nto), na_vec.begin ());

    // Create parallel worker
    OneProportionalDist one_dist (RcppParallel::RVector <int> (fromi), toi,
            edge_type, nverts, vx, vy,
            g, heap_type, is_spatial,
            RcppParallel::RMatrix <double> (dout));

    size_t chunk_size = run_sp::get_chunk_size (nfrom);
    RcppParallel::parallelFor (0, nfrom, one_dist, chunk_size);
    
    return (dout);
}
