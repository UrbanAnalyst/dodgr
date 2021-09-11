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
        const std::vector <size_t>& edge_type,
        const std::vector <T>& dist,
        const std::vector <T>& wt)
{
    for (size_t i = 0; i < nedges; ++i)
    {
        size_t fromi = vert_map.at(from [i]);
        size_t toi = vert_map.at(to [i]);
        g->addNewEdge (fromi, toi, dist [i], wt [i], edge_type [i]);
    }
}
// # nocov end

struct OneProportionalDist : public RcppParallel::Worker
{
    RcppParallel::RVector <int> dp_fromi;
    const std::vector <size_t> toi;
    const std::vector <size_t> edge_type;
    const size_t nverts;
    const std::vector <double> vx;
    const std::vector <double> vy;
    const std::shared_ptr <DGraph> g;
    const std::string heap_type;
    const size_t num_edge_types;

    RcppParallel::RMatrix <double> dout;

    // constructor
    OneProportionalDist (
            const RcppParallel::RVector <int> fromi,
            const std::vector <size_t> toi_in,
            const std::vector <size_t> edge_type_in,
            const size_t nverts_in,
            const std::vector <double> vx_in,
            const std::vector <double> vy_in,
            const std::shared_ptr <DGraph> g_in,
            const std::string & heap_type_in,
            const size_t & num_edge_types_in,
            RcppParallel::RMatrix <double> dout_in) :
        dp_fromi (fromi), toi (toi_in), edge_type (edge_type_in),
        nverts (nverts_in), vx (vx_in), vy (vy_in),
        g (g_in), heap_type (heap_type_in), 
        num_edge_types (num_edge_types_in),
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
        std::vector <double> d (nverts * num_edge_types);
        std::vector <long int> prev (nverts);

        std::vector <double> heuristic (nverts, 0.0);

        const size_t nto = toi.size ();

        for (std::size_t i = begin; i < end; i++)
        {
            size_t from_i = static_cast <size_t> (dp_fromi [i]);

            // only implemented for spatial graphs
            for (size_t j = 0; j < nverts; j++)
            {
                const double dx = vx [j] - vx [from_i],
                    dy = vy [j] - vy [from_i];
                heuristic [j] = sqrt (dx * dx + dy * dy);
            }
            pathfinder->AStarEdgeType (d, w, prev, heuristic, from_i, toi);

            for (size_t j = 0; j < toi.size (); j++)
            {
                if (w [toi [j]] < INFINITE_DOUBLE)
                {
                    for (size_t k = 0; k < num_edge_types; k++)
                    {
                        const double dto = d [toi [j] + k * nverts];
                        if (dto < INFINITE_DOUBLE)
                            dout (i, j + k * nto) = dto;
                    }
                }
            }
        }
    }
                                   
};

size_t proportional::num_edge_types (const std::vector <size_t> &edge_type)
{
    std::unordered_set <size_t> type_set;
    for (auto e: edge_type)
        type_set.emplace (e);

    return type_set.size ();
}

void PF::PathFinder::AStarEdgeType (std::vector<double>& d,
        std::vector<double>& w,
        std::vector<long int>& prev,
        const std::vector<double>& heur,
        const size_t v0,
        const std::vector <size_t> &to_index)
{
    // d is already nverts * num_edge_types, and init_arrays only fills without
    // resizing.
    const DGraphEdge *edge;

    const size_t n = m_graph->nVertices();
    const std::vector<DGraphVertex>& vertices = m_graph->vertices();

    PF::PathFinder::init_arrays (d, w, prev, m_open, m_closed, v0, n);
    m_heap->insert (v0, heur [v0]);

    size_t n_reached = 0;
    const size_t n_targets = to_index.size ();
    bool *is_target = new bool [n];
    std::fill (is_target, is_target + n, false);
    for (auto t: to_index)
        is_target [t] = true;

    while (m_heap->nItems() > 0) {
        size_t v = m_heap->deleteMin();

        m_closed [v] = true;
        m_open [v] = false;

        edge = vertices [v].outHead;
        scan_edge_types_heur (edge, d, w, prev, m_open, m_closed, v, heur);

        if (is_target [v])
            n_reached++;
        if (n_reached == n_targets)
            break;
    } // end while m_heap->nItems
    delete [] is_target;
}


void PF::PathFinder::scan_edge_types_heur (const DGraphEdge *edge,
        std::vector<double>& d,
        std::vector<double>& w,
        std::vector<long int>& prev,
        bool *m_open_vec,
        const bool *m_closed_vec,
        const size_t &v0,
        const std::vector<double> &heur)    // heuristic for A*
{
    const size_t n = w.size ();

    while (edge) {

        const size_t et = edge->target;
        const size_t edge_id = edge->edge_id;

        if (!m_closed_vec [et])
        {
            double wt = w [v0] + edge->wt;
            if (wt < w [et]) {
                d [et] = d [v0] + edge->dist;
                if (edge_id > 0L) {
                    if (d [et + edge_id * n] == INFINITE_DOUBLE)
                        d [et + edge_id * n] = edge->dist;
                    else
                        d [et + edge_id * n] += edge->dist;
                }
                w [et] = wt;
                prev [et] = static_cast <int> (v0);

                if (m_open_vec [et]) {
                    m_heap->decreaseKey(et, wt + heur [et] - heur [v0]);
                }
                else {
                    m_heap->insert (et, wt + heur [et] - heur [v0]);
                    m_open_vec [et] = true;
                }
            } else
                m_closed [et] = true;
        }
        edge = edge->nextOut;
    }
}

//' rcpp_get_sp_dists_proportional
//'
//' The `graph` must have an `edge_type` column of non-negative integers,
//' with 0 denoting edges which are not aggregated, and all other values
//' defining aggregation categories.
//'
//' Implemented in parallal form only; no single-threaded version, and
//' only for AStar (so graphs must be spatial).
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_get_sp_dists_proportional (const Rcpp::DataFrame graph,
        const Rcpp::DataFrame vert_map_in,
        Rcpp::IntegerVector fromi,
        Rcpp::IntegerVector toi_in,
        const std::string& heap_type)
{
    std::vector <size_t> toi =
        Rcpp::as <std::vector <size_t> > ( toi_in);

    size_t nfrom = static_cast <size_t> (fromi.size ());
    size_t nto = static_cast <size_t> (toi.size ());

    const std::vector <std::string> from = graph ["from"];
    const std::vector <std::string> to = graph ["to"];
    const std::vector <double> dist = graph ["d"];
    const std::vector <double> wt = graph ["d_weighted"];
    const std::vector <size_t> edge_type = graph ["edge_type"];

    const size_t num_types = proportional::num_edge_types (edge_type);

    const size_t nedges = static_cast <size_t> (graph.nrow ());
    std::map <std::string, size_t> vert_map;
    std::vector <std::string> vert_map_id = vert_map_in ["vert"];
    std::vector <size_t> vert_map_n = vert_map_in ["id"];
    const size_t nverts = run_sp::make_vert_map (vert_map_in, vert_map_id,
            vert_map_n, vert_map);

    std::vector <double> vx (nverts), vy (nverts);
    vx = Rcpp::as <std::vector <double> > (vert_map_in ["x"]);
    vy = Rcpp::as <std::vector <double> > (vert_map_in ["y"]);

    std::shared_ptr <DGraph> g = std::make_shared <DGraph> (nverts);
    inst_graph (g, nedges, vert_map, from, to, edge_type, dist, wt);

    // One standard [nfrom, nto] matrix plus one extra for each value of edge_type:
    Rcpp::NumericVector na_vec = Rcpp::NumericVector (nfrom * nto * num_types,
            Rcpp::NumericVector::get_na ());
    Rcpp::NumericMatrix dout (static_cast <int> (nfrom),
            static_cast <int> (nto * num_types), na_vec.begin ());

    // Create parallel worker
    OneProportionalDist one_dist (RcppParallel::RVector <int> (fromi), toi,
            edge_type, nverts, vx, vy,
            g, heap_type, num_types,
            RcppParallel::RMatrix <double> (dout));

    size_t chunk_size = run_sp::get_chunk_size (nfrom);
    RcppParallel::parallelFor (0, nfrom, one_dist, chunk_size);
    
    return (dout);
}
