#include "run_sp.h"
#include "pathfinders.h"
#include "heaps/heap.h"

#include <algorithm> // std::fill
#include <fstream> // file output for parallel jobs
#include <unordered_set>

/*************************************************************************
 * Direct implementation of
 * "A Faster Algorithm for Betweenness Centrality", Ulrik Brandes (2001)
 * Journal of Mathematical Sociology 25(2):163-177
 * - same algorithm as used in igraph and networkx
 *************************************************************************/

const double epsilon = DBL_MIN; // edge weight comparison == 0
// see
// https://github.com/igraph/igraph/blob/96c2cc751063bf3ba7e920e99793956013cef6b5/include/igraph_nongraph.h#L41
// which defines epsilon as 1e-10. That value is passed to `igraph_cmp_epsilon`,
// which is defined at
// https://github.com/igraph/igraph/blob/96c2cc751063bf3ba7e920e99793956013cef6b5/src/math/utils.c#L108
// and uses a comparison of fabs (diff) < (eps * DBL_MIN),
// where DBL_MIN is the standard lower, non-zero limit
// https://en.cppreference.com/w/cpp/types/numeric_limits/min
// of generally around 1e-308, so (eps * DBL_MIN) is then sub-normal. This code
// just uses an epsilon equal to DBL_MIN.

// # nocov start
template <typename T>
void inst_graph (std::shared_ptr<DGraph> g, size_t nedges,
        const std::map <std::string, size_t>& vert_map,
        const std::vector <std::string>& from,
        const std::vector <std::string>& to,
        const std::vector <T>& dist,
        const std::vector <T>& wt)
{

    std::unordered_set < std::pair <size_t, size_t>, centrality::edge_pair_hash > edge_set;

    for (size_t i = 0; i < nedges; ++i)
    {
        size_t fromi = vert_map.at(from [i]);
        size_t toi = vert_map.at(to [i]);

        std::pair <size_t, size_t> edge_pair {fromi, toi};
        if (edge_set.find (edge_pair) == edge_set.end ())
        {
            edge_set.emplace (edge_pair);
            g->addNewEdge (fromi, toi, dist [i], wt [i], i);
        }

    }
}
// # nocov end

struct OneCentralityVert : public RcppParallel::Worker
{
    size_t nverts; // can't be const because of reinterpret case
    const std::string heap_type;
    const std::vector <double> vert_wts;
    const double dist_threshold;
    std::shared_ptr <DGraph> g;

    std::vector <double> output;

    // Constructor 1: The main constructor
    OneCentralityVert (
            const size_t nverts_in,
            const std::string heap_type_in,
            const std::vector <double> vert_wts_in,
            const double dist_threshold_in,
            const std::shared_ptr <DGraph> g_in) :
        nverts (nverts_in), heap_type (heap_type_in), 
        vert_wts (vert_wts_in), dist_threshold (dist_threshold_in),
        g (g_in), output ()
    {
        output.resize (nverts, 0.0);
    }

    // Constructor 2: The Split constructor
    OneCentralityVert (
            const OneCentralityVert &oneCentralityVert,
            RcppParallel::Split) :
        nverts (oneCentralityVert.nverts),
        heap_type (oneCentralityVert.heap_type), 
        vert_wts (oneCentralityVert.vert_wts),
        dist_threshold (oneCentralityVert.dist_threshold),
        g (oneCentralityVert.g), output ()
    {
        output.resize (nverts, 0.0);
    }


    // Parallel function operator
    void operator() (size_t begin, size_t end)
    {
        std::shared_ptr<PF::PathFinder> pathfinder =
            std::make_shared <PF::PathFinder> (nverts,
                    *run_sp::getHeapImpl (heap_type), g);

        std::vector <double> cent (nverts, 0.0);

        for (size_t v = begin; v < end; v++)
        {
            if (RcppThread::isInterrupted (v % static_cast<int>(100) == 0))
                return;
            const double vwt = vert_wts [v];
            pathfinder->Centrality_vertex (cent,
                    static_cast <size_t> (v),
                    vwt, dist_threshold);
        }

        for (size_t i = 0; i < nverts; i++)
            output [i] += cent [i];
    }

    void join (const OneCentralityVert &rhs)
    {
        for (size_t i = 0; i < output.size (); i++)
            output [i] += rhs.output [i];
    }
};

struct OneCentralityEdge : public RcppParallel::Worker
{
    size_t nverts; // can't be const because of reinterpret case
    size_t nedges;
    const std::string heap_type;
    const std::vector <double> vert_wts;
    const double dist_threshold;
    std::shared_ptr <DGraph> g;

    std::vector <double> output;

    // Constructor 1: The main constructor
    OneCentralityEdge (
            const size_t nverts_in,
            const size_t nedges_in,
            const std::string heap_type_in,
            const std::vector <double> vert_wts_in,
            const double dist_threshold_in,
            const std::shared_ptr <DGraph> g_in) :
        nverts (nverts_in), nedges (nedges_in), 
        heap_type (heap_type_in), vert_wts (vert_wts_in),
        dist_threshold (dist_threshold_in),
        g (g_in), output ()
    {
        output.resize (nedges, 0.0);
    }

    // Constructor 2: The Split constructor
    OneCentralityEdge (
            const OneCentralityEdge &oneCentralityEdge,
            RcppParallel::Split) :
        nverts (oneCentralityEdge.nverts),
        nedges (oneCentralityEdge.nedges), 
        heap_type (oneCentralityEdge.heap_type),
        vert_wts (oneCentralityEdge.vert_wts),
        dist_threshold (oneCentralityEdge.dist_threshold),
        g (oneCentralityEdge.g),
        output ()
    {
        output.resize (nedges, 0.0);
    }

    // Parallel function operator
    void operator() (size_t begin, size_t end)
    {
        std::shared_ptr<PF::PathFinder> pathfinder =
            std::make_shared <PF::PathFinder> (nverts,
                    *run_sp::getHeapImpl (heap_type), g);

        std::vector <double> cent (nedges, 0.0);

        for (size_t v = begin; v < end; v++)
        {
            if (RcppThread::isInterrupted (v % static_cast<int>(100) == 0))
                return;
            const double vwt = vert_wts [v];
            pathfinder->Centrality_edge (cent, v, vwt, nedges, dist_threshold);
        }
        for (size_t i = 0; i < nedges; i++)
            output [i] += cent [i];
    }

    void join (const OneCentralityEdge &rhs)
    {
        for (size_t i = 0; i < output.size (); i++)
            output [i] += rhs.output [i];
    }
};

void PF::PathFinder::Centrality_vertex (
        std::vector <double>& cent,
        const size_t s,
        const double vert_wt,
        const double dist_threshold)
{
    const DGraphEdge *edge;

    const size_t n = m_graph->nVertices();
    const std::vector <DGraphVertex>& vertices = m_graph->vertices();

    std::deque <size_t> v_stack;

    std::vector <double> w (n, 0.0);
    w [s] = 1.0;

    m_heap->insert (s, -1.0);

    // sigma as long int because graphs can be bigger than INT_MAX
    std::vector <long int> sigma (n, 0);
    sigma [s] = 1L;

    std::vector <std::vector <size_t> > prev_vert (n);

    while (m_heap->nItems() > 0) {
        size_t v = m_heap->deleteMin();

        if (w [v] > dist_threshold)
            continue;

        v_stack.push_back (v);

        edge = vertices [v].outHead;

        // NOTE: This 'target_set', and the one below in 'centrality_edge', are
        // only needed for graphs which are submitted with (potentially)
        // duplicate edges. The set is used to ensure centrality values are only
        // aggregated once along any set of duplicated edges. The effect of not
        // doing this is documented at
        // https://github.com/ATFutures/dodgr/issues/186.
        //
        // Nevertheless, the other commits flagged in that issue add a function
        // to the internal 'inst_graph' function at the top of this file ensures
        // no duplicated edges are added, so this 'target_set' is not actually
        // needed. The issue shows that it does not decrease computational
        // efficiency here, so it's left for the moment, but can easily be
        // removed later.
        std::unordered_set <size_t> target_set;

        while (edge) {

            size_t et = edge->target;
            double wt = w [v] + edge->wt;

            if (target_set.find (et) == target_set.end ())
            {
                target_set.emplace (et);

                if (w [et] == 0.0) // first connection to et
                {
                    prev_vert [et] = std::vector <size_t> (1L, v);

                    sigma [et] = sigma [v];
                    w [et] = wt;
                    m_heap->insert (et, wt);
                }  else if (wt < w [et])
                {
                    prev_vert [et] = std::vector <size_t> (1L, v);

                    sigma [et] = sigma [v];
                    w [et] = wt;
                    m_heap->decreaseKey (et, wt);
                } else if (fabs (wt - w [et]) < epsilon)
                {
                    std::vector <size_t> vert_vec = prev_vert [et];
                    vert_vec.push_back (v);
                    prev_vert [et] = vert_vec;

                    sigma [et] += sigma [v];
                }
            }

            edge = edge->nextOut;
        }
    } // end while nItems > 0

    // Then read from the stack and count centrality paths
    std::vector <double> delta (n, 0.0);
    while (!v_stack.empty ())
    {
        const size_t v = v_stack.back ();
        v_stack.pop_back ();
        std::vector <size_t> vert_vec = prev_vert [v];
        double tempd = (1.0 + delta [v]) / static_cast <double> (sigma [v]);
        for (auto ws: vert_vec)
        {
            delta [ws] += static_cast <double> (sigma [ws]) * tempd;
        }
        if (v != s)
            cent [v] += vert_wt * delta [v];
    }
}


void PF::PathFinder::Centrality_edge (
        std::vector <double>& cent,
        const size_t s,
        const double vert_wt,
        const size_t nedges,
        const double dist_threshold)
{
    const DGraphEdge *edge;

    const size_t n = m_graph->nVertices();
    const std::vector <DGraphVertex>& vertices = m_graph->vertices();

    std::deque <size_t> v_stack;

    std::vector <double> w (n, 0.0);
    w [s] = 1.0;

    m_heap->insert (s, -1.0);

    std::vector <long int> sigma (n, 0);
    sigma [s] = 1L;
    std::vector <long int> sigma_edge (nedges, 0);

    std::vector <std::vector <size_t> > prev_vert (n), prev_edge (n);

    while (m_heap->nItems() > 0) {

        size_t v = m_heap->deleteMin();

        if (w [v] > dist_threshold)
            continue;

        v_stack.push_back (v);

        edge = vertices [v].outHead;

        // See comment in 'centrality_vertex', above, about this 'target_set'.
        std::unordered_set <size_t> target_set;

        while (edge) {

            size_t et = edge->target;
            double wt = w [v] + edge->wt;

            if (target_set.find (et) == target_set.end ())
            {
                target_set.emplace (et);

                // DGraph has no edge iterator, so prev_edge elements contains
                // pairwise elements of [from vertex, edge_id]

                if (w [et] == 0.0) // first connection to et
                {
                    prev_edge [et] = std::vector <size_t> {v, edge->edge_id};

                    sigma [et] = sigma [v];
                    sigma_edge [edge->edge_id] = sigma [v];

                    w [et] = wt;
                    m_heap->insert (et, wt);
                }  else if (wt < w [et])
                {
                    prev_edge [et] = std::vector <size_t> {v, edge->edge_id};

                    sigma [et] = sigma [v];
                    sigma_edge [edge->edge_id] = sigma [v];

                    w [et] = wt;
                    m_heap->decreaseKey (et, wt);
                } else if (fabs (wt - w [et]) < epsilon)
                {
                    std::vector <size_t> edge_vec = prev_edge [et];
                    edge_vec.push_back (v);
                    edge_vec.push_back (edge->edge_id);
                    prev_edge [et] = edge_vec;

                    sigma [et] += sigma [v];
                    sigma_edge [edge->edge_id] += sigma [v];
                }
            }

            edge = edge->nextOut;
        }
    } // end while nItems > 0

    // Then read from the stack and count centrality paths
    std::vector <double> delta (n, 0.0);
    while (!v_stack.empty ())
    {
        const size_t v = v_stack.back ();
        v_stack.pop_back ();
        std::vector <size_t> edge_vec = prev_edge [v];
        double tempd = (1.0 + delta [v]) / static_cast <double> (sigma [v]);

        std::vector <size_t>::iterator it = edge_vec.begin ();
        // The dereferenced edge_vec iterator is simply a direct index
        while (it != edge_vec.end ())
        {
            delta [*it] += static_cast <double> (sigma [*it]) * tempd;
            it = std::next (it);
            cent [*it] += static_cast <double> (sigma_edge [*it]) * tempd * vert_wt;
            it = std::next (it);
        }
    }
}


//' rcpp_centrality - parallel function
//'
//' sample is used to estimate timing, by calculating centrality from just a few
//' vertices.
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_centrality (const Rcpp::DataFrame graph,
        const Rcpp::DataFrame vert_map_in,
        const std::string& heap_type,
        const double dist_threshold,
        const bool edge_centrality,
        const int sample)
{
    std::vector <std::string> from = graph ["from"];
    std::vector <std::string> to = graph ["to"];
    std::vector <double> dist = graph ["d"];
    std::vector <double> wt = graph ["d_weighted"];

    const size_t nedges = static_cast <size_t> (graph.nrow ());
    std::map <std::string, size_t> vert_map;
    std::vector <std::string> vert_map_id = vert_map_in ["vert"];
    std::vector <size_t> vert_map_n = vert_map_in ["id"];
    size_t nverts = run_sp::make_vert_map (vert_map_in, vert_map_id,
            vert_map_n, vert_map);

    Rcpp::CharacterVector v_nms = vert_map_in.attr ("names");
    std::vector <double> vert_wts (static_cast <size_t> (vert_map_in.nrow ()), 1.0);
    for (auto n: v_nms) {
        if (n == "vert_wts") {
            std::vector <double> tempd = vert_map_in ["vert_wts"];
            std::copy (tempd.begin (), tempd.end (), vert_wts.begin ());
            break;
        }
    }

    std::shared_ptr <DGraph> g = std::make_shared <DGraph> (nverts);
    inst_graph (g, nedges, vert_map, from, to, dist, wt);

    size_t nverts_to_use = nverts;
    if (sample > 0)
        nverts_to_use = static_cast <size_t> (sample);

    std::vector <double> result;
    if (edge_centrality)
    {
        OneCentralityEdge one_centrality (nverts, nedges, heap_type,
                vert_wts, dist_threshold, g);

        RcppParallel::parallelReduce (0, nverts_to_use, one_centrality);
        result = one_centrality.output;
    } else // vertex centrality
    {
        OneCentralityVert one_centrality (nverts, heap_type, vert_wts,
                dist_threshold, g);

        RcppParallel::parallelReduce (0, nverts_to_use, one_centrality);
        result = one_centrality.output;
    }

    return Rcpp::wrap (result);
}
