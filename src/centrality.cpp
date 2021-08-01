#include "run_sp.h"
#include "pathfinders.h"
#include "heaps/heap.h"

#include <algorithm> // std::fill
#include <fstream> // file output for parallel jobs

/*************************************************************************
 * Direct implementation of
 * "A Faster Algorithm for Betweenness Centrality", Ulrik Brandes (2001)
 * Journal of Mathematical Sociology 25(2):163-177
 * - same algorithm as used in igraph and networkx
 *************************************************************************/

const double epsilon = 1.0e-10; // edge weight comparison == 0
// see https://github.com/igraph/igraph/blob/master/src/igraph_math.h#L49

// # nocov start
template <typename T>
void inst_graph (std::shared_ptr<DGraph> g, unsigned int nedges,
        const std::map <std::string, unsigned int>& vert_map,
        const std::vector <std::string>& from,
        const std::vector <std::string>& to,
        const std::vector <T>& dist,
        const std::vector <T>& wt)
{
    for (unsigned int i = 0; i < nedges; ++i)
    {
        unsigned int fromi = vert_map.at(from [i]);
        unsigned int toi = vert_map.at(to [i]);
        g->addNewEdge (fromi, toi, dist [i], wt [i], i);
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
                    static_cast <unsigned int> (v),
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
        const unsigned int s,
        const double vert_wt,
        const double dist_threshold)
{
    const DGraphEdge *edge;

    const unsigned int n = m_graph->nVertices();
    const std::vector <DGraphVertex>& vertices = m_graph->vertices();

    std::deque <unsigned int> v_stack;

    std::vector <double> w (n, 0.0);
    w [s] = 1.0;

    m_heap->insert (s, -1.0);

    std::vector <int> sigma (n, 0);
    sigma [s] = 1;

    std::vector <std::vector <unsigned int> > prev_vert (n);

    while (m_heap->nItems() > 0) {
        unsigned int v = m_heap->deleteMin();

        if (w [v] > dist_threshold)
            continue;

        v_stack.push_back (v);

        edge = vertices [v].outHead;
        while (edge) {

            unsigned int et = edge->target;
            double wt = w [v] + edge->wt;

            std::vector <unsigned int> vert_vec;

            if (w [et] == 0.0) // first connection to et
            {
                vert_vec.resize (1);
                vert_vec [0] = v;
                prev_vert [et] = vert_vec;

                sigma [et] = sigma [v];
                w [et] = wt;
                m_heap->insert (et, wt);
            }  else if (wt < w [et])
            {
                vert_vec.resize (1);
                vert_vec [0] = v;
                prev_vert [et] = vert_vec;

                sigma [et] = sigma [v];
                w [et] = wt;
                m_heap->decreaseKey (et, wt);
            } else if (fabs (wt - w [et]) < epsilon)
            {
                vert_vec = prev_vert [et];
                vert_vec.resize (vert_vec.size () + 1);
                vert_vec [vert_vec.size () - 1] = v;
                prev_vert [et] = vert_vec;

                sigma [et] += sigma [v];
            }

            edge = edge->nextOut;
        }
    } // end while nItems > 0

    // Then read from the stack and count centrality paths
    std::vector <double> delta (n, 0.0);
    while (!v_stack.empty ())
    {
        const unsigned int v = v_stack.back ();
        v_stack.pop_back ();
        std::vector <unsigned int> vert_vec = prev_vert [v];
        double tempd = (1.0 + delta [v]) / sigma [v];
        for (auto ws: vert_vec)
        {
            delta [ws] += sigma [ws] * tempd;
        }
        if (v != s)
            cent [v] += vert_wt * delta [v];
    }
}


void PF::PathFinder::Centrality_edge (
        std::vector <double>& cent,
        const unsigned int s,
        const double vert_wt,
        const unsigned int nedges,
        const double dist_threshold)
{
    const DGraphEdge *edge;

    const unsigned int n = m_graph->nVertices();
    const std::vector <DGraphVertex>& vertices = m_graph->vertices();

    std::deque <unsigned int> v_stack;

    std::vector <double> w (n, 0.0);
    w [s] = 1.0;

    m_heap->insert (s, -1.0);

    std::vector <int> sigma (n, 0);
    sigma [s] = 1;
    std::vector <int> sigma_edge (nedges, 0);

    std::vector <std::vector <unsigned int> > prev_vert (n), prev_edge (n);

    while (m_heap->nItems() > 0) {
        unsigned int v = m_heap->deleteMin();

        if (w [v] > dist_threshold)
            continue;

        v_stack.push_back (v);

        edge = vertices [v].outHead;
        while (edge) {

            unsigned int et = edge->target;
            double wt = w [v] + edge->wt;

            // DGraph has no edge iterator, so edge_vec contains pairwise
            // elements of [from vertex, edge_id]
            std::vector <unsigned int> edge_vec;

            if (w [et] == 0.0) // first connection to et
            {
                edge_vec.resize (2);
                edge_vec [0] = v;
                edge_vec [1] = edge->edge_id;
                prev_edge [et] = edge_vec;

                sigma [et] = sigma [v];
                sigma_edge [edge->edge_id] = sigma [v];

                w [et] = wt;
                m_heap->insert (et, wt);
            }  else if (wt < w [et])
            {
                edge_vec.resize (2);
                edge_vec [0] = v;
                edge_vec [1] = edge->edge_id;
                prev_edge [et] = edge_vec;

                sigma [et] = sigma [v];
                sigma_edge [edge->edge_id] = sigma [v];

                w [et] = wt;
                m_heap->decreaseKey (et, wt);
            } else if (fabs (wt - w [et]) < epsilon)
            {
                edge_vec = prev_edge [et];
                edge_vec.resize (edge_vec.size () + 2);
                edge_vec [edge_vec.size () - 2] = v;
                edge_vec [edge_vec.size () - 1] = edge->edge_id;
                prev_edge [et] = edge_vec;

                sigma [et] += sigma [v];
                sigma_edge [edge->edge_id] += sigma [v];
            }

            edge = edge->nextOut;
        }
    } // end while nItems > 0

    // Then read from the stack and count centrality paths
    std::vector <double> delta (n, 0.0);
    while (!v_stack.empty ())
    {
        const unsigned int v = v_stack.back ();
        v_stack.pop_back ();
        std::vector <unsigned int> edge_vec = prev_edge [v];
        double tempd = (1.0 + delta [v]) / sigma [v];

        std::vector <unsigned int>::iterator it = edge_vec.begin ();
        // The dereferenced edge_vec iterator is simply a direct index
        while (it != edge_vec.end ())
        {
            delta [*it] += sigma [*it] * tempd;
            it = std::next (it);
            cent [*it] += sigma_edge [*it] * tempd * vert_wt;
            it = std::next (it);
        }
    }
}


//' rcpp_centrality_vertex - parallel function
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

    const unsigned int nedges = static_cast <unsigned int> (graph.nrow ());
    std::map <std::string, unsigned int> vert_map;
    std::vector <std::string> vert_map_id = vert_map_in ["vert"];
    std::vector <unsigned int> vert_map_n = vert_map_in ["id"];
    size_t nverts = run_sp::make_vert_map (vert_map_in, vert_map_id,
            vert_map_n, vert_map);

    Rcpp::CharacterVector v_nms = vert_map_in.attr ("names");
    std::vector <double> vert_wts (vert_map_in.nrow (), 1.0);
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
