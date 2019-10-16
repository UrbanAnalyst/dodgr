#include "run_sp.h"
#include "pathfinders.h"
#include "heaps/heap.h"

#include <algorithm> // std::fill

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

void PF::PathFinder::Centrality_vertex (
        std::vector <double>& w,
        std::vector <double>& cent,
        const unsigned int s)
{
    const DGraphEdge *edge;

    const unsigned int n = m_graph->nVertices();
    const std::vector <DGraphVertex>& vertices = m_graph->vertices();

    std::deque <unsigned int> v_stack;

    std::fill (w.begin (), w.end (), 0.0);
    w [s] = 1.0;

    m_heap->insert (s, -1.0);

    std::vector <int> sigma (n, 0);
    sigma [s] = 1;

    std::vector <std::vector <unsigned int> > prev_vert (n);

    while (m_heap->nItems() > 0) {
        unsigned int v = m_heap->deleteMin();

        v_stack.push_back (v);

        edge = vertices [v].outHead;
        while (edge) {

            unsigned int et = edge->target;
            double wt = w [v] + edge->dist;

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
            cent [v] += delta [v];
    }
}


void PF::PathFinder::Centrality_edge (
        std::vector <double>& w,
        std::vector <double>& cent,
        const unsigned int s,
        const unsigned int nedges)
{
    const DGraphEdge *edge;

    const unsigned int n = m_graph->nVertices();
    const std::vector <DGraphVertex>& vertices = m_graph->vertices();

    std::deque <unsigned int> v_stack;

    std::fill (w.begin (), w.end (), 0.0);
    w [s] = 1.0;

    m_heap->insert (s, -1.0);

    std::vector <int> sigma (n, 0);
    sigma [s] = 1;
    std::vector <int> sigma_edge (nedges, 0);

    std::vector <std::vector <unsigned int> > prev_vert (n), prev_edge (n);

    while (m_heap->nItems() > 0) {
        unsigned int v = m_heap->deleteMin();

        v_stack.push_back (v);

        edge = vertices [v].outHead;
        while (edge) {

            unsigned int et = edge->target;
            double wt = w [v] + edge->dist;

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
            cent [*it] += sigma_edge [*it] * tempd;
            it = std::next (it);
        }
    }
}


//' rcpp_centrality
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_centrality (const Rcpp::DataFrame graph,
        const Rcpp::DataFrame vert_map_in,
        const std::string& heap_type,
        bool edges)
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

    std::shared_ptr <DGraph> g = std::make_shared <DGraph> (nverts);
    inst_graph (g, nedges, vert_map, from, to, dist, wt);

    Rcpp::NumericVector na_vec = Rcpp::NumericVector (nverts,
            Rcpp::NumericVector::get_na ());
    //Rcpp::NumericVector dout (static_cast <int> (nverts), na_vec.begin ());

    // Create parallel worker
    /*
    OneIso one_iso (fromi, nverts, g, dlim, heap_type, dout);

    RcppParallel::parallelFor (0, static_cast <size_t> (fromi.length ()),
            one_iso);
    */
    std::vector <double> w (nverts), cent;
    if (edges)
        cent.resize (nedges);
    else
        cent.resize (nverts);
    std::fill (cent.begin (), cent.end (), 0.0);

    for (unsigned int i = 0; i < nverts; i++)
    {
        std::shared_ptr<PF::PathFinder> pathfinder =
            std::make_shared <PF::PathFinder> (nverts,
                *run_sp::getHeapImpl(heap_type), g);
        
        pathfinder->init (g); // specify the graph

        if (edges)
            pathfinder->Centrality_edge (w, cent, i, nedges);
        else
            pathfinder->Centrality_vertex (w, cent, i);
    }

    Rcpp::NumericVector dout = Rcpp::wrap (cent);

    return (dout);
}

