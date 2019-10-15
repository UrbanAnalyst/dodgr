#include "run_sp.h"
#include "pathfinders.h"
#include "heaps/heap.h"

#include <algorithm> // std::fill

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
        g->addNewEdge (fromi, toi, dist [i], wt [i]);
    }
}
// # nocov end

void PF::PathFinder::Centrality (
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

    std::vector <std::vector <int> > prev_arr (n);

    while (m_heap->nItems() > 0) {
        unsigned int v = m_heap->deleteMin();

        v_stack.push_back (v);

        edge = vertices [v].outHead;
        while (edge) {

            unsigned int et = edge->target;
            double wt = w [v] + edge->dist;

            std::unordered_set <int> w_set;
            std::vector <int> w_vec;

            if (w [et] == 0.0) // first connection to et
            {
                w_vec.resize (1);
                w_vec [0] = v;
                prev_arr [et] = w_vec;

                sigma [et] = sigma [v];
                w [et] = wt;
                m_heap->insert (et, wt);
            }  else if (wt < w [et])
            {
                w_vec.resize (1);
                w_vec [0] = v;
                prev_arr [et] = w_vec;

                sigma [et] = sigma [v];
                w [et] = wt;
                m_heap->decreaseKey (et, wt);
            } else if (wt == w [et])
            {
                w_vec = prev_arr [et];
                w_vec.resize (w_vec.size () + 1);
                w_vec [w_vec.size () - 1] = v;
                prev_arr [et] = w_vec;

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
        std::vector <int> w_vec = prev_arr [v];
        double tempd = (1.0 + delta [v]) / sigma [v];
        for (auto ws: w_vec)
        {
            delta [ws] += sigma [ws] * tempd;
        }
        if (v != s)
            cent [v] += delta [v];
    }
}


//' rcpp_centrality
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_centrality (const Rcpp::DataFrame graph,
        const Rcpp::DataFrame vert_map_in,
        const std::string& heap_type)
{
    std::vector <std::string> from = graph ["from"];
    std::vector <std::string> to = graph ["to"];
    std::vector <double> dist = graph ["d"];
    std::vector <double> wt = graph ["d_weighted"];

    unsigned int nedges = static_cast <unsigned int> (graph.nrow ());
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
    std::vector <double> w (nverts), cent (nverts, 0.0);
    for (unsigned int i = 0; i < nverts; i++)
    {
        std::shared_ptr<PF::PathFinder> pathfinder =
            std::make_shared <PF::PathFinder> (nverts,
                *run_sp::getHeapImpl(heap_type), g);
        
        pathfinder->init (g); // specify the graph

        pathfinder->Centrality (w, cent, i);
    }

    Rcpp::NumericVector dout = Rcpp::wrap (cent);

    return (dout);
}

