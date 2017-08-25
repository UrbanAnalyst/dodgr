#include "dgraph.h"
#include "dijkstra.h"
#include "heaps/heap_lib.h"

#include <algorithm> // std::fill

#include <Rcpp.h>

const float INFINITE_FLOAT =  std::numeric_limits<float>::max ();
const float INFINITE_INT =  std::numeric_limits<int>::max ();

template <typename T>
void inst_graph (DGraph *g, unsigned int nedges,
        std::map <std::string, unsigned int> &vert_map,
        std::vector <std::string> &from,
        std::vector <std::string> &to,
        std::vector <T> &dist,
        std::vector <T> &wt)
{
    for (unsigned int i = 0; i < nedges; ++i)
    {
        unsigned int fromi = vert_map [from [i]];
        unsigned int toi = vert_map [to [i]];
        g->addNewEdge (fromi, toi, dist [i], wt [i]);
    }
}
template void inst_graph <float> (DGraph *g, unsigned int nedges,
        std::map <std::string, unsigned int> &vert_map,
        std::vector <std::string> &from,
        std::vector <std::string> &to,
        std::vector <float> &dist,
        std::vector <float> &wt);

Dijkstra * dijkstra_bheap (unsigned int nverts)
{
    HeapD<BHeap> heapD;
    Dijkstra *dijkstra = new Dijkstra (nverts, &heapD);
    return dijkstra;
}

Dijkstra * dijkstra_fheap (unsigned int nverts)
{
    HeapD<FHeap> heapD;
    Dijkstra *dijkstra = new Dijkstra (nverts, &heapD);
    return dijkstra;
}

Dijkstra * dijkstra_heap23 (unsigned int nverts)
{
    HeapD<Heap23> heapD;
    Dijkstra *dijkstra = new Dijkstra (nverts, &heapD);
    return dijkstra;
}

Dijkstra * dijkstra_triheap (unsigned int nverts)
{
    HeapD<TriHeap> heapD;
    Dijkstra *dijkstra = new Dijkstra (nverts, &heapD);
    return dijkstra;
}

Dijkstra * dijkstra_triheapext (unsigned int nverts)
{
    HeapD<TriHeapExt> heapD;
    Dijkstra *dijkstra = new Dijkstra (nverts, &heapD);
    return dijkstra;
}

Dijkstra * dijkstra_radix (unsigned int nverts)
{
    HeapD<RadixHeap> heapD;
    Dijkstra *dijkstra = new Dijkstra (nverts, &heapD);
    return dijkstra;
}

//' rcpp_get_sp
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_get_sp (Rcpp::DataFrame graph,
        Rcpp::DataFrame vert_map_in,
        std::vector <int> fromi,
        std::vector <int> toi,
        std::string heap_type)
{
    if (fromi [0] < 0) // use all vertices
    {
        Rcpp::NumericVector id_vec = vert_map_in ["id"];
        fromi = Rcpp::as <std::vector <int>> (id_vec);
    }
    if (toi [0] < 0) // use all vertices
    {
        Rcpp::NumericVector id_vec = vert_map_in ["id"];
        toi = Rcpp::as <std::vector <int>> (id_vec);
    }
    unsigned int nfrom = fromi.size (), nto = toi.size ();

    std::vector <std::string> from = graph ["from"];
    std::vector <std::string> to = graph ["to"];
    std::vector <float> dist = graph ["d"];
    std::vector <float> wt = graph ["w"];

    unsigned int nedges = graph.nrow ();
    std::map <std::string, unsigned int> vert_map;
    std::vector <std::string> vert_map_id = vert_map_in ["vert"];
    std::vector <unsigned int> vert_map_n = vert_map_in ["id"];
    for (int i = 0; i < vert_map_in.nrow (); ++i)
    {
        vert_map.emplace (vert_map_id [i], vert_map_n [i]);
    }
    unsigned int nverts = vert_map.size ();

    DGraph *g = new DGraph (nverts);
    inst_graph (g, nedges, vert_map, from, to, dist, wt);

    Dijkstra *dijkstra;

    if (heap_type == "FHeap")
        dijkstra = dijkstra_fheap (nverts);
    else if (heap_type == "BHeap")
        dijkstra = dijkstra_bheap (nverts);
    else if (heap_type == "Heap23")
        dijkstra = dijkstra_heap23 (nverts);
    else if (heap_type == "TriHeap")
        dijkstra = dijkstra_triheap (nverts);
    else if (heap_type == "TriHeapExt")
        dijkstra = dijkstra_triheapext (nverts);
    else if (heap_type == "Radix")
        dijkstra = dijkstra_radix (nverts);

    float* w = new float [nverts];
    float* d = new float [nverts];
    int* prev = new int [nverts];

    dijkstra->init (g); // specify the graph

    // initialise dout matrix to NA
    Rcpp::NumericVector na_vec = Rcpp::NumericVector (nfrom * nto,
            Rcpp::NumericVector::get_na ());
    Rcpp::NumericMatrix dout (nfrom, nto, na_vec.begin ());
    for (unsigned int v = 0; v < nfrom; v++)
    {
        std::fill (w, w + nverts, INFINITE_FLOAT);

        dijkstra->run (d, w, prev, fromi [v]);
        for (unsigned int vi = 0; vi < nto; vi++)
            if (w [toi [vi]] < INFINITE_FLOAT)
                dout (v, vi) = d [toi [vi]];
    }

    // Tracing a path:
    /*
       std::fill (w, w + nverts, INFINITE_FLOAT);
       dijkstra->run (d, w, prev, 2);
       int target = 50;
       while (target > 0)
       {
       float di = d [target];
       Rcpp::Rcout << "[" << target << " -> ";
       target = prev [target];
       Rcpp::Rcout << target << "]: " << di << std::endl;
       }
       */

    delete [] d;
    delete [] w;
    delete [] prev;

    delete dijkstra;
    delete g;

    return (dout);
}
