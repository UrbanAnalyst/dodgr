
#include "dgraph.h"
#include "dijkstra.h"
#include "heap_lib.h"

#include <Rcpp.h>

const float INFINITE_FLOAT =  std::numeric_limits<float>::max ();
const float INFINITE_INT =  std::numeric_limits<int>::max ();

unsigned int inst_vert_map (unsigned int nedges,
        std::map <std::string, unsigned int> &vert_map,
        std::vector <std::string> &from,
        std::vector <std::string> &to)
{
    unsigned int nverts = 0;
    for (unsigned int i = 0; i < nedges; i++)
    {
        if (vert_map.find (from [i]) == vert_map.end ())
            vert_map.emplace (from [i], nverts++);
        if (vert_map.find (to [i]) == vert_map.end ())
            vert_map.emplace (to [i], nverts++);
    }

    return nverts;
}

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
template void inst_graph <int> (DGraph *g, unsigned int nedges,
        std::map <std::string, unsigned int> &vert_map,
        std::vector <std::string> &from,
        std::vector <std::string> &to,
        std::vector <int> &dist,
        std::vector <int> &wt);

void inst_graph_int (DGraphInt *g, unsigned int nedges,
        std::map <std::string, unsigned int> &vert_map,
        std::vector <std::string> &from,
        std::vector <std::string> &to,
        std::vector <int> &dist,
        std::vector <int> &wt)
{
    for (unsigned int i = 0; i < nedges; ++i)
    {
        unsigned int fromi = vert_map [from [i]];
        unsigned int toi = vert_map [to [i]];
        g->addNewEdge (fromi, toi, dist [i], wt [i]);
    }
}

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

//' rcpp_get_sp
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_get_sp (Rcpp::DataFrame graph, std::string heap_type)
{
    std::vector <std::string> from = graph ["from_id"];
    std::vector <std::string> to = graph ["to_id"];
    std::vector <float> dist = graph ["d"];
    std::vector <float> wt = graph ["d_weighted"];

    unsigned int nedges = graph.nrow ();
    std::map <std::string, unsigned int> vert_map;
    unsigned int nverts = inst_vert_map (nedges, vert_map, from, to);

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

    float* w = new float [nverts];
    float* d = new float [nverts];
    for(unsigned int v = 0; v < nverts; v++)
    {
        w [v] = INFINITE_FLOAT;
        d [v] = INFINITE_FLOAT;
    }

    Rcpp::NumericMatrix dout (nverts, nverts);
    for(unsigned int v = 0; v < nverts; v++)
    {
        dijkstra->init (g); // specify the graph
        dijkstra->run (d, w, v);
        for(unsigned int vi = 0; vi < nverts; vi++)
            dout (v, vi) = d [vi];
    }

    delete dijkstra;
    delete g;

    delete [] d;
    delete [] w;

    return (dout);
}

//' rcpp_get_sp_radix
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_get_sp_radix (Rcpp::DataFrame graph)
{
    std::vector <std::string> from = graph ["from_id"];
    std::vector <std::string> to = graph ["to_id"];
    std::vector <int> dist = graph ["d"];
    std::vector <int> wt = graph ["d_weighted"];

    unsigned int nedges = graph.nrow ();
    std::map <std::string, unsigned int> vert_map;
    unsigned int nverts = inst_vert_map (nedges, vert_map, from, to);

    DGraphInt *g = new DGraphInt (nverts);
    inst_graph_int (g, nedges, vert_map, from, to, dist, wt);

    IntHeapD <RadixHeap> heapD;
    DijkstraInt *dijkstra = new DijkstraInt (nverts, &heapD);
    
    int* w = new int [nverts];
    int* d = new int [nverts];
    for(unsigned int v = 0; v < nverts; v++)
    {
        w [v] = INFINITE_INT;
        d [v] = INFINITE_INT;
    }

    Rcpp::NumericMatrix dout (nverts, nverts);
    for(unsigned int v = 0; v < nverts; v++)
    {
        dijkstra->init (g); // specify the graph
        dijkstra->run (d, w, v);
        for(unsigned int vi = 0; vi < nverts; vi++)
            dout (v, vi) = d [vi];
    }

    delete dijkstra;
    delete g;

    delete [] d;
    delete [] w;

    return (dout);
}
