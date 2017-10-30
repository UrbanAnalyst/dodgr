#include "dgraph.h"
#include "dijkstra.h"
#include "heaps/heap_lib.h"

#include "run_sp.h"

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

//' rcpp_get_sp_dists
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_get_sp_dists (Rcpp::DataFrame graph,
        Rcpp::DataFrame vert_map_in,
        std::vector <int> fromi,
        std::vector <int> toi,
        std::string heap_type)
{
    Rcpp::NumericVector id_vec;
    if (fromi [0] < 0) // use all vertices
    {
        id_vec = vert_map_in ["id"];
        fromi = Rcpp::as <std::vector <int> > (id_vec);
    }
    if (toi [0] < 0) // use all vertices
    {
        if (id_vec.size() == 0)
            id_vec = vert_map_in ["id"];
        toi = Rcpp::as <std::vector <int> > (id_vec);
    }
    size_t nfrom = fromi.size (), nto = toi.size ();

    std::vector <std::string> from = graph ["from"];
    std::vector <std::string> to = graph ["to"];
    std::vector <double> dist = graph ["d"];
    std::vector <double> wt = graph ["w"];

    unsigned int nedges = static_cast <unsigned int> (graph.nrow ());
    std::map <std::string, unsigned int> vert_map;
    std::vector <std::string> vert_map_id = vert_map_in ["vert"];
    std::vector <unsigned int> vert_map_n = vert_map_in ["id"];
    for (unsigned int i = 0;
            i < static_cast <unsigned int> (vert_map_in.nrow ()); ++i)
    {
        vert_map.emplace (vert_map_id [i], vert_map_n [i]);
    }
    unsigned int nverts = static_cast <unsigned int> (vert_map.size ());

    DGraph *g = new DGraph (nverts);
    inst_graph (g, nedges, vert_map, from, to, dist, wt);

    Dijkstra *dijkstra = std::nullptr_t ();

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

    double* w = new double [nverts];
    double* d = new double [nverts];
    int* prev = new int [nverts];

    dijkstra->init (g); // specify the graph

    // initialise dout matrix to NA
    Rcpp::NumericVector na_vec = Rcpp::NumericVector (nfrom * nto,
            Rcpp::NumericVector::get_na ());
    Rcpp::NumericMatrix dout (static_cast <int> (nfrom),
            static_cast <int> (nto), na_vec.begin ());
    for (unsigned int v = 0; v < nfrom; v++)
    {
        std::fill (w, w + nverts, INFINITE_DOUBLE);
        std::fill (d, d + nverts, INFINITE_DOUBLE);

        dijkstra->run (d, w, prev, static_cast <unsigned int> (fromi [v]));
        for (unsigned int vi = 0; vi < nto; vi++)
            if (w [toi [vi]] < INFINITE_DOUBLE)
                dout (v, vi) = d [toi [vi]];
    }

    delete [] d;
    delete [] w;
    delete [] prev;

    delete dijkstra;
    delete g;

    return (dout);
}


//' rcpp_get_paths
//'
//' @param graph The data.frame holding the graph edges
//' @param vert_map_in map from <std::string> vertex ID to (0-indexed) integer
//' index of vertices
//' @param fromi Index into vert_map_in of vertex numbers
//' @param toi Index into vert_map_in of vertex numbers
//'
//' @note The graph is constructed with 0-indexed vertex numbers contained in
//' code{vert_map_in}. Both \code{fromi} and \code{toi} already map directly
//' onto these. The graph has to be constructed by first constructing a
//' \code{std::map} object (\code{vertmap}) for \code{vert_map_in}, then
//' translating all \code{graph["from"/"to"]} values into these indices. This
//' construction is done in \code{inst_graph}.
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::List rcpp_get_paths (Rcpp::DataFrame graph,
        Rcpp::DataFrame vert_map_in,
        std::vector <int> fromi,
        std::vector <int> toi,
        std::string heap_type)
{
    Rcpp::NumericVector id_vec;
    if (fromi [0] < 0) // use all vertices
    {
        id_vec = vert_map_in ["id"];
        fromi = Rcpp::as <std::vector <int> > (id_vec);
    }
    if (toi [0] < 0) // use all vertices
    {
        if (id_vec.size() == 0)
            id_vec = vert_map_in ["id"];
        toi = Rcpp::as <std::vector <int> > (id_vec);
    }
    size_t nfrom = fromi.size (), nto = toi.size ();

    std::vector <std::string> from = graph ["from"];
    std::vector <std::string> to = graph ["to"];
    std::vector <double> dist = graph ["d"];
    std::vector <double> wt = graph ["w"];

    unsigned int nedges = static_cast <unsigned int> (graph.nrow ());
    std::map <std::string, unsigned int> vert_map;
    std::vector <std::string> vert_map_id = vert_map_in ["vert"];
    std::vector <unsigned int> vert_map_n = vert_map_in ["id"];
    for (unsigned int i = 0;
            i < static_cast <unsigned int> (vert_map_in.nrow ()); ++i)
    {
        vert_map.emplace (vert_map_id [i], vert_map_n [i]);
    }
    unsigned int nverts = static_cast <unsigned int> (vert_map.size ());

    DGraph *g = new DGraph (nverts);
    inst_graph (g, nedges, vert_map, from, to, dist, wt);

    Dijkstra *dijkstra = std::nullptr_t ();

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

    dijkstra->init (g); // specify the graph

    Rcpp::List res (nfrom);
    double* w = new double [nverts];
    double* d = new double [nverts];
    int* prev = new int [nverts];

    for (unsigned int v = 0; v < nfrom; v++)
    {
        std::fill (w, w + nverts, INFINITE_DOUBLE);
        std::fill (d, d + nverts, INFINITE_DOUBLE);

        dijkstra->run (d, w, prev, static_cast <unsigned int> (fromi [v]));

        Rcpp::List res1 (nto);
        for (unsigned int vi = 0; vi < nto; vi++)
        {
            std::vector <unsigned int> onePath;
            if (w [toi [vi]] < INFINITE_DOUBLE)
            {
                int target = toi [vi];
                while (target < INFINITE_INT)
                {
                    // Note that targets are all C++ 0-indexed and are converted
                    // directly here to R-style 1-indexes.
                    onePath.push_back (static_cast <unsigned int> (target + 1));
                    target = prev [target];
                    if (target < 0)
                        break;
                }
            }
            std::reverse (onePath.begin (), onePath.end ());
            res1 [vi] = onePath;
        }
        res [v] = res1;
    }


    delete [] d;
    delete [] w;
    delete [] prev;

    delete dijkstra;
    delete g;

    return (res);
}

//' rcpp_aggregate_flows
//'
//' @param graph The data.frame holding the graph edges
//' @param vert_map_in map from <std::string> vertex ID to (0-indexed) integer
//' index of vertices
//' @param fromi Index into vert_map_in of vertex numbers
//' @param toi Index into vert_map_in of vertex numbers
//'
//' @note The flow data to be used for aggregation is a matrix mapping flows
//' betwen each pair of from and to points.
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_aggregate_flows (Rcpp::DataFrame graph,
        Rcpp::DataFrame vert_map_in,
        std::vector <int> fromi,
        std::vector <int> toi,
        Rcpp::NumericMatrix flows,
        std::string heap_type)
{
    Rcpp::NumericVector id_vec;
    if (fromi [0] < 0) // use all vertices
    {
        id_vec = vert_map_in ["id"];
        fromi = Rcpp::as <std::vector <int> > (id_vec);
    }
    if (toi [0] < 0) // use all vertices
    {
        if (id_vec.size() == 0)
            id_vec = vert_map_in ["id"];
        toi = Rcpp::as <std::vector <int> > (id_vec);
    }
    size_t nfrom = fromi.size (), nto = toi.size ();

    std::vector <std::string> from = graph ["from"];
    std::vector <std::string> to = graph ["to"];
    std::vector <double> dist = graph ["d"];
    std::vector <double> wt = graph ["w"];

    unsigned int nedges = static_cast <unsigned int> (graph.nrow ());
    std::vector <std::string> vert_name = vert_map_in ["vert"];
    std::vector <unsigned int> vert_indx = vert_map_in ["id"];
    // Make map from vertex name to integer index
    std::map <std::string, unsigned int> vert_map_i;
    for (unsigned int i = 0;
            i < static_cast <unsigned int> (vert_map_in.nrow ()); ++i)
    {
        vert_map_i.emplace (vert_name [i], vert_indx [i]);
    }
    unsigned int nverts = static_cast <unsigned int> (vert_map_i.size ());

    /* Flows from the dijkstra output are reallocated based on matching vertex
     * pairs to edge indices. Note, however, that contracted graphs frequently
     * have duplicate vertex pairs with different distances. The following
     * therefore uses two maps, one to hold the ultimate index from vertex
     * pairs, and the other to hold minimal distances.
     */
    std::unordered_map <std::string, unsigned int> verts_to_edge_map;
    std::unordered_map <std::string, double> verts_to_dist_map;
    for (unsigned int i = 0; i < from.size (); i++)
    {
        std::string two_verts = "f" + from [i] + "t" + to [i];
        verts_to_edge_map.emplace (two_verts, i);
        if (verts_to_dist_map.find (two_verts) == verts_to_dist_map.end ())
            verts_to_dist_map.emplace (two_verts, wt [i]);
        else if (wt [i] < verts_to_dist_map.at (two_verts))
        {
            verts_to_dist_map [two_verts] = wt [i];
            verts_to_edge_map [two_verts] = i;
        }
    }

    DGraph *g = new DGraph (nverts);
    inst_graph (g, nedges, vert_map_i, from, to, dist, wt);

    Dijkstra *dijkstra = std::nullptr_t ();

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

    dijkstra->init (g); // specify the graph

    Rcpp::List res (nfrom);
    double* w = new double [nverts];
    double* d = new double [nverts];
    int* prev = new int [nverts];

    unsigned int junk = 0; // TODO: Delete that!

    Rcpp::NumericVector aggregate_flows (from.size ()); // 0-filled by default
    for (unsigned int v = 0; v < nfrom; v++)
    {
        std::fill (w, w + nverts, INFINITE_DOUBLE);
        std::fill (d, d + nverts, INFINITE_DOUBLE);

        dijkstra->run (d, w, prev, static_cast <unsigned int> (fromi [v]));

        Rcpp::List res1 (nto);
        for (unsigned int vi = 0; vi < nto; vi++)
        {
            if (vi != v) // Exclude self-flows
            {
                std::vector <unsigned int> onePath;
                double flow_ij = flows (v, vi);
                if (flow_ij > 0)
                    junk++;
                if (w [toi [vi]] < INFINITE_DOUBLE)
                {
                    // target values are int indices into vert_map_in, which means
                    // corresponding vertex IDs can be taken directly from
                    // vert_name
                    int target = toi [vi];
                    while (target < INFINITE_INT)
                    {
                        if (prev [target] >= 0 && prev [target] < INFINITE_INT)
                        {
                            std::string v2 = "f" +
                                vert_name [static_cast <size_t> (prev [target])] +
                                "t" + vert_name [static_cast <size_t> (target)];
                            aggregate_flows [verts_to_edge_map.at (v2)] += flow_ij;
                        }

                        target = prev [target];
                        // Only allocate that flow from origin vertex v to all
                        // previous vertices up until the target vi
                        if (target < 0 || target == fromi [v])
                        {
                            break;
                        }
                    } // end while target
                } // end if w < INF
            } // end if vi != v
        } // end for vi over nto
    } // end for v over nfrom

    delete [] d;
    delete [] w;
    delete [] prev;

    delete dijkstra;
    delete g;

    return (aggregate_flows);
}
