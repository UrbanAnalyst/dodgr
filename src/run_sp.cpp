
#include "run_sp.h"

#include "dgraph.h"
#include "heaps/heap_lib.h"

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


std::shared_ptr <HeapDesc> getHeapImpl(const std::string& heap_type)
{
  if (heap_type == "FHeap")
    return std::make_shared <HeapD<FHeap> >();
  else if (heap_type == "BHeap")
    return std::make_shared <HeapD<BHeap> >();
  else if (heap_type == "Heap23")
    return std::make_shared <HeapD<Heap23> >();
  else if (heap_type == "TriHeap")
    return std::make_shared <HeapD<TriHeap> >();
  else if (heap_type == "TriHeapExt")
    return std::make_shared <HeapD<TriHeapExt> >();
  else if (heap_type == "Radix")
    return std::make_shared <HeapD<RadixHeap> >();
  else 
    throw std::runtime_error("invalid heap type: " + heap_type);
}


struct oneDijkstra : public RcppParallel::Worker
{
    RcppParallel::RVector <int> dp_fromi;
    Rcpp::IntegerVector toi;
    size_t nverts;

    std::shared_ptr <DGraph> g;
    std::string heap_type;

    RcppParallel::RMatrix <double> dout;

    // constructor
    oneDijkstra (
            const Rcpp::IntegerVector fromi,
            const Rcpp::IntegerVector toi,
            const size_t nverts,
            const std::shared_ptr <DGraph> g,
            const std::string & heap_type,
            Rcpp::NumericMatrix dout) :
        dp_fromi (fromi), toi (toi), nverts (nverts),
        g (g), heap_type (heap_type),
        dout (dout)
    {
    }

    // Parallel function operator
    void operator() (std::size_t begin, std::size_t end)
    {
        std::shared_ptr<Dijkstra> dijkstra =
            std::make_shared <Dijkstra> (nverts,
                    *getHeapImpl (heap_type), g);
        std::vector <double> w (nverts);
        std::vector <double> d (nverts);
        std::vector <int> prev (nverts);

        for (std::size_t i = begin; i < end; i++)
        {
            // These have to be reserved within the parallel operator function!
            std::fill (w.begin (), w.end (), INFINITE_DOUBLE);
            std::fill (d.begin (), d.end (), INFINITE_DOUBLE);

            dijkstra->run (d, w, prev,
                    static_cast <unsigned int> (dp_fromi [i]));
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

//' rcpp_get_sp_dists_par
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_get_sp_dists_par (Rcpp::DataFrame graph,
        Rcpp::DataFrame vert_map_in,
        Rcpp::IntegerVector fromi,
        Rcpp::IntegerVector toi,
        const std::string& heap_type)
{
    Rcpp::NumericVector id_vec;
    if (fromi [0] < 0) // use all vertices
    {
        id_vec = vert_map_in ["id"];
        fromi = id_vec;
    }
    if (toi [0] < 0) // use all vertices
    {
        if (id_vec.size () == 0)
            id_vec = vert_map_in ["id"];
        toi = id_vec;
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
    size_t nverts = static_cast <size_t> (vert_map.size ());

    std::shared_ptr <DGraph> g = std::make_shared <DGraph> (nverts);
    inst_graph (g, nedges, vert_map, from, to, dist, wt);

    Rcpp::NumericVector na_vec = Rcpp::NumericVector (nfrom * nto,
            Rcpp::NumericVector::get_na ());
    Rcpp::NumericMatrix dout (static_cast <int> (nfrom),
            static_cast <int> (nto), na_vec.begin ());

    // Create parallel worker
    oneDijkstra one_dijkstra (fromi, toi, nverts, g, heap_type, dout);

    RcppParallel::parallelFor (0, fromi.length (), one_dijkstra);
    
    return (dout);
}

//' rcpp_get_sp_dists
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_get_sp_dists (Rcpp::DataFrame graph,
        Rcpp::DataFrame vert_map_in,
        std::vector <int> fromi,
        std::vector <int> toi,
        const std::string& heap_type)
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

    std::shared_ptr<DGraph> g = std::make_shared<DGraph>(nverts);
    inst_graph (g, nedges, vert_map, from, to, dist, wt);

    std::shared_ptr<Dijkstra> dijkstra = std::make_shared<Dijkstra>(nverts, *getHeapImpl(heap_type), g);

    std::vector<double> w(nverts);
    std::vector<double> d(nverts);
    std::vector<int> prev(nverts);

    //dijkstra->init (g); // specify the graph

    // initialise dout matrix to NA
    Rcpp::NumericVector na_vec = Rcpp::NumericVector (nfrom * nto,
            Rcpp::NumericVector::get_na ());
    Rcpp::NumericMatrix dout (static_cast <int> (nfrom),
            static_cast <int> (nto), na_vec.begin ());

    for (unsigned int v = 0; v < nfrom; v++)
    {
        std::fill (w.begin(), w.end(), INFINITE_DOUBLE);
        std::fill (d.begin(), d.end(), INFINITE_DOUBLE);

        dijkstra->run (d, w, prev, static_cast <unsigned int> (fromi [v]));
        for (unsigned int vi = 0; vi < nto; vi++)
            if (w [toi [vi]] < INFINITE_DOUBLE)
            {
                dout (v, vi) = d [toi [vi]];
            }
    }
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
        const std::string& heap_type)
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

    std::shared_ptr<DGraph> g = std::make_shared<DGraph>(nverts);
    inst_graph (g, nedges, vert_map, from, to, dist, wt);

    std::shared_ptr<Dijkstra> dijkstra = std::make_shared<Dijkstra>(nverts, *getHeapImpl(heap_type), g);
    
    Rcpp::List res (nfrom);
    std::vector<double> w(nverts);
    std::vector<double> d(nverts);
    std::vector<int> prev(nverts);

    for (unsigned int v = 0; v < nfrom; v++)
    {
        std::fill (w.begin(), w.end(), INFINITE_DOUBLE);
        std::fill (d.begin(), d.end(), INFINITE_DOUBLE);

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

    std::shared_ptr<DGraph> g = std::make_shared<DGraph> (nverts);
    inst_graph (g, nedges, vert_map_i, from, to, dist, wt);

    std::shared_ptr<Dijkstra> dijkstra = std::make_shared<Dijkstra>(nverts, *getHeapImpl(heap_type), g);

    Rcpp::List res (nfrom);
    std::vector<double> w(nverts);
    std::vector<double> d(nverts);
    std::vector<int> prev(nverts);

    Rcpp::NumericVector aggregate_flows (from.size ()); // 0-filled by default
    for (unsigned int v = 0; v < nfrom; v++)
    {
        std::fill (w.begin(), w.end(), INFINITE_DOUBLE);
        std::fill (d.begin(), d.end(), INFINITE_DOUBLE);

        dijkstra->run (d, w, prev, static_cast <unsigned int> (fromi [v]));

        Rcpp::List res1 (nto);
        for (unsigned int vi = 0; vi < nto; vi++)
        {
            if (fromi [v] != toi [vi]) // Exclude self-flows
            {
                std::vector <unsigned int> onePath;
                double flow_ij = flows (v, vi);
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

    return (aggregate_flows);
}
