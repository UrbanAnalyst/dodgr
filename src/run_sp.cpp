
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


struct OneDist : public RcppParallel::Worker
{
    RcppParallel::RVector <int> dp_fromi;
    Rcpp::IntegerVector toi;
    size_t nverts;

    std::shared_ptr <DGraph> g;
    std::string heap_type;

    RcppParallel::RMatrix <double> dout;

    // constructor
    OneDist (
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


size_t make_vert_map (const Rcpp::DataFrame &vert_map_in,
        const std::vector <std::string> &vert_map_id,
        const std::vector <unsigned int> &vert_map_n,
        std::map <std::string, unsigned int> &vert_map)
{
    for (unsigned int i = 0;
            i < static_cast <unsigned int> (vert_map_in.nrow ()); ++i)
    {
        vert_map.emplace (vert_map_id [i], vert_map_n [i]);
    }
    size_t nverts = static_cast <size_t> (vert_map.size ());
    return (nverts);
}

size_t get_fromi_toi (const Rcpp::DataFrame &vert_map_in,
        Rcpp::IntegerVector &fromi, Rcpp::IntegerVector &toi,
        Rcpp::NumericVector &id_vec)
{
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
    size_t nfrom = fromi.size ();
    return nfrom;
}

size_t get_fromi (const Rcpp::DataFrame &vert_map_in,
        Rcpp::IntegerVector &fromi, Rcpp::NumericVector &id_vec)
{
    if (fromi [0] < 0) // use all vertices
    {
        id_vec = vert_map_in ["id"];
        fromi = id_vec;
    }
    size_t nfrom = fromi.size ();
    return nfrom;
}

// Flows from the dijkstra output are reallocated based on matching vertex
// pairs to edge indices. Note, however, that contracted graphs frequently
// have duplicate vertex pairs with different distances. The following
// therefore uses two maps, one to hold the ultimate index from vertex
// pairs, and the other to hold minimal distances. This is used in flow routines
// only.
void make_vert_to_edge_maps (const std::vector <std::string> &from,
        const std::vector <std::string> &to, const std::vector <double> &wt,
        std::unordered_map <std::string, unsigned int> &verts_to_edge_map,
        std::unordered_map <std::string, double> &verts_to_dist_map)
{
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
}

//' rcpp_get_sp_dists_par
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_get_sp_dists_par (const Rcpp::DataFrame graph,
        const Rcpp::DataFrame vert_map_in,
        Rcpp::IntegerVector fromi,
        Rcpp::IntegerVector toi,
        const std::string& heap_type)
{
    Rcpp::NumericVector id_vec;
    size_t nfrom = get_fromi_toi (vert_map_in, fromi, toi, id_vec);
    size_t nto = toi.size ();

    std::vector <std::string> from = graph ["from"];
    std::vector <std::string> to = graph ["to"];
    std::vector <double> dist = graph ["d"];
    std::vector <double> wt = graph ["w"];

    unsigned int nedges = static_cast <unsigned int> (graph.nrow ());
    std::map <std::string, unsigned int> vert_map;
    std::vector <std::string> vert_map_id = vert_map_in ["vert"];
    std::vector <unsigned int> vert_map_n = vert_map_in ["id"];
    size_t nverts = make_vert_map (vert_map_in, vert_map_id,
            vert_map_n, vert_map);

    std::shared_ptr <DGraph> g = std::make_shared <DGraph> (nverts);
    inst_graph (g, nedges, vert_map, from, to, dist, wt);

    Rcpp::NumericVector na_vec = Rcpp::NumericVector (nfrom * nto,
            Rcpp::NumericVector::get_na ());
    Rcpp::NumericMatrix dout (static_cast <int> (nfrom),
            static_cast <int> (nto), na_vec.begin ());

    // Create parallel worker
    OneDist one_dist (fromi, toi, nverts, g, heap_type, dout);

    RcppParallel::parallelFor (0, fromi.length (), one_dist);
    
    return (dout);
}

//' rcpp_get_sp_dists
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_get_sp_dists (const Rcpp::DataFrame graph,
        const Rcpp::DataFrame vert_map_in,
        Rcpp::IntegerVector fromi,
        Rcpp::IntegerVector toi,
        const std::string& heap_type)
{
    Rcpp::NumericVector id_vec;
    size_t nfrom = get_fromi_toi (vert_map_in, fromi, toi, id_vec);
    size_t nto = toi.size ();

    std::vector <std::string> from = graph ["from"];
    std::vector <std::string> to = graph ["to"];
    std::vector <double> dist = graph ["d"];
    std::vector <double> wt = graph ["w"];

    unsigned int nedges = static_cast <unsigned int> (graph.nrow ());
    std::map <std::string, unsigned int> vert_map;
    std::vector <std::string> vert_map_id = vert_map_in ["vert"];
    std::vector <unsigned int> vert_map_n = vert_map_in ["id"];
    size_t nverts = make_vert_map (vert_map_in, vert_map_id,
            vert_map_n, vert_map);

    std::shared_ptr<DGraph> g = std::make_shared<DGraph>(nverts);
    inst_graph (g, nedges, vert_map, from, to, dist, wt);

    std::shared_ptr <Dijkstra> dijkstra = std::make_shared <Dijkstra> (nverts,
            *getHeapImpl(heap_type), g);

    std::vector<double> w (nverts);
    std::vector<double> d (nverts);
    std::vector<int> prev (nverts);

    dijkstra->init (g); // specify the graph

    // initialise dout matrix to NA
    Rcpp::NumericVector na_vec = Rcpp::NumericVector (nfrom * nto,
            Rcpp::NumericVector::get_na ());
    Rcpp::NumericMatrix dout (static_cast <int> (nfrom),
            static_cast <int> (nto), na_vec.begin ());

    for (unsigned int v = 0; v < nfrom; v++)
    {
        Rcpp::checkUserInterrupt ();
        std::fill (w.begin(), w.end(), INFINITE_DOUBLE);
        std::fill (d.begin(), d.end(), INFINITE_DOUBLE);

        dijkstra->run (d, w, prev, static_cast <unsigned int> (fromi [v]));
        for (unsigned int vi = 0; vi < nto; vi++)
        {
            if (toi [vi] < INFINITE_INT)
            {
                dout (v, vi) = d [toi [vi]];
            }
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
Rcpp::List rcpp_get_paths (const Rcpp::DataFrame graph,
        const Rcpp::DataFrame vert_map_in,
        Rcpp::IntegerVector fromi,
        Rcpp::IntegerVector toi,
        const std::string& heap_type)
{
    Rcpp::NumericVector id_vec;
    size_t nfrom = get_fromi_toi (vert_map_in, fromi, toi, id_vec);
    size_t nto = toi.size ();

    std::vector <std::string> from = graph ["from"];
    std::vector <std::string> to = graph ["to"];
    std::vector <double> dist = graph ["d"];
    std::vector <double> wt = graph ["w"];

    unsigned int nedges = static_cast <unsigned int> (graph.nrow ());
    std::map <std::string, unsigned int> vert_map;
    std::vector <std::string> vert_map_id = vert_map_in ["vert"];
    std::vector <unsigned int> vert_map_n = vert_map_in ["id"];
    size_t nverts = make_vert_map (vert_map_in, vert_map_id,
            vert_map_n, vert_map);

    std::shared_ptr<DGraph> g = std::make_shared<DGraph>(nverts);
    inst_graph (g, nedges, vert_map, from, to, dist, wt);

    std::shared_ptr<Dijkstra> dijkstra = std::make_shared <Dijkstra> (nverts,
            *getHeapImpl(heap_type), g);
    
    Rcpp::List res (nfrom);
    std::vector<double> w(nverts);
    std::vector<double> d(nverts);
    std::vector<int> prev(nverts);

    for (unsigned int v = 0; v < nfrom; v++)
    {
        Rcpp::checkUserInterrupt ();
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

//' rcpp_flows_aggregate
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
Rcpp::NumericVector rcpp_flows_aggregate (const Rcpp::DataFrame graph,
        const Rcpp::DataFrame vert_map_in,
        Rcpp::IntegerVector fromi,
        Rcpp::IntegerVector toi,
        Rcpp::NumericMatrix flows,
        std::string heap_type)
{
    Rcpp::NumericVector id_vec;
    size_t nfrom = get_fromi_toi (vert_map_in, fromi, toi, id_vec);
    size_t nto = toi.size ();

    std::vector <std::string> from = graph ["from"];
    std::vector <std::string> to = graph ["to"];
    std::vector <double> dist = graph ["d"];
    std::vector <double> wt = graph ["w"];

    unsigned int nedges = static_cast <unsigned int> (graph.nrow ());
    std::vector <std::string> vert_name = vert_map_in ["vert"];
    std::vector <unsigned int> vert_indx = vert_map_in ["id"];
    // Make map from vertex name to integer index
    std::map <std::string, unsigned int> vert_map_i;
    size_t nverts = make_vert_map (vert_map_in, vert_name,
            vert_indx, vert_map_i);

    std::unordered_map <std::string, unsigned int> verts_to_edge_map;
    std::unordered_map <std::string, double> verts_to_dist_map;
    make_vert_to_edge_maps (from, to, wt, verts_to_edge_map, verts_to_dist_map);

    std::shared_ptr<DGraph> g = std::make_shared<DGraph> (nverts);
    inst_graph (g, nedges, vert_map_i, from, to, dist, wt);

    std::shared_ptr<Dijkstra> dijkstra = std::make_shared <Dijkstra> (nverts,
            *getHeapImpl(heap_type), g);

    Rcpp::List res (nfrom);
    std::vector<double> w(nverts);
    std::vector<double> d(nverts);
    std::vector<int> prev(nverts);

    Rcpp::NumericVector aggregate_flows (from.size ()); // 0-filled by default
    for (unsigned int v = 0; v < nfrom; v++)
    {
        Rcpp::checkUserInterrupt ();
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

//' rcpp_flows_disperse
//'
//' Modified version of \code{rcpp_aggregate_flows} that aggregates flows to all
//' destinations from given set of origins, with flows attenuated by distance from
//' those origins.
//'
//' @param graph The data.frame holding the graph edges
//' @param vert_map_in map from <std::string> vertex ID to (0-indexed) integer
//' index of vertices
//' @param fromi Index into vert_map_in of vertex numbers
//' @param k Coefficient of (current proof-of-principle-only) exponential
//' distance decay function.  If value of \code{k<0} is given, a standard
//' logistic polynomial will be used.
//'
//' @note The flow data to be used for aggregation is a matrix mapping flows
//' betwen each pair of from and to points.
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_flows_disperse (const Rcpp::DataFrame graph,
        const Rcpp::DataFrame vert_map_in,
        Rcpp::IntegerVector fromi,
        double k,
        Rcpp::NumericMatrix flows,
        std::string heap_type)
{
    Rcpp::NumericVector id_vec;
    size_t nfrom = get_fromi (vert_map_in, fromi, id_vec);

    std::vector <std::string> from = graph ["from"];
    std::vector <std::string> to = graph ["to"];
    std::vector <double> dist = graph ["d"];
    std::vector <double> wt = graph ["w"];

    unsigned int nedges = static_cast <unsigned int> (graph.nrow ());
    std::vector <std::string> vert_name = vert_map_in ["vert"];
    std::vector <unsigned int> vert_indx = vert_map_in ["id"];
    // Make map from vertex name to integer index
    std::map <std::string, unsigned int> vert_map_i;
    size_t nverts = make_vert_map (vert_map_in, vert_name,
            vert_indx, vert_map_i);

    std::unordered_map <std::string, unsigned int> verts_to_edge_map;
    std::unordered_map <std::string, double> verts_to_dist_map;
    make_vert_to_edge_maps (from, to, wt, verts_to_edge_map, verts_to_dist_map);

    std::shared_ptr<DGraph> g = std::make_shared<DGraph> (nverts);
    inst_graph (g, nedges, vert_map_i, from, to, dist, wt);

    std::shared_ptr <Dijkstra> dijkstra = std::make_shared <Dijkstra> (nverts,
            *getHeapImpl(heap_type), g);

    Rcpp::List res (nfrom);
    std::vector<double> w(nverts);
    std::vector<double> d(nverts);
    std::vector<int> prev(nverts);

    Rcpp::NumericVector aggregate_flows (from.size ()); // 0-filled by default
    for (unsigned int v = 0; v < nfrom; v++)
    {
        Rcpp::checkUserInterrupt ();
        std::fill (w.begin(), w.end(), INFINITE_DOUBLE);
        std::fill (d.begin(), d.end(), INFINITE_DOUBLE);

        dijkstra->run (d, w, prev, static_cast <unsigned int> (fromi [v]));

        for (unsigned int vi = 0; vi < nverts; vi++)
        {
            if (prev [vi] > 0)
            {
                // NOTE: Critically important that these are in the right order!
                std::string vert_to = vert_name [vi],
                    vert_from = vert_name [prev [vi]];
                std::string two_verts = "f" + vert_from + "t" + vert_to;
                if (verts_to_edge_map.find (two_verts) == verts_to_edge_map.end ())
                    Rcpp::stop ("vertex pair forms no known edge");

                unsigned int indx = verts_to_edge_map [two_verts];
                if (d [vi] != INFINITE_DOUBLE)
                {
                    if (k > 0.0)
                        aggregate_flows [indx] += flows (v, 0) * exp (-d [vi] / k);
                    else // standard logistic polynomial for UK cycling models
                    {
                        double lp = -3.894 + (-0.5872 * d [vi]) +
                            (1.832 * sqrt (d [vi])) +
                            (0.007956 * d [vi] * d [vi]);
                        aggregate_flows [indx] += flows (v, 0) *
                            exp (lp) / (1.0 + exp (lp));
                    }
                }
            }
        }
    } // end for v over nfrom

    return (aggregate_flows);
}

//' rcpp_spatial_interaction
//'
//' Singly constrained spatial interaction model using exponential interaction
//' function. Given an input vector of site densities, this function maps these
//' onto full 2D spatial interaction terms (origin-destination) by calculating
//' values of \code{T_i T_j exp (-k d_jk)} for each site, \code{i}, normalised
//' by total sums for that site (\code{\sum_j}).
//'
//' @param graph The data.frame holding the graph edges
//' @param vert_map_in map from <std::string> vertex ID to (0-indexed) integer
//' index of vertices
//' @param nodes Index into vert_map_in of vertex numbers
//' @param k Coefficient of exponential spatial interaction function.
//' @param dens Vector of densities of same size as both \code{vert_map_in}
//' and \code{fromi}.
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_spatial_interaction (const Rcpp::DataFrame graph,
        const Rcpp::DataFrame vert_map_in,
        Rcpp::IntegerVector nodes,
        double k,
        Rcpp::NumericVector dens,
        std::string heap_type)
{
    Rcpp::NumericVector id_vec;
    size_t nfrom = get_fromi (vert_map_in, nodes, id_vec);

    std::vector <std::string> from = graph ["from"];
    std::vector <std::string> to = graph ["to"];
    std::vector <double> dist = graph ["d"];
    std::vector <double> wt = graph ["w"];

    unsigned int nedges = static_cast <unsigned int> (graph.nrow ());
    std::vector <std::string> vert_name = vert_map_in ["vert"];
    std::vector <unsigned int> vert_indx = vert_map_in ["id"];
    // Make map from vertex name to integer index
    std::map <std::string, unsigned int> vert_map_i;
    size_t nverts = make_vert_map (vert_map_in, vert_name,
            vert_indx, vert_map_i);

    std::unordered_map <std::string, unsigned int> verts_to_edge_map;
    std::unordered_map <std::string, double> verts_to_dist_map;
    make_vert_to_edge_maps (from, to, wt, verts_to_edge_map, verts_to_dist_map);

    std::shared_ptr<DGraph> g = std::make_shared<DGraph> (nverts);
    inst_graph (g, nedges, vert_map_i, from, to, dist, wt);

    std::shared_ptr <Dijkstra> dijkstra = std::make_shared <Dijkstra> (nverts,
            *getHeapImpl(heap_type), g);

    Rcpp::List res (nfrom);
    std::vector<double> w(nverts);
    std::vector<double> d(nverts);
    std::vector<int> prev(nverts);

    Rcpp::NumericMatrix SI (nodes.size (), nodes.size ());
    for (unsigned int v = 0; v < nfrom; v++)
    {
        Rcpp::checkUserInterrupt ();
        std::fill (w.begin(), w.end(), INFINITE_DOUBLE);
        std::fill (d.begin(), d.end(), INFINITE_DOUBLE);

        dijkstra->run (d, w, prev, static_cast <unsigned int> (nodes [v]));

        double flowsums = 0.0;
        for (unsigned int vi = 0; vi < nfrom; vi++)
        {
            // Doesn't matter if d [] == INFINITE_INT
            const double tempd = dens (v) * dens (vi) * exp (-d [nodes [vi]] / k);
            SI (v, vi) = tempd;
            flowsums += tempd;
        }
        if (flowsums == 0.0)
        {
            for (unsigned int vi = 0; vi < nfrom; vi++)
                SI (v, vi) = 0.0;
        } else
        {
            for (unsigned int vi = 0; vi < nfrom; vi++)
                SI (v, vi) = dens (v) * SI (v, vi) / flowsums;
        }
    } // end for v over nfrom

    return (SI);
}
