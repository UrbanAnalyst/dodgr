
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


std::shared_ptr <HeapDesc> run_sp::getHeapImpl(const std::string& heap_type)
{
  if (heap_type == "FHeap")
    return std::make_shared <HeapD<FHeap> >();
  else if (heap_type == "BHeap" || heap_type == "set") // heap not used for set 
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
    throw std::runtime_error("invalid heap type: " + heap_type); // # nocov
}


struct OneDist : public RcppParallel::Worker
{
    RcppParallel::RVector <int> dp_fromi;
    const std::vector <unsigned int> toi;
    const size_t nverts;
    const std::vector <double> vx;
    const std::vector <double> vy;
    const std::shared_ptr <DGraph> g;
    const std::string heap_type;
    bool is_spatial;

    RcppParallel::RMatrix <double> dout;

    // constructor
    OneDist (
            const Rcpp::IntegerVector fromi,
            //const Rcpp::IntegerVector toi_in,
            const std::vector <unsigned int> toi_in,
            const size_t nverts_in,
            const std::vector <double> vx_in,
            const std::vector <double> vy_in,
            const std::shared_ptr <DGraph> g_in,
            const std::string & heap_type_in,
            const bool & is_spatial_in,
            Rcpp::NumericMatrix dout_in) :
        dp_fromi (fromi), toi (toi_in), nverts (nverts_in),
        vx (vx_in), vy (vy_in),
        g (g_in), heap_type (heap_type_in), is_spatial (is_spatial_in),
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
        std::vector <double> d (nverts);
        std::vector <int> prev (nverts);

        std::vector <double> heuristic (nverts, 0.0);

        for (std::size_t i = begin; i < end; i++)
        {
            unsigned int from_i = static_cast <unsigned int> (dp_fromi [i]);
            // These have to be reserved within the parallel operator function!
            std::fill (w.begin (), w.end (), INFINITE_DOUBLE);
            std::fill (d.begin (), d.end (), INFINITE_DOUBLE);

            if (is_spatial)
            {
                double dmax = 0.0;
                for (size_t j = 0; j < nverts; j++)
                {
                    double dx = vx [j] - vx [from_i],
                        dy = vy [j] - vy [from_i];
                    heuristic [j] = sqrt (dx * dx + dy * dy);
                    if (heuristic [j] > dmax)
                        dmax = heuristic [j];
                }
                pathfinder->AStar (d, w, prev, heuristic, from_i, toi);
            } else if (heap_type.find ("set") == std::string::npos)
                pathfinder->Dijkstra (d, w, prev, from_i, toi);
            else
                pathfinder->Dijkstra_set (d, w, prev, from_i);

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


struct OneIso : public RcppParallel::Worker
{
    RcppParallel::RVector <int> dp_fromi;
    const size_t nverts;
    const std::shared_ptr <DGraph> g;
    const Rcpp::NumericVector dlimit;
    const std::string heap_type;

    RcppParallel::RMatrix <double> dout;

    // constructor
    OneIso (
            const Rcpp::IntegerVector fromi,
            const size_t nverts_in,
            const std::shared_ptr <DGraph> g_in,
            const Rcpp::NumericVector dlimit_in,
            const std::string & heap_type_in,
            Rcpp::NumericMatrix dout_in) :
        dp_fromi (fromi), nverts (nverts_in),
        g (g_in), dlimit (dlimit_in),
        heap_type (heap_type_in), dout (dout_in)
    {
    }

    // Parallel function operator
    void operator() (std::size_t begin, std::size_t end)
    {
        std::shared_ptr<PF::PathFinder> pathfinder =
            std::make_shared <PF::PathFinder> (nverts,
                    *run_sp::getHeapImpl (heap_type), g);
        std::vector <double> w (nverts);
        std::vector <double> d (nverts);
        std::vector <int> prev (nverts);

        double dlimit_max = *std::max_element (dlimit.begin (), dlimit.end ());

        for (std::size_t i = begin; i < end; i++)
        {
            unsigned int from_i = static_cast <unsigned int> (dp_fromi [i]);
            std::fill (w.begin (), w.end (), INFINITE_DOUBLE);
            std::fill (d.begin (), d.end (), INFINITE_DOUBLE);

            pathfinder->DijkstraLimit (d, w, prev, from_i, dlimit_max);

            // Get the set of terminal vertices.
            std::set <int> terminal_verts;
            for (int j: prev)
            {
                size_t sj = static_cast <size_t> (j);
                if (w [sj] > INFINITE_INT && w [sj] < dlimit_max)
                {
                    // # nocov start
                    if (prev [sj] < INFINITE_INT)
                    {
                        terminal_verts.erase (j);
                        terminal_verts.emplace (prev [sj]);
                    }
                    // # nocov end
                }
            }


            for (size_t j = 0; j < nverts; j++)
            {
                if (terminal_verts.find (static_cast <int> (j)) !=
                        terminal_verts.end ())
                {
                    // Flag terminal verts with -d
                    dout (i, j) = -w [j]; // # nocov
                } else if (prev [j] < INFINITE_INT && w [j] < dlimit_max)
                {
                    for (auto k: dlimit)
                    {
                        size_t sp = static_cast <size_t> (prev [j]);
                        if (w [j] > k && w[sp] < k)
                            dout (i, sp) = -k; // flag isohull verts with -k
                        else
                            dout (i, j) = w [j]; // distance of other internal verts
                    }
                }
            }
        }
    }
                                   
};


size_t run_sp::make_vert_map (const Rcpp::DataFrame &vert_map_in,
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

// Flows from the pathfinder output are reallocated based on matching vertex
// pairs to edge indices. Note, however, that contracted graphs frequently
// have duplicate vertex pairs with different distances. The following
// therefore uses two maps, one to hold the ultimate index from vertex
// pairs, and the other to hold minimal distances. This is used in flow routines
// only.
void run_sp::make_vert_to_edge_maps (const std::vector <std::string> &from,
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
        Rcpp::IntegerVector toi_in,
        const std::string& heap_type,
        const bool is_spatial)
{
    std::vector <unsigned int> toi =
        Rcpp::as <std::vector <unsigned int> > ( toi_in);

    Rcpp::NumericVector id_vec;
    size_t nfrom = static_cast <size_t> (fromi.size ());
    size_t nto = static_cast <size_t> (toi.size ());

    std::vector <std::string> from = graph ["from"];
    std::vector <std::string> to = graph ["to"];
    std::vector <double> dist = graph ["d"];
    std::vector <double> wt = graph ["w"];

    unsigned int nedges = static_cast <unsigned int> (graph.nrow ());
    std::map <std::string, unsigned int> vert_map;
    std::vector <std::string> vert_map_id = vert_map_in ["vert"];
    std::vector <unsigned int> vert_map_n = vert_map_in ["id"];
    size_t nverts = run_sp::make_vert_map (vert_map_in, vert_map_id,
            vert_map_n, vert_map);

    std::vector <double> vx (nverts), vy (nverts);
    if (is_spatial)
    {
        vx = Rcpp::as <std::vector <double> > (vert_map_in ["x"]);
        vy = Rcpp::as <std::vector <double> > (vert_map_in ["y"]);
    }

    std::shared_ptr <DGraph> g = std::make_shared <DGraph> (nverts);
    inst_graph (g, nedges, vert_map, from, to, dist, wt);

    Rcpp::NumericVector na_vec = Rcpp::NumericVector (nfrom * nto,
            Rcpp::NumericVector::get_na ());
    Rcpp::NumericMatrix dout (static_cast <int> (nfrom),
            static_cast <int> (nto), na_vec.begin ());

    // Create parallel worker
    OneDist one_dist (fromi, toi, nverts, vx, vy,
            g, heap_type, is_spatial, dout);

    RcppParallel::parallelFor (0, static_cast <size_t> (fromi.length ()),
            one_dist);
    
    return (dout);
}

//' rcpp_get_iso
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_get_iso (const Rcpp::DataFrame graph,
        const Rcpp::DataFrame vert_map_in,
        Rcpp::IntegerVector fromi,
        Rcpp::NumericVector dlim,
        const std::string& heap_type)
{
    Rcpp::NumericVector id_vec;
    const size_t nfrom = static_cast <size_t> (fromi.size ());

    std::vector <std::string> from = graph ["from"];
    std::vector <std::string> to = graph ["to"];
    std::vector <double> dist = graph ["d"];
    std::vector <double> wt = graph ["w"];

    unsigned int nedges = static_cast <unsigned int> (graph.nrow ());
    std::map <std::string, unsigned int> vert_map;
    std::vector <std::string> vert_map_id = vert_map_in ["vert"];
    std::vector <unsigned int> vert_map_n = vert_map_in ["id"];
    size_t nverts = run_sp::make_vert_map (vert_map_in, vert_map_id,
            vert_map_n, vert_map);

    std::vector <double> vx (nverts), vy (nverts);

    std::shared_ptr <DGraph> g = std::make_shared <DGraph> (nverts);
    inst_graph (g, nedges, vert_map, from, to, dist, wt);

    Rcpp::NumericVector na_vec = Rcpp::NumericVector (nfrom * nverts,
            Rcpp::NumericVector::get_na ());
    Rcpp::NumericMatrix dout (static_cast <int> (nfrom),
            static_cast <int> (nverts), na_vec.begin ());

    // Create parallel worker
    OneIso one_iso (fromi, nverts, g, dlim, heap_type, dout);

    RcppParallel::parallelFor (0, static_cast <size_t> (fromi.length ()),
            one_iso);
    
    return (dout);
}

//' rcpp_get_sp_dists
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_get_sp_dists (const Rcpp::DataFrame graph,
        const Rcpp::DataFrame vert_map_in,
        Rcpp::IntegerVector fromi,
        Rcpp::IntegerVector toi_in,
        const std::string& heap_type)
{
    std::vector <unsigned int> toi =
        Rcpp::as <std::vector <unsigned int> > ( toi_in);
    Rcpp::NumericVector id_vec;
    size_t nfrom = static_cast <size_t> (fromi.size ());
    size_t nto = static_cast <size_t> (toi.size ());

    std::vector <std::string> from = graph ["from"];
    std::vector <std::string> to = graph ["to"];
    std::vector <double> dist = graph ["d"];
    std::vector <double> wt = graph ["w"];

    unsigned int nedges = static_cast <unsigned int> (graph.nrow ());
    std::map <std::string, unsigned int> vert_map;
    std::vector <std::string> vert_map_id = vert_map_in ["vert"];
    std::vector <unsigned int> vert_map_n = vert_map_in ["id"];
    size_t nverts = run_sp::make_vert_map (vert_map_in, vert_map_id,
            vert_map_n, vert_map);

    std::shared_ptr<DGraph> g = std::make_shared<DGraph>(nverts);
    inst_graph (g, nedges, vert_map, from, to, dist, wt);

    std::shared_ptr <PF::PathFinder> pathfinder =
        std::make_shared <PF::PathFinder> (
            nverts, *run_sp::getHeapImpl(heap_type), g);

    std::vector<double> w (nverts);
    std::vector<double> d (nverts);
    std::vector<int> prev (nverts);

    pathfinder->init (g); // specify the graph

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

        pathfinder->Dijkstra (d, w, prev,
                static_cast <unsigned int> (fromi [v]), toi);
        for (unsigned int vi = 0; vi < nto; vi++)
        {
            if (w [static_cast <size_t> (toi [vi])] < INFINITE_DOUBLE)
            {
                dout (v, vi) = d [static_cast <size_t> (toi [vi])];
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
//' @note Returns 1-indexed values indexing directly into the R input
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::List rcpp_get_paths (const Rcpp::DataFrame graph,
        const Rcpp::DataFrame vert_map_in,
        Rcpp::IntegerVector fromi,
        Rcpp::IntegerVector toi_in,
        const std::string& heap_type)
{
    std::vector <unsigned int> toi =
        Rcpp::as <std::vector <unsigned int> > ( toi_in);
    Rcpp::NumericVector id_vec;
    size_t nfrom = static_cast <size_t> (fromi.size ());
    size_t nto = static_cast <size_t> (toi.size ());

    std::vector <std::string> from = graph ["from"];
    std::vector <std::string> to = graph ["to"];
    std::vector <double> dist = graph ["d"];
    std::vector <double> wt = graph ["w"];

    unsigned int nedges = static_cast <unsigned int> (graph.nrow ());
    std::map <std::string, unsigned int> vert_map;
    std::vector <std::string> vert_map_id = vert_map_in ["vert"];
    std::vector <unsigned int> vert_map_n = vert_map_in ["id"];
    size_t nverts = run_sp::make_vert_map (vert_map_in, vert_map_id,
            vert_map_n, vert_map);

    std::shared_ptr<DGraph> g = std::make_shared<DGraph>(nverts);
    inst_graph (g, nedges, vert_map, from, to, dist, wt);

    std::shared_ptr<PF::PathFinder> pathfinder =
        std::make_shared <PF::PathFinder> (nverts,
            *run_sp::getHeapImpl(heap_type), g);
    
    Rcpp::List res (nfrom);
    std::vector<double> w (nverts);
    std::vector<double> d (nverts);
    std::vector<int> prev (nverts);

    pathfinder->init (g); // specify the graph

    for (unsigned int v = 0; v < nfrom; v++)
    {
        Rcpp::checkUserInterrupt ();
        std::fill (w.begin(), w.end(), INFINITE_DOUBLE);
        std::fill (d.begin(), d.end(), INFINITE_DOUBLE);

        pathfinder->Dijkstra (d, w, prev,
                static_cast <unsigned int> (fromi [v]), toi);

        Rcpp::List res1 (nto);
        for (unsigned int vi = 0; vi < nto; vi++)
        {
            std::vector <int> onePath;
            if (w [toi [vi]] < INFINITE_DOUBLE)
            {
                int target = toi_in [vi]; // target can be -1!
                while (target < INFINITE_INT)
                {
                    // Note that targets are all C++ 0-indexed and are converted
                    // directly here to R-style 1-indexes.
                    onePath.push_back (target + 1);
                    target = static_cast <int> (prev [static_cast <unsigned int> (target)]);
                    if (target < 0 || target == fromi [v])
                        break;
                }
            }
            if (onePath.size () >= 1)
            {
                onePath.push_back (fromi [v] + 1);
                std::reverse (onePath.begin (), onePath.end ());
                res1 [vi] = onePath;
            }
        }
        res [v] = res1;
    }
    return (res);
}
