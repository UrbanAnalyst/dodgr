
#include "run_sp.h"

#include "dgraph.h"
#include "heaps/heap_lib.h"

// # nocov start
template <typename T>
void inst_graph (std::shared_ptr<DGraph> g, size_t nedges,
        const std::map <std::string, size_t>& vert_map,
        const std::vector <std::string>& from,
        const std::vector <std::string>& to,
        const std::vector <T>& dist,
        const std::vector <T>& wt)
{
    for (size_t i = 0; i < nedges; ++i)
    {
        size_t fromi = vert_map.at(from [i]);
        size_t toi = vert_map.at(to [i]);
        g->addNewEdge (fromi, toi, dist [i], wt [i], i);
    }
}
// # nocov end

// RcppParallel jobs can be chunked to a specified "grain size"; see
// https://rcppcore.github.io/RcppParallel/#grain_size
// This function determines chunk size such that there are at least 100 chunks
// for a given `nfrom`.
size_t run_sp::get_chunk_size (const size_t nfrom)
{
    size_t chunk_size;

    if (nfrom > 1000)
        chunk_size = 100;
    else if (nfrom > 100)
        chunk_size = 10;
    else
        chunk_size = 1;

    return chunk_size;
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
  else
    throw std::runtime_error("invalid heap type: " + heap_type); // # nocov
}


struct OneDist : public RcppParallel::Worker
{
    RcppParallel::RVector <int> dp_fromi;
    const std::vector <size_t> toi;
    const size_t nverts;
    const std::vector <double> vx;
    const std::vector <double> vy;
    const std::shared_ptr <DGraph> g;
    const std::string heap_type;
    bool is_spatial;

    RcppParallel::RMatrix <double> dout;

    // constructor
    OneDist (
            const RcppParallel::RVector <int> fromi,
            const std::vector <size_t> toi_in,
            const size_t nverts_in,
            const std::vector <double> vx_in,
            const std::vector <double> vy_in,
            const std::shared_ptr <DGraph> g_in,
            const std::string & heap_type_in,
            const bool & is_spatial_in,
            RcppParallel::RMatrix <double> dout_in) :
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
        std::vector <long int> prev (nverts);

        std::vector <double> heuristic (nverts, 0.0);

        for (std::size_t i = begin; i < end; i++)
        {
            size_t from_i = static_cast <size_t> (dp_fromi [i]);

            if (is_spatial)
            {
                for (size_t j = 0; j < nverts; j++)
                {
                    const double dx = vx [j] - vx [from_i],
                        dy = vy [j] - vy [from_i];
                    heuristic [j] = sqrt (dx * dx + dy * dy);
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

struct OneDistPaired : public RcppParallel::Worker
{
    RcppParallel::RVector <int> dp_fromtoi;
    const size_t nverts;
    const size_t nfrom;
    const std::vector <double> vx;
    const std::vector <double> vy;
    const std::shared_ptr <DGraph> g;
    const std::string heap_type;
    bool is_spatial;

    RcppParallel::RMatrix <double> dout;

    // constructor
    OneDistPaired (
            const RcppParallel::RVector <int> fromtoi,
            const size_t nverts_in,
            const size_t nfrom_in,
            const std::vector <double> vx_in,
            const std::vector <double> vy_in,
            const std::shared_ptr <DGraph> g_in,
            const std::string & heap_type_in,
            const bool & is_spatial_in,
            RcppParallel::RMatrix <double> dout_in) :
        dp_fromtoi (fromtoi), nverts (nverts_in), nfrom (nfrom_in),
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
        std::vector <long int> prev (nverts);

        std::vector <double> heuristic (nverts, 0.0);

        for (std::size_t i = begin; i < end; i++)
        {
            const size_t from_i = static_cast <size_t> (dp_fromtoi [i]);
            const std::vector <size_t> to_i = {static_cast <size_t> (dp_fromtoi [nfrom + i])};

            if (is_spatial)
            {
                // need to set an additional target vertex that is somewhat
                // beyond the single actual target vertex. Default here is max
                // heuristic, but reduced in following loop.
                long int max_h_index = -1;
                double max_h_value = -1.0;
                for (size_t j = 0; j < nverts; j++)
                {
                    const double dx = vx [j] - vx [from_i],
                        dy = vy [j] - vy [from_i];
                    heuristic [j] = sqrt (dx * dx + dy * dy);
                    if (heuristic [j] > max_h_value) {
                        max_h_value = heuristic [j];
                        max_h_index = static_cast <long int> (j);
                    }
                }
                const double htemp = heuristic [static_cast <size_t> (dp_fromtoi [nfrom + i])];
                double min_h_value = max_h_value;
                long int min_h_index = max_h_index;
                // Arbitrary relative distance threshold
                // TODO: Are there likely to be cases where this might need to
                // be adjusted?
                const double thr = 0.1;
                for (size_t j = 0; j < nverts; j++) {
                    if ((heuristic [j] < (thr * htemp)) && (heuristic [j] > min_h_value)) {
                        min_h_value = heuristic [j];
                        min_h_index = static_cast <long int> (j);
                    }
                }
                const std::vector <size_t> to_i2 = {to_i [0], static_cast <size_t> (min_h_index)};
                pathfinder->AStar (d, w, prev, heuristic, from_i, to_i2);
            } else if (heap_type.find ("set") == std::string::npos)
                pathfinder->Dijkstra (d, w, prev, from_i, to_i);
            else
                pathfinder->Dijkstra_set (d, w, prev, from_i);

            if (w [to_i [0]] < INFINITE_DOUBLE)
                dout (i, 0) = d [to_i [0]];
        }
    }
                                   
};


struct OneIso : public RcppParallel::Worker
{
    RcppParallel::RVector <int> dp_fromi;
    const size_t nverts;
    const std::shared_ptr <DGraph> g;
    const RcppParallel::RVector <double> dlimit;
    const std::string heap_type;

    RcppParallel::RMatrix <double> dout;

    // constructor
    OneIso (
            const RcppParallel::RVector <int> fromi,
            const size_t nverts_in,
            const std::shared_ptr <DGraph> g_in,
            const RcppParallel::RVector <double> dlimit_in,
            const std::string & heap_type_in,
            RcppParallel::RMatrix <double> dout_in) :
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
        std::vector <long int> prev (nverts);

        const double dlimit_max = *std::max_element (dlimit.begin (), dlimit.end ());

        for (std::size_t i = begin; i < end; i++)
        {
            std::fill (w.begin (), w.end (), INFINITE_DOUBLE);
            std::fill (d.begin (), d.end (), INFINITE_DOUBLE);
            std::fill (prev.begin (), prev.end (), INFINITE_INT);

            size_t from_i = static_cast <size_t> (dp_fromi [i]);

            pathfinder->DijkstraLimit (d, w, prev, from_i, dlimit_max);

            // Get the set of terminal vertices: those with w < dlimit_max but 
            // with no previous (outward-going) nodes
            std::unordered_set <int> terminal_verts;
            for (size_t j = 0; j < nverts; j++)
            {
                if (prev [j] == INFINITE_INT && d [j] < dlimit_max)
                {
                    terminal_verts.emplace (static_cast <int> (j));
                }
            }


            for (size_t j = 0; j < nverts; j++)
            {
                if (terminal_verts.find (static_cast <int> (j)) !=
                        terminal_verts.end ())
                {
                    // Flag terminal verts with -d
                    dout (i, j) = -d [j]; // # nocov
                } else if (prev [j] < INFINITE_INT && d [j] < dlimit_max)
                {
                    size_t st_prev = static_cast <size_t> (prev [j]);
                    for (auto k: dlimit)
                    {
                        if (d [j] > k && d [st_prev] < k)
                            dout (i, st_prev) = -k; // flag isohull verts with -k
                        else
                            dout (i, j) = d [j]; // distance of other internal verts
                    }
                }
            }
        }
    }
                                   
};


size_t run_sp::make_vert_map (const Rcpp::DataFrame &vert_map_in,
        const std::vector <std::string> &vert_map_id,
        const std::vector <size_t> &vert_map_n,
        std::map <std::string, size_t> &vert_map)
{
    for (size_t i = 0;
            i < static_cast <size_t> (vert_map_in.nrow ()); ++i)
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
        std::unordered_map <std::string, size_t> &verts_to_edge_map,
        std::unordered_map <std::string, double> &verts_to_dist_map)
{
    for (size_t i = 0; i < from.size (); i++)
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
    std::vector <size_t> toi =
        Rcpp::as <std::vector <size_t> > ( toi_in);

    size_t nfrom = static_cast <size_t> (fromi.size ());
    size_t nto = static_cast <size_t> (toi.size ());

    const std::vector <std::string> from = graph ["from"];
    const std::vector <std::string> to = graph ["to"];
    const std::vector <double> dist = graph ["d"];
    const std::vector <double> wt = graph ["d_weighted"];

    const size_t nedges = static_cast <size_t> (graph.nrow ());
    std::map <std::string, size_t> vert_map;
    std::vector <std::string> vert_map_id = vert_map_in ["vert"];
    std::vector <size_t> vert_map_n = vert_map_in ["id"];
    const size_t nverts = run_sp::make_vert_map (vert_map_in, vert_map_id,
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
    OneDist one_dist (RcppParallel::RVector <int> (fromi), toi,
            nverts, vx, vy, g, heap_type, is_spatial,
            RcppParallel::RMatrix <double> (dout));

    size_t chunk_size = run_sp::get_chunk_size (nfrom);
    RcppParallel::parallelFor (0, nfrom, one_dist, chunk_size);
    
    return (dout);
}

//' rcpp_get_sp_dists_par
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_get_sp_dists_paired_par (const Rcpp::DataFrame graph,
        const Rcpp::DataFrame vert_map_in,
        Rcpp::IntegerVector fromi,
        Rcpp::IntegerVector toi,
        const std::string& heap_type,
        const bool is_spatial)
{
    if (fromi.size () != toi.size ())
        Rcpp::stop ("pairwise dists must have from.size == to.size");
    long int n = fromi.size ();
    size_t n_st = static_cast <size_t> (n);

    const std::vector <std::string> from = graph ["from"];
    const std::vector <std::string> to = graph ["to"];
    const std::vector <double> dist = graph ["d"];
    const std::vector <double> wt = graph ["d_weighted"];

    const size_t nedges = static_cast <size_t> (graph.nrow ());
    std::map <std::string, size_t> vert_map;
    std::vector <std::string> vert_map_id = vert_map_in ["vert"];
    std::vector <size_t> vert_map_n = vert_map_in ["id"];
    const size_t nverts = run_sp::make_vert_map (vert_map_in, vert_map_id,
            vert_map_n, vert_map);

    std::vector <double> vx (nverts), vy (nverts);
    if (is_spatial)
    {
        vx = Rcpp::as <std::vector <double> > (vert_map_in ["x"]);
        vy = Rcpp::as <std::vector <double> > (vert_map_in ["y"]);
    }

    std::shared_ptr <DGraph> g = std::make_shared <DGraph> (nverts);
    inst_graph (g, nedges, vert_map, from, to, dist, wt);

    Rcpp::NumericVector na_vec = Rcpp::NumericVector (n_st,
            Rcpp::NumericVector::get_na ());
    Rcpp::NumericMatrix dout (static_cast <int> (n), 1, na_vec.begin ());

    // Paired (fromi, toi) in a single vector
    Rcpp::IntegerVector fromto (2 * n_st);
    for (int i = 0; i < n; i++)
    {
        size_t i_t = static_cast <size_t> (i);
        fromto [i] = fromi (i_t);
        fromto [i + n] = toi (i_t);
    }

    // Create parallel worker
    OneDistPaired one_dist_paired (RcppParallel::RVector <int> (fromto),
            nverts, n_st, vx, vy, g, heap_type, is_spatial,
            RcppParallel::RMatrix <double> (dout));

    size_t chunk_size = run_sp::get_chunk_size (n_st);
    RcppParallel::parallelFor (0, n_st, one_dist_paired, chunk_size);
    
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
    const size_t nfrom = static_cast <size_t> (fromi.size ());

    std::vector <std::string> from = graph ["from"];
    std::vector <std::string> to = graph ["to"];
    std::vector <double> dist = graph ["d"];
    std::vector <double> wt = graph ["d_weighted"];

    size_t nedges = static_cast <size_t> (graph.nrow ());
    std::map <std::string, size_t> vert_map;
    std::vector <std::string> vert_map_id = vert_map_in ["vert"];
    std::vector <size_t> vert_map_n = vert_map_in ["id"];
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
    OneIso one_iso (RcppParallel::RVector <int> (fromi), nverts, g,
            RcppParallel::RVector <double> (dlim), heap_type,
            RcppParallel::RMatrix <double> (dout));

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
    std::vector <size_t> toi =
        Rcpp::as <std::vector <size_t> > ( toi_in);
    size_t nfrom = static_cast <size_t> (fromi.size ());
    size_t nto = static_cast <size_t> (toi.size ());

    std::vector <std::string> from = graph ["from"];
    std::vector <std::string> to = graph ["to"];
    std::vector <double> dist = graph ["d"];
    std::vector <double> wt = graph ["d_weighted"];

    size_t nedges = static_cast <size_t> (graph.nrow ());
    std::map <std::string, size_t> vert_map;
    std::vector <std::string> vert_map_id = vert_map_in ["vert"];
    std::vector <size_t> vert_map_n = vert_map_in ["id"];
    size_t nverts = run_sp::make_vert_map (vert_map_in, vert_map_id,
            vert_map_n, vert_map);

    std::shared_ptr<DGraph> g = std::make_shared<DGraph>(nverts);
    inst_graph (g, nedges, vert_map, from, to, dist, wt);

    std::vector<double> w (nverts);
    std::vector<double> d (nverts);
    std::vector<long int> prev (nverts);

    // initialise dout matrix to NA
    Rcpp::NumericVector na_vec = Rcpp::NumericVector (nfrom * nto,
            Rcpp::NumericVector::get_na ());
    Rcpp::NumericMatrix dout (static_cast <int> (nfrom),
            static_cast <int> (nto), na_vec.begin ());


    for (size_t i = 0; i < nfrom; i++)
    {
        // These lines (re-)initialise the heap, so have to be called for each v
        std::shared_ptr <PF::PathFinder> pathfinder =
            std::make_shared <PF::PathFinder> (
                nverts, *run_sp::getHeapImpl(heap_type), g);

        pathfinder->init (g); // specify the graph

        Rcpp::checkUserInterrupt ();
        std::fill (w.begin(), w.end(), INFINITE_DOUBLE);
        std::fill (d.begin(), d.end(), INFINITE_DOUBLE);
        size_t fromi_i = static_cast <size_t> (fromi [static_cast <R_xlen_t> (i)]);
        d [fromi_i] = w [fromi_i] = 0.0;

        pathfinder->Dijkstra (d, w, prev, fromi_i, toi);
        for (size_t j = 0; j < nto; j++)
        {
            if (w [static_cast <size_t> (toi [j])] < INFINITE_DOUBLE)
            {
                dout (i, j) = d [static_cast <size_t> (toi [j])];
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
    std::vector <size_t> toi =
        Rcpp::as <std::vector <size_t> > ( toi_in);
    size_t nfrom = static_cast <size_t> (fromi.size ());
    size_t nto = static_cast <size_t> (toi.size ());

    std::vector <std::string> from = graph ["from"];
    std::vector <std::string> to = graph ["to"];
    std::vector <double> dist = graph ["d"];
    std::vector <double> wt = graph ["d_weighted"];

    size_t nedges = static_cast <size_t> (graph.nrow ());
    std::map <std::string, size_t> vert_map;
    std::vector <std::string> vert_map_id = vert_map_in ["vert"];
    std::vector <size_t> vert_map_n = vert_map_in ["id"];
    size_t nverts = run_sp::make_vert_map (vert_map_in, vert_map_id,
            vert_map_n, vert_map);

    std::shared_ptr<DGraph> g = std::make_shared<DGraph>(nverts);
    inst_graph (g, nedges, vert_map, from, to, dist, wt);

    Rcpp::List res (nfrom);
    std::vector<double> w (nverts);
    std::vector<double> d (nverts);
    std::vector<long int> prev (nverts);

    for (size_t i = 0; i < nfrom; i++)
    {
        const R_xlen_t i_R = static_cast <R_xlen_t> (i);

        // These lines (re-)initialise the heap, so have to be called for each i
        std::shared_ptr<PF::PathFinder> pathfinder =
            std::make_shared <PF::PathFinder> (nverts,
                *run_sp::getHeapImpl(heap_type), g);
        
        pathfinder->init (g); // specify the graph

        Rcpp::checkUserInterrupt ();
        std::fill (w.begin(), w.end(), INFINITE_DOUBLE);
        std::fill (d.begin(), d.end(), INFINITE_DOUBLE);
        std::fill (prev.begin(), prev.end(), INFINITE_INT);
        d [static_cast <size_t> (fromi [i_R])] =
            w [static_cast <size_t> (fromi [i_R])] = 0.0;

        pathfinder->Dijkstra (d, w, prev,
                static_cast <size_t> (fromi [i_R]), toi);

        Rcpp::List res1 (nto);
        for (size_t j = 0; j < nto; j++)
        {
            std::vector <long int> onePath;
            if (w [toi [j]] < INFINITE_DOUBLE)
            {
                long int target = toi_in [static_cast <R_xlen_t> (j)]; // target can be -1!
                while (target < INFINITE_INT)
                {
                    // Note that targets are all C++ 0-indexed and are converted
                    // directly here to R-style 1-indexes.
                    onePath.push_back (target + 1L);
                    target = prev [static_cast <size_t> (target)];
                    if (target < 0L || target == fromi [i_R])
                        break;
                }
            }
            if (onePath.size () >= 1)
            {
                onePath.push_back (static_cast <R_xlen_t> (fromi [i_R] + 1L));
                std::reverse (onePath.begin (), onePath.end ());
                res1 [static_cast <R_xlen_t> (j)] = onePath;
            }
        }
        res [static_cast <R_xlen_t> (i)] = res1;
    }
    return (res);
}
