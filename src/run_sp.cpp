
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
    throw std::runtime_error("invalid heap type: " + heap_type);
}


struct OneDist : public RcppParallel::Worker
{
    RcppParallel::RVector <int> dp_fromi;
    const Rcpp::IntegerVector toi;
    const size_t nverts;
    const std::shared_ptr <DGraph> g;
    const std::string heap_type;

    RcppParallel::RMatrix <double> dout;

    // constructor
    OneDist (
            const Rcpp::IntegerVector fromi,
            const Rcpp::IntegerVector toi_in,
            const size_t nverts_in,
            const std::shared_ptr <DGraph> g_in,
            const std::string & heap_type_in,
            Rcpp::NumericMatrix dout_in) :
        dp_fromi (fromi), toi (toi_in), nverts (nverts_in),
        g (g_in), heap_type (heap_type_in), dout (dout_in)
    {
    }

    // Parallel function operator
    void operator() (std::size_t begin, std::size_t end)
    {
        std::shared_ptr<Dijkstra> dijkstra =
            std::make_shared <Dijkstra> (nverts,
                    *run_sp::getHeapImpl (heap_type), g);
        std::vector <double> w (nverts);
        std::vector <double> d (nverts);
        std::vector <int> prev (nverts);

        for (std::size_t i = begin; i < end; i++)
        {
            // These have to be reserved within the parallel operator function!
            std::fill (w.begin (), w.end (), INFINITE_DOUBLE);
            std::fill (d.begin (), d.end (), INFINITE_DOUBLE);

            if (heap_type.find ("set") == std::string::npos)
                dijkstra->run (d, w, prev,
                        static_cast <unsigned int> (dp_fromi [i]));
            else
                dijkstra->run_set (d, w, prev,
                        static_cast <unsigned int> (dp_fromi [i]));
            for (long int j = 0; j < toi.size (); j++)
            {
                if (w [static_cast <size_t> (toi [j])] < INFINITE_DOUBLE)
                {
                    dout (i, static_cast <size_t> (j)) =
                        d [static_cast <size_t> (toi [j])];
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

size_t run_sp::get_fromi_toi (const Rcpp::DataFrame &vert_map_in,
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
    return static_cast <size_t> (fromi.size ());
}

size_t run_sp::get_fromi (const Rcpp::DataFrame &vert_map_in,
        Rcpp::IntegerVector &fromi, Rcpp::NumericVector &id_vec)
{
    if (fromi [0] < 0) // use all vertices
    {
        id_vec = vert_map_in ["id"];
        fromi = id_vec;
    }
    return static_cast <size_t> (fromi.size ());
}

// Flows from the dijkstra output are reallocated based on matching vertex
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
        Rcpp::IntegerVector toi,
        const std::string& heap_type)
{
    Rcpp::NumericVector id_vec;
    size_t nfrom = run_sp::get_fromi_toi (vert_map_in, fromi, toi, id_vec);
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

    std::shared_ptr <DGraph> g = std::make_shared <DGraph> (nverts);
    inst_graph (g, nedges, vert_map, from, to, dist, wt);

    Rcpp::NumericVector na_vec = Rcpp::NumericVector (nfrom * nto,
            Rcpp::NumericVector::get_na ());
    Rcpp::NumericMatrix dout (static_cast <int> (nfrom),
            static_cast <int> (nto), na_vec.begin ());

    // Create parallel worker
    OneDist one_dist (fromi, toi, nverts, g, heap_type, dout);

    RcppParallel::parallelFor (0, static_cast <size_t> (fromi.length ()),
            one_dist);
    
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
    size_t nfrom = run_sp::get_fromi_toi (vert_map_in, fromi, toi, id_vec);
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

    std::shared_ptr <Dijkstra> dijkstra = std::make_shared <Dijkstra> (nverts,
            *run_sp::getHeapImpl(heap_type), g);

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
            //if (toi [vi] < INFINITE_INT)
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
        Rcpp::IntegerVector toi,
        const std::string& heap_type)
{
    Rcpp::NumericVector id_vec;
    size_t nfrom = run_sp::get_fromi_toi (vert_map_in, fromi, toi, id_vec);
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

    std::shared_ptr<Dijkstra> dijkstra = std::make_shared <Dijkstra> (nverts,
            *run_sp::getHeapImpl(heap_type), g);
    
    Rcpp::List res (nfrom);
    std::vector<double> w (nverts);
    std::vector<double> d (nverts);
    std::vector<int> prev (nverts);

    dijkstra->init (g); // specify the graph

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
            if (w [static_cast <size_t> (toi [vi])] < INFINITE_DOUBLE)
            {
                int target = toi [vi]; // target can be -1!
                while (target < INFINITE_INT)
                {
                    // Note that targets are all C++ 0-indexed and are converted
                    // directly here to R-style 1-indexes.
                    onePath.push_back (static_cast <unsigned int> (target + 1));
                    target = prev [static_cast <size_t> (target)];
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

struct OneFlow : public RcppParallel::Worker
{
    RcppParallel::RVector <int> dp_fromi;
    const Rcpp::IntegerVector toi;
    const Rcpp::NumericMatrix flows;
    const std::vector <std::string> vert_name;
    const std::unordered_map <std::string, unsigned int> verts_to_edge_map;
    size_t nverts; // can't be const because of reinterpret cast
    size_t nedges;
    const std::string dirtxt;
    const std::string heap_type;

    std::shared_ptr <DGraph> g;

    // constructor
    OneFlow (
            const Rcpp::IntegerVector fromi,
            const Rcpp::IntegerVector toi_in,
            const Rcpp::NumericMatrix flows_in,
            const std::vector <std::string>  vert_name_in,
            const std::unordered_map <std::string, unsigned int> verts_to_edge_map_in,
            const size_t nverts_in,
            const size_t nedges_in,
            const std::string dirtxt_in,
            const std::string &heap_type_in,
            const std::shared_ptr <DGraph> g_in) :
        dp_fromi (fromi), toi (toi_in), flows (flows_in), vert_name (vert_name_in),
        verts_to_edge_map (verts_to_edge_map_in),
        nverts (nverts_in), nedges (nedges_in), dirtxt (dirtxt_in),
        heap_type (heap_type_in), g (g_in)
    {
    }

    // Function to generate random file names
    std::string random_name(size_t len) {
        auto randchar = []() -> char
        {
            const char charset[] = \
               "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
            const size_t max_index = (sizeof(charset) - 1);
            //return charset [ rand() % max_index ];
            size_t i = static_cast <size_t> (floor (unif_rand () * max_index));
            return charset [i];
        };
        std::string str (len, 0);
        std::generate_n (str.begin(), len, randchar);
        return str;
    }

    // Parallel function operator
    void operator() (size_t begin, size_t end)
    {
        std::shared_ptr<Dijkstra> dijkstra =
            std::make_shared <Dijkstra> (nverts,
                    *run_sp::getHeapImpl (heap_type), g);
        std::vector <double> w (nverts);
        std::vector <double> d (nverts);
        std::vector <int> prev (nverts);

        std::vector <double> flowvec (nedges, 0.0);

        for (size_t i = begin; i < end; i++)
        {
            // These have to be reserved within the parallel operator function!
            std::fill (w.begin (), w.end (), INFINITE_DOUBLE);
            std::fill (d.begin (), d.end (), INFINITE_DOUBLE);

            dijkstra->run (d, w, prev,
                    static_cast <unsigned int> (dp_fromi [i]));
            for (size_t j = 0; j < static_cast <size_t> (toi.size ()); j++)
            {
                long int ltj = static_cast <long int> (j);
                if (dp_fromi [i] != toi [ltj]) // Exclude self-flows
                {
                    double flow_ij = flows (i, j);
                    if (w [static_cast <size_t> (toi [ltj])] < INFINITE_DOUBLE &&
                            flow_ij > 0.0)
                    {
                        int target = toi [ltj]; // can equal -1
                        while (target < INFINITE_INT)
                        {
                            size_t stt = static_cast <size_t> (target);
                            if (prev [stt] >= 0 && prev [stt] < INFINITE_INT)
                            {
                                std::string v2 = "f" +
                                    vert_name [static_cast <size_t> (prev [stt])] +
                                    "t" + vert_name [stt];
                                // multiple flows can aggregate to same edge, so
                                // this has to be +=, not just =!
                                flowvec [verts_to_edge_map.at (v2)] += flow_ij;
                            }

                            target = prev [stt];
                            // Only allocate that flow from origin vertex v to all
                            // previous vertices up until the target vi
                            if (target < 0 || target == dp_fromi [i])
                            {
                                break;
                            }
                        }
                    }
                }
            }
        } // end for i
        // dump flowvec to a file; chance of re-generating same file name is
        // 61^10, so there's no check for re-use of same
        std::string file_name = dirtxt + "flow_" + random_name (10) + ".dat";
        std::ofstream out_file;
        out_file.open (file_name, std::ios::binary | std::ios::out);
        out_file.write (reinterpret_cast <char *>(&nedges), sizeof (size_t));
        out_file.write (reinterpret_cast <char *>(&flowvec [0]),
                static_cast <std::streamsize> (nedges * sizeof (double)));
        out_file.close ();
    } // end parallel function operator
};

//' rcpp_aggregate_files
//'
//' @param file_names List of fill names of files (that is, with path) provided
//' from R, coz otherwise this is C++17 with an added library flag.
//' @param len Length of flows, which is simply the number of edges in the
//' graph.
//'
//' Each parallel flow aggregation worker dumps results to a randomly-named
//' file. This routine reassembles those results into a single aggregate vector.
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_aggregate_files (const Rcpp::CharacterVector file_names,
        const int len)
{
    Rcpp::NumericVector flows (len, 0.0);

    for (int i = 0; i < file_names.size (); i++)
    {
        size_t nedges;
        std::ifstream in_file (file_names [i], std::ios::binary | std::ios::in);
        in_file.read (reinterpret_cast <char *>(&nedges), sizeof (size_t));
        std::vector <double> flows_i (nedges);
        in_file.read (reinterpret_cast <char *>(&flows_i [0]),
                static_cast <std::streamsize> (nedges * sizeof (double)));
        in_file.close ();

        if (nedges != static_cast <size_t> (len))
            Rcpp::stop ("aggregate flows have inconsistent sizes");
        
        for (size_t j = 0; j < nedges; j++)
            flows [static_cast <long> (j)] += flows_i [j];
    }
    return flows;
}

//' rcpp_flows_aggregate_par
//'
//' @param graph The data.frame holding the graph edges
//' @param vert_map_in map from <std::string> vertex ID to (0-indexed) integer
//' index of vertices
//' @param fromi Index into vert_map_in of vertex numbers
//' @param toi Index into vert_map_in of vertex numbers
//'
//' @note The parallelisation is achieved by dumping the results of each thread
//' to a file, with aggregation performed at the end by simply reading back and
//' aggregating all files. There is no way to aggregate into a single vector
//' because threads have to be independent. The only danger with this approach
//' is that multiple threads may generate the same file names, but with names 10
//' characters long, that chance should be 1 / 62 ^ 10.
//'
//' @noRd
// [[Rcpp::export]]
void rcpp_flows_aggregate_par (const Rcpp::DataFrame graph,
        const Rcpp::DataFrame vert_map_in,
        Rcpp::IntegerVector fromi,
        Rcpp::IntegerVector toi,
        const Rcpp::NumericMatrix flows,
        const std::string dirtxt,
        const std::string heap_type)
{
    Rcpp::NumericVector id_vec;
    const size_t nfrom = run_sp::get_fromi_toi (vert_map_in, fromi, toi, id_vec);

    const std::vector <std::string> from = graph ["from"];
    const std::vector <std::string> to = graph ["to"];
    const std::vector <double> dist = graph ["d"];
    const std::vector <double> wt = graph ["w"];

    const unsigned int nedges = static_cast <unsigned int> (graph.nrow ());
    const std::vector <std::string> vert_name = vert_map_in ["vert"];
    const std::vector <unsigned int> vert_indx = vert_map_in ["id"];
    // Make map from vertex name to integer index
    std::map <std::string, unsigned int> vert_map_i;
    const size_t nverts = run_sp::make_vert_map (vert_map_in, vert_name,
            vert_indx, vert_map_i);

    std::unordered_map <std::string, unsigned int> verts_to_edge_map;
    std::unordered_map <std::string, double> verts_to_dist_map;
    run_sp::make_vert_to_edge_maps (from, to, wt, verts_to_edge_map, verts_to_dist_map);

    std::shared_ptr <DGraph> g = std::make_shared <DGraph> (nverts);
    inst_graph (g, nedges, vert_map_i, from, to, dist, wt);

    // Create parallel worker
    OneFlow one_flow (fromi, toi, flows, vert_name, verts_to_edge_map,
            nverts, nedges, dirtxt, heap_type, g);

    GetRNGstate (); // Initialise R random seed
    RcppParallel::parallelFor (0, nfrom, one_flow);
    PutRNGstate ();
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
    size_t nfrom = run_sp::get_fromi (vert_map_in, fromi, id_vec);

    std::vector <std::string> from = graph ["from"];
    std::vector <std::string> to = graph ["to"];
    std::vector <double> dist = graph ["d"];
    std::vector <double> wt = graph ["w"];

    unsigned int nedges = static_cast <unsigned int> (graph.nrow ());
    std::vector <std::string> vert_name = vert_map_in ["vert"];
    std::vector <unsigned int> vert_indx = vert_map_in ["id"];
    // Make map from vertex name to integer index
    std::map <std::string, unsigned int> vert_map_i;
    size_t nverts = run_sp::make_vert_map (vert_map_in, vert_name,
            vert_indx, vert_map_i);

    std::unordered_map <std::string, unsigned int> verts_to_edge_map;
    std::unordered_map <std::string, double> verts_to_dist_map;
    run_sp::make_vert_to_edge_maps (from, to, wt, verts_to_edge_map, verts_to_dist_map);

    std::shared_ptr<DGraph> g = std::make_shared<DGraph> (nverts);
    inst_graph (g, nedges, vert_map_i, from, to, dist, wt);

    std::shared_ptr <Dijkstra> dijkstra = std::make_shared <Dijkstra> (nverts,
            *run_sp::getHeapImpl(heap_type), g);

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
                    vert_from = vert_name [static_cast <size_t> (prev [vi])];
                std::string two_verts = "f" + vert_from + "t" + vert_to;
                if (verts_to_edge_map.find (two_verts) == verts_to_edge_map.end ())
                    Rcpp::stop ("vertex pair forms no known edge");

                unsigned int indx = verts_to_edge_map [two_verts];
                if (d [vi] < INFINITE_DOUBLE)
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
