
#include "run_sp.h"
#include "flows.h"

#include "dgraph.h"
#include "heaps/heap_lib.h"

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

struct OneAggregate : public RcppParallel::Worker
{
    RcppParallel::RVector <int> dp_fromi;
    const std::vector <size_t> toi;
    const RcppParallel::RMatrix <double> flows;
    const std::vector <std::string> vert_name;
    const std::unordered_map <std::string, size_t> verts_to_edge_map;
    size_t nverts; // can't be const because of reinterpret cast
    size_t nedges;
    const bool norm_sums;
    const double tol;
    const std::string heap_type;
    std::shared_ptr <DGraph> g;

    std::vector <double> output;

    // Constructor 1: The main constructor
    OneAggregate (
            const RcppParallel::RVector <int> fromi,
            const std::vector <size_t> toi_in,
            const RcppParallel::RMatrix <double> flows_in,
            const std::vector <std::string>  vert_name_in,
            const std::unordered_map <std::string, size_t> verts_to_edge_map_in,
            const size_t nverts_in,
            const size_t nedges_in,
            const bool norm_sums_in,
            const double tol_in,
            const std::string &heap_type_in,
            const std::shared_ptr <DGraph> g_in) :
        dp_fromi (fromi), toi (toi_in), flows (flows_in),
        vert_name (vert_name_in),
        verts_to_edge_map (verts_to_edge_map_in),
        nverts (nverts_in), nedges (nedges_in), norm_sums (norm_sums_in),
        tol (tol_in), heap_type (heap_type_in), g (g_in), output ()
    {
        output.resize (nedges, 0.0);
    }

    // Constructor 2: The Split constructor
    OneAggregate (
            const OneAggregate& oneAggregate,
            RcppParallel::Split) :
        dp_fromi (oneAggregate.dp_fromi),
        toi (oneAggregate.toi),
        flows (oneAggregate.flows),
        vert_name (oneAggregate.vert_name),
        verts_to_edge_map (oneAggregate.verts_to_edge_map),
        nverts (oneAggregate.nverts),
        nedges (oneAggregate.nedges),
        norm_sums (oneAggregate.norm_sums),
        tol (oneAggregate.tol),
        heap_type (oneAggregate.heap_type),
        g (oneAggregate.g), output ()
    {
        output.resize (nedges, 0.0);
    }

    // Parallel function operator
    void operator() (size_t begin, size_t end)
    {
        std::shared_ptr<PF::PathFinder> pathfinder =
            std::make_shared <PF::PathFinder> (nverts,
                    *run_sp::getHeapImpl (heap_type), g);
        std::vector <double> w (nverts);
        std::vector <double> d (nverts);
        std::vector <long int> prev (nverts);

        for (size_t i = begin; i < end; i++)
        {
            //if (RcppThread::isInterrupted (i % static_cast<int>(100) == 0))
            //if (RcppThread::isInterrupted ())
            //    return;

            // These have to be reserved within the parallel operator function!
            std::fill (w.begin (), w.end (), INFINITE_DOUBLE);
            std::fill (d.begin (), d.end (), INFINITE_DOUBLE);
            std::fill (prev.begin (), prev.end (), INFINITE_INT);

            size_t from_i = static_cast <size_t> (dp_fromi [i]);
            d [from_i] = w [from_i] = 0.0;

            // reduce toi to only those within tolerance limt
            double flim = 0.0;
            size_t nto = toi.size ();
            if (tol > 0.0) {
                double fmax = 0.0;
                for (size_t j = 0; j < static_cast <size_t> (flows.ncol ()); j++)
                    if (flows (i, j) > fmax)
                        fmax = flows (i, j);
                flim = fmax * tol;
                size_t nto = 0;
                for (size_t j = 0; j < static_cast <size_t> (flows.ncol ()); j++)
                    if (flows (i, j) > flim)
                        nto++;

                if (nto == 0)
                    continue;
            }

            // toi_index is into cols of flow matrix
            // toi_reduced is into the vertex vectors (d, w, prev)
            std::vector <size_t> toi_reduced, toi_index;
            toi_reduced.reserve (nto);
            toi_index.reserve (nto);
            for (size_t j = 0; j < static_cast <size_t> (flows.ncol ()); j++)
                if (flows (i, j) > flim)
                {
                    toi_index.push_back (static_cast <size_t> (j));
                    toi_reduced.push_back (toi [j]);
                }

            pathfinder->Dijkstra (d, w, prev, from_i, toi_reduced);
            for (size_t j = 0; j < toi_reduced.size (); j++)
            {
                if (from_i != toi_reduced [j]) // Exclude self-flows
                {
                    double flow_ij = flows (i, toi_index [j]);
                    if (w [toi_reduced [j]] < INFINITE_DOUBLE && flow_ij > 0.0)
                    {
                        // count how long the path is, so flows on
                        // each edge can be divided by this length
                        int path_len = 1;
                        if (norm_sums)
                        {
                            path_len = 0;
                            long int target_t = static_cast <long int> (toi_reduced [j]);
                            size_t from_t = static_cast <size_t> (dp_fromi [i]);
                            while (target_t < INFINITE_INT)
                            {
                                path_len++;
                                size_t target_size_t = static_cast <size_t> (target_t);
                                target_t = prev [target_size_t];
                                if (target_t < 0 || target_size_t == from_t)
                                    break;
                            }
                        }

                        long int target = static_cast <int> (toi_reduced [j]); // can equal -1
                        while (target < INFINITE_INT)
                        {
                            size_t stt = static_cast <size_t> (target);
                            if (prev [stt] >= 0 && prev [stt] < INFINITE_INT)
                            {
                                std::string v2 = "f" +
                                    vert_name [static_cast <size_t> (prev [stt])] +
                                    "t" + vert_name [stt];
                                output [verts_to_edge_map.at (v2)] +=
                                    flow_ij / static_cast <double> (path_len);
                            }

                            target = prev [stt];
                            // Only allocate that flow from origin vertex v to all
                            // previous vertices up until the target vi
                            if (target < 0L || target == dp_fromi [i])
                            {
                                break;
                            }
                        }
                    }
                }
            }
        } // end for i
    } // end parallel function operator

    void join (const OneAggregate &rhs)
    {
        for (size_t i = 0; i < output.size (); i++)
            output [i] += rhs.output [i];
    }
};


struct OneDisperse : public RcppParallel::Worker
{
    RcppParallel::RVector <int> dp_fromi;
    const RcppParallel::RVector <double> dens;
    const std::vector <std::string> vert_name;
    const std::unordered_map <std::string, size_t> verts_to_edge_map;
    size_t nverts; // can't be const because of reinterpret cast
    size_t nedges;
    const RcppParallel::RVector <double> kfrom;
    const double tol;
    const std::string heap_type;
    std::shared_ptr <DGraph> g;

    std::vector <double> output;

    // Constructor 1: The main constructor
    OneDisperse (
            const RcppParallel::RVector <int> fromi,
            const RcppParallel::RVector <double> dens_in,
            const std::vector <std::string>  vert_name_in,
            const std::unordered_map <std::string, size_t> verts_to_edge_map_in,
            const size_t nverts_in,
            const size_t nedges_in,
            const RcppParallel::RVector <double> kfrom_in,
            const double tol_in,
            const std::string &heap_type_in,
            const std::shared_ptr <DGraph> g_in) :
        dp_fromi (fromi), dens (dens_in), vert_name (vert_name_in),
        verts_to_edge_map (verts_to_edge_map_in),
        nverts (nverts_in), nedges (nedges_in), kfrom (kfrom_in),
        tol (tol_in), heap_type (heap_type_in), g (g_in), output ()
    {
        const R_xlen_t nfrom = static_cast <R_xlen_t> (dens.size ());
        const R_xlen_t nk = static_cast <R_xlen_t> (kfrom.size ()) / nfrom;
        const size_t out_size = nedges * static_cast <size_t> (nk);
        output.resize (out_size, 0.0);
    }

    // Constructor 2: The Split constructor
    OneDisperse (
            const OneDisperse& oneDisperse,
            RcppParallel::Split) :
        dp_fromi (oneDisperse.dp_fromi), dens (oneDisperse.dens),
        vert_name (oneDisperse.vert_name),
        verts_to_edge_map (oneDisperse.verts_to_edge_map),
        nverts (oneDisperse.nverts), nedges (oneDisperse.nedges),
        kfrom (oneDisperse.kfrom), tol (oneDisperse.tol),
        heap_type (oneDisperse.heap_type), g (oneDisperse.g), output ()
    {
        const R_xlen_t nfrom = static_cast <R_xlen_t> (dens.size ());
        const R_xlen_t nk = static_cast <R_xlen_t> (kfrom.size ()) / nfrom;
        size_t out_size = nedges * static_cast <size_t> (nk);
        output.resize (out_size, 0.0);
    }


    // Parallel function operator
    void operator() (size_t begin, size_t end)
    {
        std::shared_ptr<PF::PathFinder> pathfinder =
            std::make_shared <PF::PathFinder> (nverts,
                    *run_sp::getHeapImpl (heap_type), g);
        std::vector <double> w (nverts);
        std::vector <double> d (nverts);
        std::vector <long int> prev (nverts);

        const size_t nfrom = dens.size ();
        const size_t nk = kfrom.size () / nfrom;

        for (size_t i = begin; i < end; i++) // over the from vertices
        {
            //if (RcppThread::isInterrupted (i % static_cast<int>(10) == 0))
            if (RcppThread::isInterrupted ())
                return;

            // translate k-value to distance limit based on tol
            // exp(-d / k) = tol -> d = -k * log (tol)
            // k_from holds nk vectors of different k-values, each of length
            // nedges.
            double dlim = 0.0;
            for (size_t k = 0; k < nk; k++)
                if (kfrom [i + k * nfrom] > dlim)
                    dlim = kfrom [i + k * nfrom]; // dlim is max k-value
            dlim = -dlim * log (tol); // converted to actual dist limit.

            std::fill (w.begin (), w.end (), INFINITE_DOUBLE);
            std::fill (d.begin (), d.end (), INFINITE_DOUBLE);
            std::fill (prev.begin (), prev.end (), INFINITE_INT);

            const size_t from_i = static_cast <size_t> (dp_fromi [i]);
            d [from_i] = w [from_i] = 0.0;

            pathfinder->DijkstraLimit (d, w, prev, from_i, dlim);

            std::vector <double> flows_i (nedges * nk, 0.0);
            std::vector <double> expsum (nk, 0.0);

            for (size_t j = 0; j < nverts; j++)
            {
                if (prev [j] >= 0 && prev [j] < INFINITE_INT)
                {
                    const std::string vert_to = vert_name [j],
                        vert_from = vert_name [static_cast <size_t> (prev [j])];
                    const std::string two_verts = "f" + vert_from + "t" + vert_to;

                    size_t index = verts_to_edge_map.at (two_verts);
                    if (d [j] < INFINITE_DOUBLE)
                    {
                        for (size_t k = 0; k < nk; k++)
                        {
                            double exp_jk;
                            if (kfrom [i + k * nfrom] > 0.0)
                                exp_jk = exp (-d [j] / kfrom [i + k * nfrom]);
                            else
                            {
                                // standard logistic polygonal for UK cycling
                                // models
                                double lp = -3.894 + (-0.5872 * d [j]) +
                                    (1.832 * sqrt (d [j])) +
                                    (0.007956 * d [j] * d [j]);
                                exp_jk = exp (lp) / (1.0 + exp (lp));
                            }
                            const size_t k_st = static_cast <size_t> (k);
                            expsum [k_st] += exp_jk;
                            flows_i [index + k_st * nedges] += dens [i] * exp_jk;
                        }
                    }
                }
            } // end for j
            for (size_t k = 0; k < nk; k++)
                if (expsum [k] > tol)
                    for (size_t j = 0; j < nedges; j++)
                        output [j + k * nedges] +=
                            flows_i [j + k * nedges] / expsum [k];
        } // end for i
    } // end parallel function operator

    void join (const OneDisperse &rhs)
    {
        for (size_t i = 0; i < output.size (); i++)
            output [i] += rhs.output [i];
    }
};

struct OneSI : public RcppParallel::Worker
{
    RcppParallel::RVector <int> dp_fromi;
    const std::vector <size_t> toi;
    const RcppParallel::RVector <double> k_from;
    const RcppParallel::RVector <double> dens_from;
    const RcppParallel::RVector <double> dens_to;
    const std::vector <std::string> vert_name;
    const std::unordered_map <std::string, size_t> verts_to_edge_map;
    size_t nverts; // can't be const because of reinterpret cast
    size_t nedges;
    const bool norm_sums;
    const double tol;
    const std::string heap_type;
    std::shared_ptr <DGraph> g;

    std::vector <double> output;

    // Constructor 1: The main constructor
    OneSI (
            const RcppParallel::RVector <int> fromi,
            const std::vector <size_t> toi_in,
            const RcppParallel::RVector <double> k_from_in,
            const RcppParallel::RVector <double> dens_from_in,
            const RcppParallel::RVector <double> dens_to_in,
            const std::vector <std::string>  vert_name_in,
            const std::unordered_map <std::string, size_t> verts_to_edge_map_in,
            const size_t nverts_in,
            const size_t nedges_in,
            const bool norm_sums_in,
            const double tol_in,
            const std::string &heap_type_in,
            const std::shared_ptr <DGraph> g_in) :
        dp_fromi (fromi), toi (toi_in), k_from (k_from_in),
        dens_from (dens_from_in), dens_to (dens_to_in),
        vert_name (vert_name_in), verts_to_edge_map (verts_to_edge_map_in),
        nverts (nverts_in), nedges (nedges_in), norm_sums (norm_sums_in),
        tol (tol_in), heap_type (heap_type_in), g (g_in), output ()
    {
        const size_t nfrom = dens_from.size ();
        const size_t nk = k_from.size () / nfrom;
        const size_t out_size = nedges * nk;
        output.resize (out_size, 0.0);
    }

    // Constructor 2: The Split constructor
    OneSI (
            const OneSI& oneSI,
            RcppParallel::Split) :
        dp_fromi (oneSI.dp_fromi), toi (oneSI.toi), k_from (oneSI.k_from),
        dens_from (oneSI.dens_from), dens_to (oneSI.dens_to),
        vert_name (oneSI.vert_name), verts_to_edge_map (oneSI.verts_to_edge_map),
        nverts (oneSI.nverts), nedges (oneSI.nedges),
        norm_sums (oneSI.norm_sums), tol (oneSI.tol),
        heap_type (oneSI.heap_type), g (oneSI.g), output ()
    {
        const size_t nfrom = dens_from.size ();
        const size_t nk = k_from.size () / nfrom;
        const size_t out_size = nedges * nk;
        output.resize (out_size, 0.0);
    }

    // Parallel function operator
    void operator() (size_t begin, size_t end)
    {
        std::shared_ptr<PF::PathFinder> pathfinder =
            std::make_shared <PF::PathFinder> (nverts,
                    *run_sp::getHeapImpl (heap_type), g);
        std::vector <double> w (nverts);
        std::vector <double> d (nverts);
        std::vector <long int> prev (nverts);

        // k_from can have multiple vectors of k-values, each equal in length to
        // the number of 'from' points. The output is then a single vector of
        // 'nedges' wrapped 'nk' times.
        const size_t nfrom = dens_from.size ();
        const size_t nk = k_from.size () / nfrom;

        for (size_t i = begin; i < end; i++)
        {
            //if (RcppThread::isInterrupted (i % static_cast<int>(100) == 0))
            if (RcppThread::isInterrupted ())
                return;

            // These have to be reserved within the parallel operator function!
            std::fill (w.begin (), w.end (), INFINITE_DOUBLE);
            std::fill (d.begin (), d.end (), INFINITE_DOUBLE);
            std::fill (prev.begin (), prev.end (), INFINITE_INT);

            size_t from_i = static_cast <size_t> (dp_fromi [i]);
            d [from_i] = w [from_i] = 0.0;

            // translate k-value to distance limit based on tol
            // exp(-d / k) = tol -> d = -k * log (tol)
            // const double dlim = -k_from [i] * log (tol);
            // k_from holds nk vectors of different k-values, each of length
            // nedges.
            double dlim = 0.0;
            for (size_t k = 0; k < nk; k++)
                if (k_from [i + k * nfrom] > dlim)
                    dlim = k_from [i + k * nfrom];
            dlim = -dlim * log (tol);

            pathfinder->DijkstraLimit (d, w, prev, from_i, dlim);

            std::vector <double> flows_i (nedges * nk, 0.0);
            std::vector <double> expsum (nk, 0.0);
            for (size_t j = 0; j < toi.size (); j++)
            {
                if (from_i != toi [j]) // Exclude self-flows
                {
                    if (d [toi [j]] < INFINITE_DOUBLE)
                    {
                        /* Flow between i and j is based on d (i, j), but is
                         * subsequently allocated to all intermediate edges.
                         * The `norm_sums` option re-scales this intermediate
                         * edge allocation such that SI value to one terminal
                         * vertex from the origin vertex of:
                         * N_i N_j exp (-d_{ij}) / sum_j (N_{ij} exp (-d_{ij}))
                         * is divided by the number of intermediate edges
                         * between i and j. This is achieved by initially only
                         * allocating the values of the numerator divided by
                         * numbers of intermediate edges, while accumulating the
                         * denominator independent of that value in order to
                         * divide all final values.
                         */
                        std::vector <double> exp_jk (nk, 0.0),
                            flow_ijk (nk, 0.0);
                        for (size_t k = 0; k < nk; k++)
                        {
                            exp_jk [k] = dens_to [j] *
                                exp (-d [toi [j]] / k_from [i + k * nfrom]);
                            flow_ijk [k] = dens_from [i] * exp_jk [k];
                            expsum [k] += exp_jk [k];
                        }

                        // need to count how long the path is, so flows on
                        // each edge can be divided by this length
                        int path_len = 1;
                        if (norm_sums)
                        {
                            path_len = 0;
                            long int target_t = static_cast <long int> (toi [j]);
                            size_t from_t = static_cast <size_t> (dp_fromi [i]);
                            while (target_t < INFINITE_INT)
                            {
                                size_t target_size_t = static_cast <size_t> (target_t);
                                path_len++;
                                target_t = prev [target_size_t];
                                if (target_t < 0 || target_size_t == from_t)
                                    break;
                            }
                        } 

                        long int target = static_cast <long int> (toi [j]); // can equal -1
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
                                size_t index = verts_to_edge_map.at (v2);
                                for (size_t k = 0; k < nk; k++)
                                {
                                    flows_i [index + k * nedges] +=
                                        flow_ijk [k] / 
                                        static_cast <double> (path_len);
                                }
                            }

                            target = static_cast <long int> (prev [stt]);
                            // Only allocate that flow from origin vertex v to all
                            // previous vertices up until the target vi
                            if (target < 0L || target == dp_fromi [i])
                            {
                                break;
                            }
                        }
                    }
                }
            } // end for j
            for (size_t k = 0; k < nk; k++)
                if (expsum [k] > tol)
                    for (size_t j = 0; j < nedges; j++)
                        output [j + k * nedges] +=
                            flows_i [j + k * nedges] / expsum [k];
        } // end for i
    } // end parallel function operator

    void join (const OneSI& rhs)
    {
        for (size_t i = 0; i < output.size (); i++)
            output [i] += rhs.output [i];
    }
};

//' rcpp_flows_aggregate_par
//'
//' @param graph The data.frame holding the graph edges
//' @param vert_map_in map from <std::string> vertex ID to (0-indexed) integer
//' index of vertices
//' @param fromi Index into vert_map_in of vertex numbers
//' @param toi Index into vert_map_in of vertex numbers
//' @param tol Relative tolerance in terms of flows below which targets
//' (to-vertices) are not considered.
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
Rcpp::NumericVector rcpp_flows_aggregate_par (const Rcpp::DataFrame graph,
        const Rcpp::DataFrame vert_map_in,
        Rcpp::IntegerVector fromi,
        Rcpp::IntegerVector toi_in,
        Rcpp::NumericMatrix flows,
        const bool norm_sums,
        const double tol,
        const std::string heap_type)
{
    std::vector <size_t> toi =
        Rcpp::as <std::vector <size_t> > (toi_in);
    const size_t nfrom = static_cast <size_t> (fromi.size ());

    const std::vector <std::string> from = graph ["from"];
    const std::vector <std::string> to = graph ["to"];
    const std::vector <double> dist = graph ["d"];
    const std::vector <double> wt = graph ["d_weighted"];

    const size_t nedges = static_cast <size_t> (graph.nrow ());
    const std::vector <std::string> vert_name = vert_map_in ["vert"];
    const std::vector <size_t> vert_indx = vert_map_in ["id"];
    // Make map from vertex name to integer index
    std::map <std::string, size_t> vert_map_i;
    const size_t nverts = run_sp::make_vert_map (vert_map_in, vert_name,
            vert_indx, vert_map_i);

    std::unordered_map <std::string, size_t> verts_to_edge_map;
    std::unordered_map <std::string, double> verts_to_dist_map;
    run_sp::make_vert_to_edge_maps (from, to, wt, verts_to_edge_map, verts_to_dist_map);

    std::shared_ptr <DGraph> g = std::make_shared <DGraph> (nverts);
    inst_graph (g, nedges, vert_map_i, from, to, dist, wt);

    // Create parallel worker
    OneAggregate oneAggregate (RcppParallel::RVector <int> (fromi), toi,
            RcppParallel::RMatrix <double> (flows), vert_name, verts_to_edge_map,
            nverts, nedges, norm_sums, tol, heap_type, g);

    size_t chunk_size = run_sp::get_chunk_size (nfrom);
    RcppParallel::parallelReduce (0, nfrom, oneAggregate, chunk_size);

    return Rcpp::wrap (oneAggregate.output);
}



//' rcpp_flows_disperse_par
//'
//' Modified version of \code{rcpp_flows_aggregate} that aggregates flows to all
//' destinations from given set of origins, with flows attenuated by distance
//' from those origins.
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
//' between each pair of from and to points.
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_flows_disperse_par (const Rcpp::DataFrame graph,
        const Rcpp::DataFrame vert_map_in,
        Rcpp::IntegerVector fromi,
        Rcpp::NumericVector k,
        Rcpp::NumericVector dens,
        const double &tol,
        std::string heap_type)
{
    const size_t nfrom = static_cast <size_t> (fromi.size ());

    std::vector <std::string> from = graph ["from"];
    std::vector <std::string> to = graph ["to"];
    std::vector <double> dist = graph ["d"];
    std::vector <double> wt = graph ["d_weighted"];

    size_t nedges = static_cast <size_t> (graph.nrow ());
    std::vector <std::string> vert_name = vert_map_in ["vert"];
    std::vector <size_t> vert_indx = vert_map_in ["id"];
    // Make map from vertex name to integer index
    std::map <std::string, size_t> vert_map_i;
    size_t nverts = run_sp::make_vert_map (vert_map_in, vert_name,
            vert_indx, vert_map_i);

    std::unordered_map <std::string, size_t> verts_to_edge_map;
    std::unordered_map <std::string, double> verts_to_dist_map;
    run_sp::make_vert_to_edge_maps (from, to, wt, verts_to_edge_map, verts_to_dist_map);

    std::shared_ptr<DGraph> g = std::make_shared<DGraph> (nverts);
    inst_graph (g, nedges, vert_map_i, from, to, dist, wt);

    // Create parallel worker
    OneDisperse oneDisperse (RcppParallel::RVector <int> (fromi),
            RcppParallel::RVector <double> (dens),
            vert_name, verts_to_edge_map,
            nverts, nedges,
            RcppParallel::RVector <double> (k), tol, heap_type, g);

    size_t chunk_size = run_sp::get_chunk_size (nfrom);
    RcppParallel::parallelReduce (0, nfrom, oneDisperse, chunk_size);

    return Rcpp::wrap (oneDisperse.output);
}


//' rcpp_flows_si
//'
//' @param graph The data.frame holding the graph edges
//' @param vert_map_in map from <std::string> vertex ID to (0-indexed) integer
//' index of vertices
//' @param fromi Index into vert_map_in of vertex numbers
//' @param toi Index into vert_map_in of vertex numbers
//' @param kvec Vector of k-values for each fromi
//' @param nvec Vector of density-values for each fromi
//' @param tol Relative tolerance in terms of flows below which targets
//' (to-vertices) are not considered.
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_flows_si (const Rcpp::DataFrame graph,
        const Rcpp::DataFrame vert_map_in,
        Rcpp::IntegerVector fromi,
        Rcpp::IntegerVector toi_in,
        Rcpp::NumericVector kvec,
        Rcpp::NumericVector dens_from,
        Rcpp::NumericVector dens_to,
        const bool norm_sums,
        const double tol,
        const std::string heap_type)
{
    std::vector <size_t> toi =
        Rcpp::as <std::vector <size_t> > ( toi_in);
    const size_t nfrom = static_cast <size_t> (fromi.size ());

    const std::vector <std::string> from = graph ["from"];
    const std::vector <std::string> to = graph ["to"];
    const std::vector <double> dist = graph ["d"];
    const std::vector <double> wt = graph ["d_weighted"];

    const size_t nedges = static_cast <size_t> (graph.nrow ());
    const std::vector <std::string> vert_name = vert_map_in ["vert"];
    const std::vector <size_t> vert_indx = vert_map_in ["id"];
    // Make map from vertex name to integer index
    std::map <std::string, size_t> vert_map_i;
    const size_t nverts = run_sp::make_vert_map (vert_map_in, vert_name,
            vert_indx, vert_map_i);

    std::unordered_map <std::string, size_t> verts_to_edge_map;
    std::unordered_map <std::string, double> verts_to_dist_map;
    run_sp::make_vert_to_edge_maps (from, to, wt, verts_to_edge_map, verts_to_dist_map);

    std::shared_ptr <DGraph> g = std::make_shared <DGraph> (nverts);
    inst_graph (g, nedges, vert_map_i, from, to, dist, wt);

    // Create parallel worker
    OneSI oneSI (RcppParallel::RVector <int> (fromi), toi, 
            RcppParallel::RVector <double> (kvec),
            RcppParallel::RVector <double> (dens_from),
            RcppParallel::RVector <double> (dens_to),
            vert_name, verts_to_edge_map,
            nverts, nedges, norm_sums, tol, heap_type, g);

    size_t chunk_size = run_sp::get_chunk_size (nfrom);
    RcppParallel::parallelReduce (0, nfrom, oneSI, chunk_size);

    return Rcpp::wrap (oneSI.output);
}
