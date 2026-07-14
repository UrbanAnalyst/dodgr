#pragma once

// Modified from code by Shane Saunders

#include <vector>
#include <memory>
#include <set>
#include <cmath> // fabs
#include <unordered_set> // used in bidirected search

// inf_dbl used only in AStar2:
#include <limits>

#include "dgraph.h"

constexpr long int INFINITE_INT =  std::numeric_limits <long int>::max ();
constexpr double INFINITE_DOUBLE =  std::numeric_limits <double>::max ();

class Heap;
class HeapDesc;
class DGraph;

namespace PF{

// Lazily-evaluated Euclidean A* heuristic: straight-line distances from the
// search origin are computed on demand as each vertex is first encountered,
// rather than pre-computed for every vertex in the graph up front. A* with
// early termination on reaching all targets often visits only a small
// fraction of vertices, so eager pre-computation of the full array is wasted
// work on large graphs.
struct AStarHeuristic
{
    const std::vector <double> &vx;
    const std::vector <double> &vy;
    const double x0;
    const double y0;

    AStarHeuristic (const std::vector <double> &vx_in,
            const std::vector <double> &vy_in,
            const size_t v0) :
        vx (vx_in), vy (vy_in), x0 (vx_in [v0]), y0 (vy_in [v0])
    {
    }

    double operator[] (const size_t v) const
    {
        const double dx = vx [v] - x0, dy = vy [v] - y0;
        return sqrt (dx * dx + dy * dy);
    }
};

/* DijkstraEdge is used for the std::set implementation, everything else is for
 * Shane Saunders's heap sort versions */
struct DijkstraEdge
{
    DijkstraEdge (double wt, size_t i): _wt (wt), _i (i) {}
    
    double _wt;
    size_t _i;

    size_t geti () const { return _i;   }
    double getw () const { return _wt;   }
};

struct by_wt
{
    bool operator () (const DijkstraEdge& lhs, const DijkstraEdge& rhs) const
    {
        if (fabs (lhs._wt - rhs._wt) < 1.0e-12)
            return lhs._i < rhs._i;
        else
            return lhs._wt < rhs._wt;
    }
}; 

typedef std::set <DijkstraEdge, by_wt> EdgeSet;


class PathFinder {
    public:
        PathFinder (size_t n,
                const HeapDesc& heapD,
                std::shared_ptr<const DGraph> g);
        ~PathFinder();
        
        // Disable copy/assign as will crash
        PathFinder(const PathFinder&) = delete;
        PathFinder& operator=(const PathFinder&) = delete;

        void init (std::shared_ptr<const DGraph> g);
        void init_arrays (
                std::vector<double>& d,
                std::vector<double>& w,
                std::vector<long int>& prev,
                bool *m_open_vec,
                bool *m_closed_vec,
                const size_t v,
                const size_t n);
        void scan_edges (
                const DGraphEdge *edge,
                std::vector<double>& d,
                std::vector<double>& w,
                std::vector<long int>& prev,
                bool *m_open_vec,
                const bool *m_closed_vec,
                const size_t &v0);
        void scan_edges_heur ( // with A* heuristic
                const DGraphEdge *edge,
                std::vector<double>& d,
                std::vector<double>& w,
                std::vector<long int>& prev,
                bool *m_open_vec,
                const bool *m_closed_vec,
                const size_t &v0,
                const AStarHeuristic &heur);
        // with A* heuristic for dists-categorical
        void scan_edge_types_heur (
                const DGraphEdge *edge,
                std::vector<double>& d,
                std::vector<double>& w,
                std::vector<long int>& prev,
                bool *m_open_vec,
                const bool *m_closed_vec,
                const size_t &v0,
                const AStarHeuristic &heur);
        // run_sp_categorical for threshold dists
        void scan_edge_types (
                const DGraphEdge *edge,
                std::vector<double>& d,
                std::vector<double>& w,
                std::vector<long int>& prev,
                bool *m_open_vec,
                const bool *m_closed_vec,
                const size_t &v0);

        void Dijkstra (
                std::vector<double>& d,
                std::vector<double>& w,
                std::vector<long int>& prev,
                const size_t v0,
                const std::vector <size_t> &to_index);
        void DijkstraNearest (
                std::vector<double>& d,
                std::vector<double>& w,
                std::vector<long int>& prev,
                const size_t v0,
                const std::vector <size_t> &to_index);
        void DijkstraLimit ( // run_sp_categorical
                std::vector<double>& d,
                std::vector<double>& w,
                std::vector<long int>& prev,
                const size_t v0,
                const double &dlim);
        void DijkstraLimitEdgeType (
                std::vector<double>& d,
                std::vector<double>& w,
                std::vector<long int>& prev,
                const size_t v0,
                const double &dlim);
        void AStar (std::vector<double>& d,
                std::vector<double>& w,
                std::vector<long int>& prev,
                const AStarHeuristic& heur,
                const size_t v0,
                const std::vector <size_t> &to_index);
        void AStarEdgeType (std::vector<double>& d,
                std::vector<double>& w,
                std::vector<long int>& prev,
                const AStarHeuristic& heur,
                const size_t v0,
                const std::vector <size_t> &to_index);
        void Dijkstra_set (std::vector <double>& d,
                std::vector<double>& w,
                std::vector<long int>& prev,
                size_t v0);
        void Centrality_vertex (
                std::vector <double>& cent,
                const size_t s,
                const double vert_wt,
                const double dist_threshold);
        void Centrality_edge (
                std::vector <double>& cent,
                const size_t s,
                const double vert_wt,
                const size_t nedges,
                const double dist_threshold);

    private:
        Heap *m_heap;        // pointer: heap
        // Convert to vector<bool>? (save memory, might be a performance hit though)
        bool *m_open;           // array: frontier set state of vertices
        bool *m_closed;         // also for bi-dir

        std::shared_ptr<const DGraph> m_graph;    // pointer: directed graph    

        EdgeSet edge_set;
};

} // end namespace PF
