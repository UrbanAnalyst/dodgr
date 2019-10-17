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

constexpr double INFINITE_INT =  std::numeric_limits<int>::max ();
constexpr double INFINITE_DOUBLE =  std::numeric_limits<double>::max ();

class Heap;
class HeapDesc;
class DGraph;

namespace PF{

/* DijkstraEdge is used for the std::set implementation, everything else is for
 * Shane Saunders's heap sort versions */
struct DijkstraEdge
{
    DijkstraEdge (double wt, unsigned int i): _wt (wt), _i (i) {}
    
    double _wt;
    unsigned int _i;

    unsigned int geti () const { return _i;   }
    double getw () const { return _wt;   }
};

struct by_wt
{
    bool operator () (const DijkstraEdge& lhs, const DijkstraEdge& rhs)
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
        PathFinder (unsigned int n,
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
                std::vector<int>& prev,
                bool *m_open_vec,
                bool *m_closed_vec,
                const unsigned int v,
                const size_t n);
        void scan_edges (
                const DGraphEdge *edge,
                std::vector<double>& d,
                std::vector<double>& w,
                std::vector<int>& prev,
                bool *m_open_vec,
                const bool *m_closed_vec,
                const unsigned int &v0);
        void scan_edges_heur ( // with A* heuristic
                const DGraphEdge *edge,
                std::vector<double>& d,
                std::vector<double>& w,
                std::vector<int>& prev,
                bool *m_open_vec,
                const bool *m_closed_vec,
                const unsigned int &v0,
                const std::vector<double> &heur);

        void Dijkstra (
                std::vector<double>& d,
                std::vector<double>& w,
                std::vector<int>& prev,
                const unsigned int v0,
                const std::vector <unsigned int> &to_index);
        void DijkstraLimit (
                std::vector<double>& d,
                std::vector<double>& w,
                std::vector<int>& prev,
                const unsigned int v0,
                const double &dlim);
        void AStar (std::vector<double>& d,
                std::vector<double>& w,
                std::vector<int>& prev,
                const std::vector<double>& heur,
                const unsigned int v0,
                const std::vector <unsigned int> &to_index);
        void Dijkstra_set (std::vector <double>& d,
                std::vector<double>& w,
                std::vector<int>& prev,
                unsigned int v0);
        void Centrality_vertex (
                std::vector <double>& cent,
                const unsigned int s,
                const double dist_threshold);
        void Centrality_edge (
                std::vector <double>& cent,
                const unsigned int s,
                const unsigned int nedges,
                const double dist_threshold);

    private:
        Heap *m_heap;        // pointer: heap
        // Convert to vector<bool>? (save memory, might be a performace hit though)
        bool *m_open;           // array: frontier set state of vertices
        bool *m_closed;         // also for bi-dir

        std::shared_ptr<const DGraph> m_graph;    // pointer: directed graph    

        EdgeSet edge_set;
};

} // end namespace PF
