#pragma once

// Modified from code by Shane Saunders

#include <vector>
#include <memory>
#include <set>
#include <unordered_set> // used in bidirected search

// inf_dbl used only in AStar2:
#include <limits>

#include "dgraph.h"

constexpr double INFINITE_INT =  std::numeric_limits<int>::max ();
constexpr double INFINITE_DOUBLE =  std::numeric_limits<double>::max ();

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
        if (lhs._wt == rhs._wt)
            return lhs._i < rhs._i;
        else
            return lhs._wt < rhs._wt;
    }
}; 

typedef std::set <DijkstraEdge, by_wt> EdgeSet;


class Heap;
class HeapDesc;
class DGraph;

class PathFinder {
    public:
        PathFinder (unsigned int n,
                const HeapDesc& heapD,
                std::shared_ptr<const DGraph> g,
                bool twoheap);
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
        void relax (
                const DGraphEdge *edge,
                std::vector<double>& d,
                std::vector<double>& w,
                std::vector<int>& prev,
                bool *m_open_vec,
                const unsigned int &v0,
                const unsigned int &target);
        void relax_heur (
                const DGraphEdge *edge,
                std::vector<double>& d,
                std::vector<double>& w,
                std::vector<int>& prev,
                bool *m_open_vec,
                const unsigned int &v0,
                const unsigned int &target,
                const std::vector<double> &heur,
                const double &dmax,
                const bool &reverse);
        void scan_edges (
                const DGraphEdge *edge,
                std::vector<double>& d,
                std::vector<double>& w,
                std::vector<int>& prev,
                bool *m_open_vec,
                const bool *m_closed_vec,
                const unsigned int &v0,
                const bool &use_heur,
                const std::vector<double> &heur,    // heuristic for A*
                const double &dmax,                 // used for reverse bidirectional
                const bool &reverse);               // reverse dir of bidirectional 

        void Dijkstra (
                std::vector<double>& d,
                std::vector<double>& w,
                std::vector<int>& prev,
                unsigned int s = 0);
        void AStar (std::vector<double>& d,
                std::vector<double>& w,
                std::vector<int>& prev,
                const std::vector<double>& heur,
                unsigned int v0);
        void AStar2 (std::vector<double>& d,
                std::vector<double>& w,
                std::vector<int>& prev,
                const std::vector<double>& heur,
                unsigned int v0, unsigned int v1);
        void Dijkstra_set (std::vector <double>& d,
                std::vector<double>& w,
                std::vector<int>& prev,
                unsigned int s = 0);

    private:
        bool _twoheap;
        Heap *m_heap;        // pointer: heap
        Heap *m_heap_rev;    // for reverse direction in bi-directional search
        // Convert to vector<bool>? (save memory, might be a performace hit though)
        bool *m_open;           // array: frontier set state of vertices
        bool *m_open_rev;       // for bi-directional search
        bool *m_closed;         // also for bi-dir
        bool *m_closed_rev;     // also for bi-dir

        std::shared_ptr<const DGraph> m_graph;    // pointer: directed graph    

        EdgeSet edge_set;
};
