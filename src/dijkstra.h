#ifndef DIJKSTRA_H
#define DIJKSTRA_H

#include <vector>
#include <memory>
#include <set>
#include <unordered_set> // used in bidirected search

// inf_dbl used only in astar2:
#include <limits>
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

/* Dijkstra's Algorithm
 * ----------------------------------------------------------------------------
 * Author:  Shane Saunders
 */

class Heap;      // Heap
class HeapDesc;  // Heap descriptor
class DGraph;    // Graph

/* --- Dijkstra ---
 * Dijkstra's single-source algorithm.
 */

class Dijkstra {
    public:
        Dijkstra (unsigned int n,
                const HeapDesc& heapD,
                std::shared_ptr<const DGraph> g,
                bool twoheap);
        ~Dijkstra();
        
        // Disable copy/assign as will crash
        Dijkstra(const Dijkstra&) = delete;
        Dijkstra& operator=(const Dijkstra&) = delete;

        void init (std::shared_ptr<const DGraph> g);
        void run (std::vector<double>& d,
                std::vector<double>& w,
                std::vector<int>& prev,
                unsigned int s = 0);
        void astar (std::vector<double>& d,
                std::vector<double>& w,
                std::vector<int>& prev,
                const std::vector<double>& heur,
                unsigned int v0);
        void astar2 (std::vector<double>& d,
                std::vector<double>& w,
                std::vector<int>& prev,
                const std::vector<double>& heur,
                unsigned int v0, unsigned int v1);
        void run_set (std::vector <double>& d,
                std::vector<double>& w,
                std::vector<int>& prev,
                unsigned int s = 0);

    private:
        bool _twoheap;
        Heap *m_heap;        // pointer: heap
        Heap *m_heap_rev;    // for reverse direction in bi-directional search
        // Convert to vector<bool>? (save memory, might be a performace hit though)
        bool *m_settled;        // array: solution set state of vertices
        bool *m_open;           // array: frontier set state of vertices
        bool *m_open2;          // for bi-directional search

        std::shared_ptr<const DGraph> m_graph;    // pointer: directed graph    

        EdgeSet edge_set;
};

#endif
