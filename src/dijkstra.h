#ifndef DIJKSTRA_H
#define DIJKSTRA_H

#include <vector>
#include <memory>

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
        Dijkstra(unsigned int n, const HeapDesc& heapD, std::shared_ptr<const DGraph> g);
        ~Dijkstra();
        
        // Disable copy/assign as will crash
        Dijkstra(const Dijkstra&) = delete;
        Dijkstra& operator=(const Dijkstra&) = delete;

        void init(std::shared_ptr<const DGraph> g);

        void run(std::vector<double>& d, std::vector<double>& w, std::vector<int>& prev, unsigned int s = 0);

    private:
        Heap *m_heap;        // pointer: heap
        // Convert to vector<bool>? (save memory, might be a performace hit though)
        bool *m_s;           // array: solution set state of vertices
        bool *m_f;           // array: frontier set state of vertices

        std::shared_ptr<const DGraph> m_graph;    // pointer: directed graph    
};

#endif
