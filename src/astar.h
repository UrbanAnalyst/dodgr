#pragma once

#include <vector>
#include <memory>
#include <set>

class Heap;      // Heap
class HeapDesc;  // Heap descriptor
class DGraph;    // Graph

class Astar {
    public:
        Astar(unsigned int n, const HeapDesc& heapD, std::shared_ptr<const DGraph> g);
        ~Astar();
        
        // Disable copy/assign as will crash
        Astar(const Astar&) = delete;
        Astar& operator=(const Astar&) = delete;

        void init (std::shared_ptr<const DGraph> g);
        void run (std::vector<double>& d,
                std::vector<double>& w,
                std::vector<int>& prev,
                unsigned int s = 0);

    private:
        Heap *m_heap;        // pointer: heap
        // Convert to vector<bool>? (save memory, might be a performace hit though)
        bool *m_s;           // array: solution set state of vertices
        bool *m_f;           // array: frontier set state of vertices

        std::shared_ptr<const DGraph> m_graph;    // pointer: directed graph    
};
