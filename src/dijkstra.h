#ifndef DIJKSTRA_H
#define DIJKSTRA_H
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
        Dijkstra(int n, HeapDesc *heapD);
        ~Dijkstra();

        void init(const DGraph *g);
        void run(float *d, float *w, int *prev, unsigned int s = 0);

    private:
        Heap *heap;        // pointer: heap
        bool *s;           // array: solution set state of vertices
        bool *f;           // array: frontier set state of vertices

        const DGraph *graph;    // pointer: directed graph    
};

#endif
