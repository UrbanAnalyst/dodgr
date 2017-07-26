#ifndef DIJKSTRA_H
#define DIJKSTRA_H
/* Dijkstra's Algorithm
 * ----------------------------------------------------------------------------
 * Author:  Shane Saunders
 */

class Heap;      // Heap
class HeapDesc;  // Heap descriptor
class DGraph;    // Graph

class IntHeap;      // Heap
class IntHeapDesc;  // Heap descriptor
class DGraphInt; // Graph

/* --- Dijkstra ---
 * Dijkstra's single-source algorithm.
 */
class Dijkstra {
  public:
    Dijkstra(int n, HeapDesc *heapD);
    ~Dijkstra();

    void init(const DGraph *g);
    //void run(float *d, unsigned int s = 0);
    void run(float *d, float *w, unsigned int s = 0);

  private:
    Heap *heap;        // pointer: heap
    bool *s;           // array: solution set state of vertices
    bool *f;           // array: frontier set state of vertices

    const DGraph *graph;    // pointer: directed graph    
};

// Implementation for Radix heaps with integer edge weights
class DijkstraInt {
  public:
    DijkstraInt(int n, IntHeapDesc *heapD);
    ~DijkstraInt();

    void init(const DGraphInt *g);
    //void run(float *d, unsigned int s = 0);
    void run(int *d, int *w, unsigned int s = 0);

  private:
    IntHeap *heap;        // pointer: heap
    bool *s;           // array: solution set state of vertices
    bool *f;           // array: frontier set state of vertices

    const DGraphInt *graph;    // pointer: directed graph    
};

#endif
