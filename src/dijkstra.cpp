#include "dijkstra.h"
#include "dgraph.h"
#include "heap.h"
/* Dijkstra's Algorithm
 * ----------------------------------------------------------------------------
 * Author:  Shane Saunders
 */

/*--- Dijkstra -----------------------------------------------------------*/

/* - Constructor -
 * Allocate the algorithm for use on graphs of n vertices.  The parameter heapD
 * (points to a heap descriptor object) specifies then heap to be used by
 * Dijkstra's algorithm.
 */
Dijkstra::Dijkstra(int n, HeapDesc *heapD)
{
    heap = heapD->newInstance(n);    
    s = new bool[n];
    f = new bool[n];
}

/* - Destructor - */
Dijkstra::~Dijkstra() {
    delete [] s;
    delete [] f;
    delete heap;
}

/* - init() -
 * Initialise the algorithm for use with the graph pointed to by g.
 */
void Dijkstra::init(const DGraph *g) {
    graph = g;
}

/* - run() -
 * Run the algorithm, computing single-source from the starting vertex v0.
 * This assumes that the array d has been initialised with d[v] = INFINITE_DIST
 * for all vertices v != v0.
 */
void Dijkstra::run(float *d, float *w, unsigned int v0)
{
    /* indexes, counters, pointers */
    const DGraphEdge *edge;

    
  /*** initialisation ***/

    /* optimise access to the data structures allocated for the algorithm */
    const unsigned int n = graph->nVertices;
    const DGraphVertex *vertices = graph->vertices;
    
  /*** algorithm ***/

    /* initialise all vertices as unexplored */
    for (unsigned int v = 0; v < n; v++)
    {
        s [v] = false;
        f [v] = false; 
    }

    /* place v0 into the frontier set with a distance of zero */
    w [v0] = 0.0;
    d [v0] = 0.0;
    heap->insert(v0, 0.0);
    f [v0] = true;

    /* repeatedly update distances from the minimum remaining trigger vertex */
    while (heap->nItems() > 0) {
        /* delete the vertex in frontier that has minimum distance */
        unsigned int v = heap->deleteMin();

        /* the selected vertex moves from the frontier to the solution set */
        s [v] = true;
        f [v] = false;

        /* explore the OUT set of v */
        edge = vertices [v].outHead;
        while (edge) {
            unsigned int et = edge->target;

            if (!s[et]) {
                float wt = w [v] + edge->wt;
                if (wt < w [et]) {
                    d [et] = d [v] + edge->dist;
                    w [et] = wt;
                    if (f [et]) {
                        heap->decreaseKey(et, wt);
                    }
                    else {
                        heap->insert (et, wt);
                        f [et] = true;
                    }
                }
            }

            edge = edge->nextOut;
        } /* while */

    } /* while */
}

/*--- Dijkstra Int -------------------------------------------------------*/
// for integer edge weights with the Radix Heap

DijkstraInt::DijkstraInt(int n, HeapDesc *heapD)
{
    heap = heapD->newInstance(n);    
    s = new bool[n];
    f = new bool[n];
}

/* - Destructor - */
DijkstraInt::~DijkstraInt() {
    delete [] s;
    delete [] f;
    delete heap;
}

/* - init() -
 * Initialise the algorithm for use with the graph pointed to by g.
 */
void DijkstraInt::init(const DGraph *g) {
    graph = g;
}

/* - run() -
 * Run the algorithm, computing single-source from the starting vertex v0.
 * This assumes that the array d has been initialised with d[v] = INFINITE_DIST
 * for all vertices v != v0.
 */
void DijkstraInt::run(int *d, int *w, unsigned int v0)
{
    /* indexes, counters, pointers */
    const DGraphEdge *edge;

    
  /*** initialisation ***/

    /* optimise access to the data structures allocated for the algorithm */
    const unsigned int n = graph->nVertices;
    const DGraphVertex *vertices = graph->vertices;
    
  /*** algorithm ***/

    /* initialise all vertices as unexplored */
    for (unsigned int v = 0; v < n; v++)
    {
        s [v] = false;
        f [v] = false; 
    }

    /* place v0 into the frontier set with a distance of zero */
    w [v0] = 0.0;
    d [v0] = 0.0;
    heap->insert(v0, 0.0);
    f [v0] = true;

    /* repeatedly update distances from the minimum remaining trigger vertex */
    while (heap->nItems() > 0) {
        /* delete the vertex in frontier that has minimum distance */
        unsigned int v = heap->deleteMin();

        /* the selected vertex moves from the frontier to the solution set */
        s [v] = true;
        f [v] = false;

        /* explore the OUT set of v */
        edge = vertices [v].outHead;
        while (edge) {
            unsigned int et = edge->target;

            if (!s[et]) {
                int wt = w [v] + edge->wt;
                if (wt < w [et]) {
                    d [et] = d [v] + edge->dist;
                    w [et] = wt;
                    if (f [et]) {
                        heap->decreaseKey(et, wt);
                    }
                    else {
                        heap->insert (et, wt);
                        f [et] = true;
                    }
                }
            }

            edge = edge->nextOut;
        } /* while */

    } /* while */
}
