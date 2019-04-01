#include "dijkstra.h"
#include "dgraph.h"
#include "heaps/heap.h"

#include <algorithm> // std::fill

/* Dijkstra's Algorithm
 * ----------------------------------------------------------------------------
 * Author:  Shane Saunders, modified to dual-weighted graphs by Mark Padgham
 */

/*--- Dijkstra -----------------------------------------------------------*/

/* - Constructor -
 * Allocate the algorithm for use on graphs of n vertices.  The parameter heapD
 * (points to a heap descriptor object) specifies the heap to be used by
 * Dijkstra's algorithm.
 */
Dijkstra::Dijkstra(unsigned int n,
        const HeapDesc& heapD,
        std::shared_ptr<const DGraph> g,
        bool twoheap) 
{
    m_heap = heapD.newInstance(n);
    m_heap_rev = heapD.newInstance(n);
    m_settled = new bool[n];
    m_open = new bool[n];
    _twoheap = twoheap;
    if (twoheap)
    {
        m_open2 = new bool[n];
        m_settled2 = new bool[n];
    }
    init(g);
}

/* - Destructor - */
Dijkstra::~Dijkstra() {
    delete [] m_settled;
    delete [] m_open;
    if (_twoheap)
    {
        delete [] m_open2;
        delete [] m_settled2;
    }
    delete m_heap;
    delete m_heap_rev;
}

/* - init() -
 * Initialise the algorithm for use with the graph pointed to by g.
 */
void Dijkstra::init(std::shared_ptr<const DGraph> g) {
    m_graph = g;
}

/* - run() -
 * Run the algorithm, computing single-source from the starting vertex v0.
 * This assumes that the array d has been initialised with d[v] = INFINITE_DIST
 * for all vertices v != v0.
 */

void Dijkstra::run (std::vector<double>& d,
        std::vector<double>& w,
        std::vector<int>& prev,
        unsigned int v0)
{
    /* indexes, counters, pointers */
    const DGraphEdge *edge;

    /*** initialisation ***/

    /* optimise access to the data structures allocated for the algorithm */
    const unsigned int n = m_graph->nVertices();
    const std::vector<DGraphVertex>& vertices = m_graph->vertices();

    /* initialise all vertices as unexplored */
    std::fill (m_settled, m_settled + n, false);
    std::fill (m_open, m_open + n, false);

    /* place v0 into the frontier set with a distance of zero */
    w [v0] = 0.0;
    d [v0] = 0.0;
    prev [v0] = -1;
    m_heap->insert(v0, 0.0);
    m_open [v0] = true;

    /* repeatedly update distances from the minimum remaining trigger vertex */
    while (m_heap->nItems() > 0) {
        /* delete the vertex in frontier that has minimum distance */
        unsigned int v = m_heap->deleteMin();

        /* the selected vertex moves from the frontier to the solution set */
        m_settled [v] = true;
        m_open [v] = false;

        /* explore the OUT set of v */
        edge = vertices [v].outHead;
        while (edge) {
            unsigned int et = edge->target;

            if (!m_settled [et]) {
                double wt = w [v] + edge->wt;
                if (wt < w [et]) {
                    d [et] = d [v] + edge->dist;
                    w [et] = wt;
                    prev [et] = static_cast <int> (v);
                    if (m_open [et]) {
                      m_heap->decreaseKey(et, wt);
                    }
                    else {
                      m_heap->insert (et, wt);
                      m_open [et] = true;
                    }
                }
            }

            edge = edge->nextOut;
        } /* while */
    } /* while */
}

void Dijkstra::astar (std::vector<double>& d,
        std::vector<double>& w,
        std::vector<int>& prev,
        const std::vector<double>& heur,
        unsigned int v0)
{
    /* indexes, counters, pointers */
    const DGraphEdge *edge;

    /*** initialisation ***/

    /* optimise access to the data structures allocated for the algorithm */
    const unsigned int n = m_graph->nVertices();
    const std::vector<DGraphVertex>& vertices = m_graph->vertices();

    /* initialise all vertices as unexplored */
    std::fill (m_settled, m_settled + n, false);
    std::fill (m_open, m_open + n, false);

    /* place v0 into the frontier set with a distance of zero */
    w [v0] = 0.0;
    d [v0] = 0.0;
    prev [v0] = -1;
    //m_heap->insert(v0, 0.0);
    m_heap->insert(v0, heur [v0]);
    m_open [v0] = true;

    /* repeatedly update distances from the minimum remaining trigger vertex */
    while (m_heap->nItems() > 0) {
        /* delete the vertex in frontier that has minimum distance */
        unsigned int v = m_heap->deleteMin();

        /* the selected vertex moves from the frontier to the solution set */
        m_settled [v] = true;
        m_open [v] = false;

        /* explore the OUT set of v */
        edge = vertices [v].outHead;
        while (edge) {
            unsigned int et = edge->target;

            if (!m_settled [et]) {
                double wt = w [v] + edge->wt;
                if (wt < w [et]) {
                    d [et] = d [v] + edge->dist;
                    w [et] = wt;
                    prev [et] = static_cast <int> (v);
                    if (m_open [et]) {
                      m_heap->decreaseKey (et, wt + heur [et] - heur [v]);
                    }
                    else {
                      m_heap->insert (et, wt + heur [et] - heur [v]);
                      m_open [et] = true;
                    }
                }
            }

            edge = edge->nextOut;
        } /* while */
    } /* while */
}

// bi-directional A*
void Dijkstra::astar2 (std::vector<double>& d,
        std::vector<double>& w,
        std::vector<int>& prev,
        const std::vector<double>& heur,
        unsigned int v0, unsigned int v1)
{
    /* indexes, counters, pointers */
    const DGraphEdge *edge;

    /*** initialisation ***/

    /* optimise access to the data structures allocated for the algorithm */
    const unsigned int n = m_graph->nVertices();
    const std::vector<DGraphVertex>& vertices = m_graph->vertices();

    /* initialise all vertices as unexplored */
    std::fill (m_settled, m_settled + n, false);
    std::fill (m_open, m_open + n, false);

    // new vectors for reverse-direction scans
    bool *frontier = new bool [n];
    std::vector <double> d_rev (n), w_rev (n);

    /* place v0 into the frontier set with a distance of zero */
    w [v0] = 0.0;
    d [v0] = 0.0;
    prev [v0] = -1;
    //m_heap->insert(v0, 0.0);
    m_heap->insert(v0, heur [v0]);
    m_open [v0] = true;

    double dmax = *std::max_element (heur.begin (), heur.end ());
    m_heap_rev->insert (v1, dmax - heur [v1]);
    m_open2 [v1] = true;

    /* repeatedly update distances from the minimum remaining trigger vertex */
    while (m_heap->nItems() > 0) {
        /* delete the vertex in frontier that has minimum distance */
        unsigned int v = m_heap->deleteMin();

        /* the selected vertex moves from the frontier to the solution set */
        m_settled [v] = true;
        m_open [v] = false;

        /* explore the OUT set of v */
        edge = vertices [v].outHead;
        while (edge) {
            unsigned int et = edge->target;

            if (!m_settled [et]) {
                double wt = w [v] + edge->wt;
                if (wt < w [et]) {
                    d [et] = d [v] + edge->dist;
                    w [et] = wt;
                    prev [et] = static_cast <int> (v);
                    if (m_open [et]) {
                      m_heap->decreaseKey (et, wt + heur [et] - heur [v]);
                    }
                    else {
                      m_heap->insert (et, wt + heur [et] - heur [v]);
                      m_open [et] = true;
                    }
                }
            }

            edge = edge->nextOut;
        } /* while */
    } /* while */

    delete [] frontier;
}

// Only sightly modified from the above, to use EdgeSet edge_set instead of
// m_heap.
void Dijkstra::run_set (std::vector<double>& d,
        std::vector<double>& w,
        std::vector<int>& prev,
        unsigned int v0)
{
    /* indexes, counters, pointers */
    const DGraphEdge *edge;

    /*** initialisation ***/

    /* optimise access to the data structures allocated for the algorithm */
    const unsigned int n = m_graph->nVertices();
    const std::vector<DGraphVertex>& vertices = m_graph->vertices();

    /* initialise all vertices as unexplored */
    std::fill (m_settled, m_settled + n, false);
    std::fill (m_open, m_open + n, false);

    /* place v0 into the frontier set with a distance of zero */
    w [v0] = 0.0;
    d [v0] = 0.0;
    prev [v0] = -1;
    m_heap->insert(v0, 0.0);
    m_open [v0] = true;

    edge_set.insert (DijkstraEdge (0.0, v0));

    /* repeatedly update distances from the minimum remaining trigger vertex */
    while (edge_set.size () > 0) {
        /* delete the vertex in frontier that has minimum distance */
        EdgeSet::iterator ei = edge_set.begin();
        unsigned int v = ei->geti();
        //double weight = ei->getw();
        edge_set.erase (ei);

        /* the selected vertex moves from the frontier to the solution set */
        m_settled [v] = true;
        m_open [v] = false;

        /* explore the OUT set of v */
        edge = vertices [v].outHead;
        while (edge) {
            unsigned int et = edge->target;

            if (!m_settled [et]) {
                double wt = w [v] + edge->wt;
                if (wt < w [et]) {
                    d [et] = d [v] + edge->dist;
                    double wold = w [et];
                    w [et] = wt;
                    prev [et] = static_cast <int> (v);
                    if (m_open [et])
                    {
                        DijkstraEdge de (wold, et);
                        if (edge_set.find (de) != edge_set.end ())
                            edge_set.erase (edge_set.find (de));
                    } else
                        m_open [et] = true;
                    edge_set.insert (DijkstraEdge (w [et], et));
                }
            }

            edge = edge->nextOut;
        } /* while */
    } /* while */
}
