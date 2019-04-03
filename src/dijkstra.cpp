#include "dijkstra.h"
#include "dgraph.h"
#include "heaps/heap.h"

#include <algorithm> // std::fill

#include <Rcpp.h> // TODO: Delete!

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
    m_closed = new bool[n];
    m_open = new bool[n];
    _twoheap = twoheap;
    if (_twoheap)
    {
        m_heap_rev = heapD.newInstance(n);
        m_open2 = new bool[n];
        m_closed2 = new bool[n];
    }
    init(g);
}

/* - Destructor - */
Dijkstra::~Dijkstra() {
    delete [] m_open;
    if (_twoheap)
    {
        delete [] m_open2;
        delete [] m_closed;
        delete [] m_closed2;
        delete m_heap_rev;
    }
    delete m_heap;
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
    std::fill (m_closed, m_closed + n, false);
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
        m_closed [v] = true;
        m_open [v] = false;

        /* explore the OUT set of v */
        edge = vertices [v].outHead;
        while (edge) {
            unsigned int et = edge->target;

            if (!m_closed [et]) {
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
    const DGraphEdge *edge;

    const unsigned int n = m_graph->nVertices();
    const std::vector<DGraphVertex>& vertices = m_graph->vertices();

    std::fill (m_closed, m_closed + n, false);
    std::fill (m_open, m_open + n, false);

    w [v0] = 0.0;
    d [v0] = 0.0;
    prev [v0] = -1;
    m_heap->insert(v0, heur [v0]);
    m_open [v0] = true;

    while (m_heap->nItems() > 0) {
        unsigned int v = m_heap->deleteMin();

        m_closed [v] = true;
        m_open [v] = false;

        edge = vertices [v].outHead;
        while (edge) {
            unsigned int et = edge->target;

            if (!m_closed [et]) {
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
        } // end while edge
    } // end while m_heap->nItems
}

// bi-directional A*
void Dijkstra::astar2 (std::vector<double>& d,
        std::vector<double>& w,
        std::vector<int>& prev,
        const std::vector<double>& heur,
        unsigned int v0, unsigned int v1)
{
    const DGraphEdge *edge;

    const unsigned int n = m_graph->nVertices();
    const std::vector<DGraphVertex>& vertices = m_graph->vertices();

    std::fill (m_open, m_open + n, false);
    std::fill (m_open2, m_open2 + n, false);
    std::fill (m_closed, m_closed + n, false);
    std::fill (m_closed2, m_closed2 + n, false);

    std::unordered_set <unsigned int> frontier, backward;

    std::vector <double> d_rev (n, INFINITE_DOUBLE),
        w_rev (n, INFINITE_DOUBLE),
        prev_rev (n, INFINITE_DOUBLE);

    // m_heap holds heuristic estimates from source v0 to target vertex v1
    // m_heap_rev holds reverse: from target v1 to source v0
    w [v0] = 0.0;
    d [v0] = 0.0;
    prev [v0] = -1;
    m_heap->insert(v0, heur [v1]);
    m_open [v0] = true;

    double dmax = heur [v1];
    m_heap_rev->insert (v1, dmax - heur [v0]); // = heur [v1]
    m_open2 [v1] = true;
    w_rev [v1] = 0.0;
    d_rev [v1] = 0.0;

    while (m_heap->nItems() > 0 && m_heap_rev->nItems() > 0) {
        //if (m_heap->getmin () <= m_heap_rev->getmin ())
        if (m_heap->nItems () >= m_heap_rev->nItems ())
        {
            unsigned int v = m_heap->deleteMin();
            m_open [v] = false;
            m_closed [v] = true;

            if (m_closed2 [v])
            {
                frontier.emplace (v);
            } else
            {
                edge = vertices [v].outHead;
                while (edge) {
                    unsigned int et = edge->target;

                    if (!m_closed [et]) {
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
                } // end while edge
            } // end else !m_closed2
        } else // else m_heap_rev.nItems > m_heap.nItems
        {
            unsigned int v = m_heap_rev->deleteMin();
            m_open2 [v] = false;
            m_closed2 [v] = true;

            backward.emplace (v);

            if (m_closed [v])
            {
                //frontier.emplace (v);
            } else
            {
                edge = vertices [v].inHead;
                while (edge) {
                    unsigned int et = edge->source;

                    if (!m_closed2 [et]) {
                        double wt = w_rev [v] + edge->wt;
                        if (wt < w_rev [et]) {
                            d_rev [et] = d_rev [v] + edge->dist;
                            w_rev [et] = wt;
                            prev [et] = static_cast <int> (v);
                            if (m_open2 [et]) {
                                m_heap_rev->decreaseKey (et,
                                        wt + dmax - (heur [et] - heur [v]));
                            }
                            else {
                                m_heap_rev->insert (et,
                                        wt + dmax - (heur [et] - heur [v]));
                                m_open2 [et] = true;
                            }
                        }
                    }

                    edge = edge->nextIn;
                } // end while edge
            } // end else !m_closed [v]
        } // end else m_heap_rev.nItems > m_heap.nItems
    } // end while m_heap->nItems

    /*
    int count1 = 0, count2 = 0;
    for (int i = 0; i < n; i++)
    {
        if (m_open [i])
        {
            count1++;
        }
        if (m_open2 [i])
        {
            count2++;
        }
    }
    Rcpp::Rcout << "(f, b) have (" << count1 << ", " << count2 <<
        ") nodes still open" << std::endl;

    // Reconstruct all distances from the frontier
    Rcpp::Rcout << "(frontier, backward) = (" << frontier.size () <<
        ", " << backward.size () << ") / " << n << std::endl;
    */

    // TODO: This is wrong
    for (auto fr:frontier)
        for (auto b:backward)
        {
            double wtemp = w [fr] + w_rev [b];
            if (wtemp < w [b])
            {
                if (b == 1148)
                    Rcpp::Rcout << "d [" << b << "] <- " <<
                        d [fr] << " + " << d_rev [b] << " = " <<
                        d [fr] + d_rev [b] << std::endl;
                w [b] = wtemp;
                d [b] = d [fr] + d_rev [b];
            }
        }
}

// Only sightly modified from the above, to use EdgeSet edge_set instead of
// m_heap.
void Dijkstra::run_set (std::vector<double>& d,
        std::vector<double>& w,
        std::vector<int>& prev,
        unsigned int v0)
{
    const DGraphEdge *edge;

    const unsigned int n = m_graph->nVertices();
    const std::vector<DGraphVertex>& vertices = m_graph->vertices();

    std::fill (m_closed, m_closed + n, false);
    std::fill (m_open, m_open + n, false);

    w [v0] = 0.0;
    d [v0] = 0.0;
    prev [v0] = -1;
    m_heap->insert(v0, 0.0);
    m_open [v0] = true;

    edge_set.insert (DijkstraEdge (0.0, v0));

    while (edge_set.size () > 0) {
        EdgeSet::iterator ei = edge_set.begin();
        unsigned int v = ei->geti();
        edge_set.erase (ei);

        m_closed [v] = true;
        m_open [v] = false;

        edge = vertices [v].outHead;
        while (edge) {
            unsigned int et = edge->target;

            if (!m_closed [et]) {
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
        } // end while edge
    } // end while edge_set.size
}
