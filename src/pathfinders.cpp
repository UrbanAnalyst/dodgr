#include "pathfinders.h"
#include "heaps/heap.h"

#include <algorithm> // std::fill

// Modified from code by Shane Saunders

// @param n Number of vertices in graph
// @param heapD The type of heap used
// @param g A DGraph object
// @param twoheap If `TRUE`, uses a bi-directional search.
PF::PathFinder::PathFinder(size_t n,
        const HeapDesc& heapD,
        std::shared_ptr<const DGraph> g)
{
    m_heap = heapD.newInstance (n);
    m_closed = new bool [n];
    m_open = new bool [n];
    init (g);
}

PF::PathFinder::~PathFinder() {
    delete [] m_open;
    delete [] m_closed;
    delete m_heap;
}

void PF::PathFinder::init(std::shared_ptr<const DGraph> g) {
    m_graph = g;
}

void PF::PathFinder::init_arrays (
        std::vector <double>& d,
        std::vector <double>& w,
        std::vector <long int>& prev,
        bool *m_open_vec,
        bool *m_closed_vec,
        const size_t v,
        const size_t n)
{
    std::fill (w.begin (), w.end (), INFINITE_DOUBLE);
    std::fill (d.begin (), d.end (), INFINITE_DOUBLE);
    std::fill (prev.begin (), prev.end (), INFINITE_INT);
    w [v] = 0.0;
    d [v] = 0.0;
    prev [v] = -1;

    std::fill (m_open_vec, m_open_vec + n, false);
    std::fill (m_closed_vec, m_closed_vec + n, false);
    m_open_vec [v] = true;
}

void PF::PathFinder::scan_edges (const DGraphEdge *edge,
        std::vector<double>& d,
        std::vector<double>& w,
        std::vector<long int>& prev,
        bool *m_open_vec,
        const bool *m_closed_vec,
        const size_t &v0)
{
    while (edge) {
        size_t et = edge->target;
        if (!m_closed_vec [et])
        {
            double wt = w [v0] + edge->wt;
            if (wt < w [et]) {
                d [et] = d [v0] + edge->dist;
                w [et] = wt;
                prev [et] = static_cast <int> (v0);

                if (m_open_vec [et]) {
                    m_heap->decreaseKey(et, wt);
                }
                else {
                    m_heap->insert (et, wt);
                    m_open_vec [et] = true;
                }
            } else
                m_closed [et] = true;
        }
        edge = edge->nextOut;
    }
}

void PF::PathFinder::scan_edges_heur (const DGraphEdge *edge,
        std::vector<double>& d,
        std::vector<double>& w,
        std::vector<long int>& prev,
        bool *m_open_vec,
        const bool *m_closed_vec,
        const size_t &v0,
        const std::vector<double> &heur)    // heuristic for A*
{
    while (edge) {
        size_t et = edge->target;
        if (!m_closed_vec [et])
        {
            double wt = w [v0] + edge->wt;
            if (wt < w [et]) {
                d [et] = d [v0] + edge->dist;
                w [et] = wt;
                prev [et] = static_cast <int> (v0);

                if (m_open_vec [et]) {
                    m_heap->decreaseKey(et, wt + heur [et] - heur [v0]);
                }
                else {
                    m_heap->insert (et, wt + heur [et] - heur [v0]);
                    m_open_vec [et] = true;
                }
            } else
                m_closed [et] = true;
        }
        edge = edge->nextOut;
    }
}

/**********************************************************************
 ************************  PATH ALGORITHMS    *************************
 **********************************************************************/

void PF::PathFinder::Dijkstra (
        std::vector<double>& d,
        std::vector<double>& w,
        std::vector<long int>& prev,
        const size_t v0,
        const std::vector <size_t> &to_index)
{
    const DGraphEdge *edge;

    const size_t n = m_graph->nVertices();
    const std::vector<DGraphVertex>& vertices = m_graph->vertices();

    PF::PathFinder::init_arrays (d, w, prev, m_open, m_closed, v0, n);
    m_heap->insert (v0, 0.0);

    size_t n_reached = 0;
    const size_t n_targets = to_index.size ();
    bool *is_target = new bool [n];
    std::fill (is_target, is_target + n, false);
    for (auto t: to_index)
        is_target [t] = true;

    while (m_heap->nItems() > 0) {
        size_t v = m_heap->deleteMin();

        m_closed [v] = true;
        m_open [v] = false;

        edge = vertices [v].outHead;
        scan_edges (edge, d, w, prev, m_open, m_closed, v);

        if (is_target [v])
            n_reached++;
        if (n_reached == n_targets)
            break;
    } // end while nItems > 0
    delete [] is_target;
}

// Modified Dijkstra that stops as soon as first target is reached. Used to find
// minimal distance to nearest one of set of possible targets.
void PF::PathFinder::DijkstraNearest (
        std::vector<double>& d,
        std::vector<double>& w,
        std::vector<long int>& prev,
        const size_t v0,
        const std::vector <size_t> &to_index)
{
    const DGraphEdge *edge;

    const size_t n = m_graph->nVertices();
    const std::vector<DGraphVertex>& vertices = m_graph->vertices();

    PF::PathFinder::init_arrays (d, w, prev, m_open, m_closed, v0, n);
    m_heap->insert (v0, 0.0);

    bool *is_target = new bool [n];
    std::fill (is_target, is_target + n, false);
    for (auto t: to_index)
        is_target [t] = true;

    while (m_heap->nItems() > 0) {
        size_t v = m_heap->deleteMin();

        m_closed [v] = true;
        m_open [v] = false;

        edge = vertices [v].outHead;
        scan_edges (edge, d, w, prev, m_open, m_closed, v);

        if (is_target [v])
            break;
    } // end while nItems > 0
    delete [] is_target;
}

// Modified pathfinder only out to specified distance limit. Done as separate
// routine to avoid costs of the `if` clause in general non-limited case.
void PF::PathFinder::DijkstraLimit (
        std::vector<double>& d,
        std::vector<double>& w,
        std::vector<long int>& prev,
        const size_t v0,
        const double &dlim)
{
    const DGraphEdge *edge;

    const size_t n = m_graph->nVertices();
    const std::vector<DGraphVertex>& vertices = m_graph->vertices();

    PF::PathFinder::init_arrays (d, w, prev, m_open, m_closed, v0, n);
    m_heap->insert (v0, 0.0);

    while (m_heap->nItems() > 0) {
        size_t v = m_heap->deleteMin();

        m_closed [v] = true;
        m_open [v] = false;

        // explore the OUT set of v only if distances are < threshold
        bool explore = false;
        edge = vertices [v].outHead;
        while (edge) {
            if ((d [v] + edge->dist) <= dlim)
            {
                explore = true;
                break;
            }
            edge = edge->nextOut;
        }

        if (explore)
        {
            edge = vertices [v].outHead;
            scan_edges (edge, d, w, prev, m_open, m_closed, v);
        }
    } // end while nItems > 0
}

void PF::PathFinder::AStar (std::vector<double>& d,
        std::vector<double>& w,
        std::vector<long int>& prev,
        const std::vector<double>& heur,
        const size_t v0,
        const std::vector <size_t> &to_index)
{
    const DGraphEdge *edge;

    const size_t n = m_graph->nVertices();
    const std::vector<DGraphVertex>& vertices = m_graph->vertices();

    PF::PathFinder::init_arrays (d, w, prev, m_open, m_closed, v0, n);
    m_heap->insert (v0, heur [v0]);

    size_t n_reached = 0;
    const size_t n_targets = to_index.size ();
    bool *is_target = new bool [n];
    std::fill (is_target, is_target + n, false);
    for (auto t: to_index)
        is_target [t] = true;

    while (m_heap->nItems() > 0) {
        size_t v = m_heap->deleteMin();

        m_closed [v] = true;
        m_open [v] = false;

        edge = vertices [v].outHead;
        scan_edges_heur (edge, d, w, prev, m_open, m_closed, v, heur);

        if (is_target [v])
            n_reached++;
        if (n_reached == n_targets)
            break;
    } // end while m_heap->nItems
    delete [] is_target;
}

// Dijkstra with C++ std::set, modified from the above to use EdgeSet edge_set
// instead of m_heap, and so all routines hard-coded here.
void PF::PathFinder::Dijkstra_set (std::vector<double>& d,
        std::vector<double>& w,
        std::vector<long int>& prev,
        size_t v0)
{
    const DGraphEdge *edge;

    const size_t n = m_graph->nVertices();
    const std::vector<DGraphVertex>& vertices = m_graph->vertices();

    PF::PathFinder::init_arrays (d, w, prev, m_open, m_closed, v0, n);
    m_heap->insert(v0, 0.0);

    edge_set.insert (DijkstraEdge (0.0, v0));

    while (edge_set.size () > 0) {
        EdgeSet::iterator ei = edge_set.begin();
        size_t v = ei->geti();
        edge_set.erase (ei);

        m_closed [v] = true;
        m_open [v] = false;

        edge = vertices [v].outHead;
        while (edge) {
            size_t et = edge->target;

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
