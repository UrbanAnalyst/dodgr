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
        std::shared_ptr<const DGraph> g,
        bool twoheap) : m_is_bidirectional(twoheap)
{
    m_heap = heapD.newInstance (n);
    m_closed = new bool [n];
    m_open = new bool [n];
    if (twoheap) {
        m_heap_rev = heapD.newInstance (n);
        m_closed2 = new bool [n];
        m_open2 = new bool [n];
    } else {
        m_heap_rev = nullptr;
        m_closed2 = nullptr;
        m_open2 = nullptr;
    }
    init (g);
}

PF::PathFinder::~PathFinder() {
    delete [] m_open;
    delete [] m_closed;
    delete m_heap;
    if (m_is_bidirectional) {
        delete [] m_open2;
        delete [] m_closed2;
        delete m_heap_rev;
    }
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
                    m_heap->decreaseKey(et, wt + heur [et]);
                }
                else {
                    m_heap->insert (et, wt + heur [et]);
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
void PF::PathFinder::scan_edges_heur_rev (const DGraphEdge *edge,
        std::vector<double>& d,
        std::vector<double>& w,
        std::vector<long int>& prev,
        bool *m_open_vec,
        const bool *m_closed_vec,
        const size_t &v0,
        const std::vector<double> &heur,
        const double h_max)
{
    while (edge) {
        size_t et = edge->source;
        if (!m_closed_vec [et])
        {
            double wt = w [v0] + edge->wt;
            if (wt < w [et]) {
                d [et] = d [v0] + edge->dist;
                w [et] = wt;
                prev [et] = static_cast <int> (v0);

                double heur_et = h_max - heur[et];
                double priority = wt + heur_et;

                if (m_open_vec [et]) {
                    m_heap_rev->decreaseKey(et, priority);
                }
                else {
                    m_heap_rev->insert (et, priority);
                    m_open_vec [et] = true;
                }
            } else
                m_closed2 [et] = true;
        }
        edge = edge->nextIn;
    }
}

void PF::PathFinder::AStar2 (std::vector<double>& d,
        std::vector<double>& w,
        std::vector<long int>& prev,
        const std::vector<double>& heur,
        const size_t v0,
        const size_t v1,
        std::vector<double>& d_rev,
        std::vector<double>& w_rev,
        std::vector<long int>& prev_rev)
{
    const size_t n = m_graph->nVertices();
    const std::vector<DGraphVertex>& vertices = m_graph->vertices();

    if (v0 == v1) {
        d[v0] = 0.0;
        w[v0] = 0.0;
        return;
    }

    PF::PathFinder::init_arrays (d, w, prev, m_open, m_closed, v0, n);
    PF::PathFinder::init_arrays (d_rev, w_rev, prev_rev, m_open2, m_closed2, v1, n);

    m_heap->insert (v0, heur [v0]);
    double h_max = heur[v0];
    m_heap_rev->insert (v1, h_max - heur [v1]);

    size_t meeting_vertex = static_cast<size_t>(-1);
    double meeting_distance = INFINITE_DOUBLE;

    while (m_heap->nItems() > 0 || m_heap_rev->nItems() > 0) {
        double min_f = (m_heap->nItems() > 0 ? m_heap->getmin() : INFINITE_DOUBLE);
        double min_r = (m_heap_rev->nItems() > 0 ? m_heap_rev->getmin() : INFINITE_DOUBLE);

        // Use early termination
        if (min_f + min_r >= meeting_distance + h_max) {
            break;
        }

        if (m_heap->nItems() > 0 && (m_heap_rev->nItems() == 0 || min_f <= min_r)) {
            size_t v = m_heap->deleteMin();
            m_closed [v] = true;
            m_open [v] = false;

            if (w_rev[v] < INFINITE_DOUBLE && w[v] + w_rev[v] < meeting_distance) {
                meeting_vertex = v;
                meeting_distance = w[v] + w_rev[v];
            }

            const DGraphEdge *edge = vertices [v].outHead;
            scan_edges_heur (edge, d, w, prev, m_open, m_closed, v, heur);
        } else {
            size_t v = m_heap_rev->deleteMin();
            m_closed2 [v] = true;
            m_open2 [v] = false;

            if (w[v] < INFINITE_DOUBLE && w[v] + w_rev[v] < meeting_distance) {
                meeting_vertex = v;
                meeting_distance = w[v] + w_rev[v];
            }

            const DGraphEdge *edge = vertices [v].inHead;
            scan_edges_heur_rev (edge, d_rev, w_rev, prev_rev, m_open2, m_closed2, v, heur, h_max);
        }
    }

    // Safety check: find the absolute best meeting point over all visited nodes
    double best_w = INFINITE_DOUBLE;
    size_t best_meet = static_cast<size_t>(-1);
    for (size_t i = 0; i < n; i++) {
        if (w[i] < INFINITE_DOUBLE && w_rev[i] < INFINITE_DOUBLE) {
            if (w[i] + w_rev[i] < best_w) {
                best_w = w[i] + w_rev[i];
                best_meet = i;
            }
        }
    }

    if (best_meet != static_cast<size_t>(-1)) {
        w[v1] = best_w;
        d[v1] = d[best_meet] + d_rev[best_meet];

        size_t curr = best_meet;
        while (curr != v1 && curr != static_cast<size_t>(-1)) {
            long int p_long = prev_rev[curr];
            if (p_long == -1 || p_long == INFINITE_INT) break;
            size_t p = static_cast<size_t>(p_long);

            prev[p] = static_cast<long int>(curr);
            curr = p;
        }
    }
}
