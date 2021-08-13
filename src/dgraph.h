#ifndef DGRAPH_H
#define DGRAPH_H

#include <bits/c++config.h>
#include <vector>
#include <limits>

typedef std::size_t size_t;

/* Directed Graphs
 * ----------------------------------------------------------------------------
 * Author: Mark Padgham, modified from code by Shane Saunders
 */


/*--- Directed Graph Classes ------------------------------------------------*/

/* --- Directed graph edge class ---
 * Each edge object represents an outgoing edge of some vertex and is stored in
 * that vertexes linked list of edges.   The member 'target' is the edge's
 * target vertex number, whereas 'source' is the edge's source vertex number.
 * The member 'dist' is the associated edge distance.  The pointers 'nextIn'
 * and 'nextOut' are used to form a linked lists of a vertices incoming and
 * outgoing edges respectively.  Such linked lists are terminated with a null
 * pointer.
 */
class DGraphEdge {
    public:
        size_t source, target, edge_id; // edge_id only used in centrality
        double dist, wt;
        DGraphEdge *nextOut, *nextIn;
};

/* --- Directed graph vertex class ---
 * Each vertex object has an associated linked lists of edge objects
 * representing the outgoing and incoming edges of that vertex.  The member
 * pointers outHead and inHead points to the first edge object in the linked
 * list of outgoing, and incoming edges respectively.  Similarly, outTail and
 * inTail point to the last edge of each linked list.  The number of outgoing
 * and incoming edges are stored in outSize and inSize respectively.
 */
class DGraphVertex {
    public:
        DGraphEdge *outHead, *outTail;
        DGraphEdge *inHead, *inTail;
        int outSize, inSize;
};

/* --- Directed graph class ---
 * Vertices in the graph are stored as an array of vertex objects, pointed to
 * by the member variable 'vertices'.  Each vertex is identified by a number
 * corresponding to its index in the vertices[] array.  The member
 * nVertices is the number of vertices in the graph.
 *
 * clear()      - Remove all edges from graph.
 *
 * addNewEdge() - Adds a new edge to the edge to the graph.
 *
 * print()      - Prints a text representation of the graph to the standard
 *                output.
 */
class DGraph {
    public:
        DGraph(size_t n);
        ~DGraph();
        
        size_t nVertices() const;
        const std::vector<DGraphVertex>& vertices() const;
        
        // disable copy/assign as will crash (double-delete)
        DGraph(const DGraph&) = delete;
        DGraph& operator=(const DGraph&) = delete;
    
        void clear();
        void addNewEdge(size_t srcVertexNo, size_t destVertexNo,
                double dist, double wt, size_t edge_id);
        bool edgeExists(size_t v, size_t w) const;
        bool reachable(size_t s) const;
        void print() const;
    private:
        void initVertices();
    
        std::vector<DGraphVertex> m_vertices;
  
};


/*---------------------------------------------------------------------------*/
#endif
