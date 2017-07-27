/* Directed Graphs
 * ----------------------------------------------------------------------------
 * Author:  Shane Saunders
 */
#include <cstdio>
#include "dgraph.h"

#include <Rcpp.h>

/*--- DGraph ----------------------------------------------------------------*/

/* --- Constructor ---
 * Creates a DGraph object containing n vertices.
 */
DGraph::DGraph(int n)
{
    nVertices = n;

    vertices = new DGraphVertex[n];
    initVertices();
}

/* --- Destructor ---
*/
DGraph::~DGraph()
{
    clear();
    delete [] vertices;
}

/* --- clear() ---
 * Clears all edges from the graph.
 */
void DGraph::clear()
{
    DGraphEdge *edge, *nextEdge;
    for(int i = 0; i < nVertices; i++) {
        edge = vertices[i].outHead;

        while(edge) {
            nextEdge = edge->nextOut;
            delete edge;
            edge = nextEdge;
        }
    }
    initVertices();
}

void DGraph::initVertices()
{
    for(int i = 0; i < nVertices; i++) {
        vertices[i].outHead = vertices[i].outTail = 0;
        vertices[i].inHead = vertices[i].inTail = 0;
        vertices[i].outSize = vertices[i].inSize = 0;
    }
}

/* --- addNewEdge() ---
 * Adds a new edge from vertex 'source' to vertex 'target' with
 * with a corresponding distance of dist.
 */
void DGraph::addNewEdge(int source, int target, float dist, float wt)
{
    DGraphEdge *newEdge = new DGraphEdge;
    newEdge->source = source;
    newEdge->target = target;
    newEdge->dist = dist;
    newEdge->wt = wt;
    newEdge->nextOut = NULL;
    newEdge->nextIn = NULL;

    DGraphVertex *vertex = &vertices[source];
    if(vertex->outTail) {
        vertex->outTail->nextOut = newEdge;
    }
    else {
        vertex->outHead = newEdge;
    }
    vertex->outTail = newEdge;
    vertex->outSize++;

    vertex = &vertices[target];
    if(vertex->inTail) {
        vertex->inTail->nextIn = newEdge;
    }
    else {
        vertex->inHead = newEdge;
    }
    vertex->inTail = newEdge;
    vertex->inSize++;
}

bool DGraph::edgeExists(int v, int w) const
{
    /* Scan all existing edges from v to determine whether an edge to w exists.
    */
    const DGraphEdge *edge = vertices[v].outHead;
    while(edge) {
        if(edge->target == w) return true;
        edge = edge->nextOut;
    }
    return false;
}

/* --- reachable() ---
 * Test whether all vertices are reachable from the source vertex s.
 */
bool DGraph::reachable (int s) const
{
    int *stack = new int [nVertices];
    int tos = 0;

    int *visited = new int [nVertices];
    for(int i = 0; i < nVertices; i++)
        visited [i] = 0;

    int vertexCount = 0;
    visited [s] = 1;
    stack [tos++] = s;
    DGraphEdge *edge;
    int v, w;
    while (tos) {
        v = stack [--tos];
        vertexCount++;
        edge = vertices [v].outHead;
        while (edge) {
            w = edge->target;
            if (!visited [w]) {
                visited [w] = 1;
                stack [tos++] = w;
            }
            edge = edge->nextOut;
        }
    }

    delete [] stack;
    delete [] visited;

    return vertexCount == nVertices;
}


/* --- print() ---
 * Prints a text representation of the graph to the standard output.
 */
void DGraph::print() const
{
    const DGraphEdge *edge;

    Rcpp::Rcout << "Graph (vertex: edge{dist} list) = " << std::endl;

    for(int i = 0; i < nVertices; i++) {
        Rcpp::Rcout << i << ": ";
        edge = vertices[i].outHead;
        while(edge) {
            Rcpp::Rcout << edge->target << "{" << edge->dist << "} ";
            edge = edge->nextOut;
        }
        Rcpp::Rcout << std::endl;
    }
}
