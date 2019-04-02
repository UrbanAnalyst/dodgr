#ifndef HEAP23_H
#define HEAP23_H
/* File heap23.h - 2-3 Heap
 * ----------------------------------------------------------------------------
 *  Shane Saunders
 */
#include "heap.h"  /* Defines the base class for heaps. */

/* This 2-3 heap implementation uses the same kind of pointer structure for
 * nodes as the Fibonacci heap does.
 *
 * This 2-3 heap implementation was designed for use with Dijkstra's single
 * source shortest path algorithm.  Nodes have a vertex number and a key.
 * This source code could easily be modified to use a pointer (void *) to any
 * structure, rather than using vertex numbers as an index.
 */


/*** Option to print debugging information.  Use 1 for yes, or 0 for no. ***/
#define HEAP23_DUMP 0



/* --- Heap23Node ---
 * Class for node objects in the 2-3 heap.
 *
 * The pointer 'parent' points to a nodes parent, and 'child' points to the
 * highest dimension child in a circular doubly linked list of child nodes.
 * The circular doubly linked list of child nodes is maintained using the
 * sibling pointers 'left' and 'right'.  The parent pointer of root nodes is
 * NULL.
 *
 * The remaining structure fields are:
 *  dim       - the nodes dimension.
 *  key       - the nodes key.
 *  item      - the number of the item that the node is associated with.
 *
 * In this implementation, dimensioning of nodes begins at zero, so the
 * dimension of a single node with no children is zero.
 */
class Heap23Node {
    public:
        Heap23Node *parent;
        Heap23Node *child;
        Heap23Node *left, *right;
        unsigned int dim;
        double key;
        unsigned int item;
};

/* The structure type for a 2-3 heap.
 *
 * trees - An array of pointers to trees at root level in the heap.  Entry i
 *         in the array points to the root node of a tree that has nodes of
 *         dimension i on the main trunk.
 * nodes - An array of pointers to nodes in the heap.  Nodes are indexed
 *         according to their vertex number.  This array can then be used to
 *         look up the node for corresponding to a vertex number, and is
 *         useful when freeing space taken up by the heap.
 * maxNodes - The maximum number of nodes allowed in the heap.
 * maxTrees - The maximum number of trees allowed in the heap (calculated from
 *             maxNodes).
 * itemCount - The current number of nodes in the heap.
 * treeSum - The binary value represented by trees in the heap.
 *           By maintaining this it is easy to keep track of the maximum rank
 *           tree in the heap.
 * compCount - can be used for experimental purposes when counting the number
 *             of key comparisons.
 */
class Heap23 : public Heap {
    public:
        Heap23(unsigned int n);
        ~Heap23();

        void insert(unsigned int item, double k);
        unsigned int deleteMin();
        void decreaseKey(unsigned int item, double newValue);
        unsigned int nItems() const { return itemCount; }

        long int nComps() const { return compCount; }

        void dump() const;

        double getmin() {
            return 0.0; // MP: dummy value not implemented yet
        }

    private:
        Heap23Node **trees;
        Heap23Node **nodes;
        unsigned int maxNodes, maxTrees, itemCount, treeSum;
        long compCount;

        void meld(Heap23Node *tree_list);
        void removeNode(Heap23Node *cutNode);

        static unsigned int merge(Heap23Node **a, Heap23Node **b);
        static void trimExtraNode(Heap23Node *x);
        static void addChild(Heap23Node *p, Heap23Node *c);
        static void replaceNode(Heap23Node *oldNode, Heap23Node *newNode);
        static void swapTrunks(Heap23Node *lowNode, Heap23Node *highNode);    

        static void dumpNodes(Heap23Node *ptr, unsigned int level);
};

/*---------------------------------------------------------------------------*/
#endif
