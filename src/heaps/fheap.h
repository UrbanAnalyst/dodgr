#ifndef FHEAP_H
#define FHEAP_H
/* File fheap.h - Fibonacci Heap
 * ----------------------------------------------------------------------------
 *  Mark Padgham, adapted from code by Shane Saunders
 */
#include "heap.h"  /* Defines the base class for heaps. */


/* Option to allow printing of debugging information.  Use 1 for yes, or 0 for
 * no.
 */
#define FHEAP_DUMP 0


/* --- FHeapNode ---
 * Fibonacci heap node class.
 *
 * A nodes has the following pointers:
 * parent      - a pointer to the nodes parent node (if any).
 * child       - a pointer to a child node (typically the highest rank child).
 * left, right - sibling pointers which provide a circular doubly linked list
 *               containing all the parents nodes children.
 *
 * The remaining fields are:
 * rank        - the nodes rank, that is, the number of children it has.
 * key         - the nodes key.
 * item        - the number of the item that the node is associated with.
 */
class FHeapNode {
    public:
        FHeapNode *parent;
        FHeapNode *left, *right;
        FHeapNode *child;
        size_t rank;
        size_t marked;
        double key;
        size_t item;
};

/* --- FHeap ---
 * Fibonacci heap class.
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
class FHeap: public Heap {
    public:
        FHeap(size_t n);
        ~FHeap();

        size_t deleteMin();
        void insert(size_t item, double k);
        void decreaseKey(size_t item, double newValue);
        size_t nItems() const { return itemCount; }

        long int nComps() const { return compCount; }
        void dump() const;

        double getmin() {
            return 0.0; // MP: dummy value not implemented yet
        }

    private:
        FHeapNode **trees;
        FHeapNode **nodes;
        size_t maxNodes, maxTrees, itemCount, treeSum;
        long int compCount;

        void meld(FHeapNode *treeList);
        static void dumpNodes(FHeapNode *node, size_t level);
};

/*---------------------------------------------------------------------------*/
#endif
