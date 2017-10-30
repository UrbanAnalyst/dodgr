#ifndef TRIHEAP_H
#define TRIHEAP_H
/* File triheap.h - Trinomial Heap
 * ----------------------------------------------------------------------------
 *  Shane Saunders
 */
#include "heap.h"  /* Defines the universal Heap class. */

/* This implementation uses a node-pair structure. */

/*** Option to print debugging information.  Use 1 for yes, or 0 for no. ***/
#define SHOW_trih 0


/* --- TriHeapNode ---
 * Trinomial heap node pairs objects are very similar to Fibonacci heap node
 * objects.  The main difference is that a partner pointer is used which points
 * to a nodes partner in a (2nd, 3rd) node-pair, from the nodes
 * (head, 2nd, 3rd) on a trunk.  Another naming convention is to call the
 * nodes on a trunk (parent, 1st child, 2nd child), so the 1st and 2nd child
 * nodes are paired.
 *
 * We use the boolean variable `extra' to identify whether a node is the extra
 * node (2nd child) in the node-pair.  An extra node only uses the `partner'
 * and `child' pointer, and relies on the left, right and parent pointers of
 * its partner (1st child).  Also, extra nodes are only pointed to by a
 * `partner' pointer of a 1st child node.
 *
 * Note that while a node is an extra node its left, right, and parent pointers
 * are undefined.  When an extra node becomes the main node of the pair, the
 * correct left, right and parent pointers should be assigned to it.
 *
 * The remaining structure fields are:
 * dim       - the nodes dimension.
 * key       - the nodes key.
 * item      - the number of the item that the node is associated with.
 *
 * In this implementation, dimensioning of nodes begins at zero, so the
 * dimension of a single node with no children is zero.
 */
class TriHeapNode {
    public:
        TriHeapNode *parent;
        TriHeapNode *left, *right;
        TriHeapNode *child;
        TriHeapNode *partner;

        unsigned int extra;    
        unsigned int dim;

        double key;
        unsigned int item;
};

/* --- TriHeap ---
 * Trinomial heap class.
 *
 * trees     - An array of pointers to trees at root level in the heap.
 *             Entry i in the array points to the root node of a tree that has
 *             nodes of dimension i on the main trunk.
 * active    - An array of pointers to active nodes in the heap.
 * nodes     - An array of pointers to nodes in the heap.  Nodes are indexed
 *             according to their vertex number.  This array can then be used
 *             to look up the node corresponding to a vertex number, and is
 *             useful when freeing space taken up by the heap.
 * maxNodes  - The maximum number of nodes allowed in the heap.
 * maxTrees  - The maximum number of trees allowed in the heap (calculated from
 *             maxNodes).
 * itemCount - The current number of nodes in the heap.
 * treeSum   - The binary value represented by trees in the heap.
 *             By maintaining this it is easy to keep track of the maximum rank
 *             tree in the heap.
 * compCount - can be used for experimental purposes when counting the number
 *             of key comparisons.
 */
class TriHeap : public Heap {
    public:
        TriHeap(unsigned int maxNodes);
        ~TriHeap();

        unsigned int deleteMin();
        void insert(unsigned int item, double k);
        void decreaseKey(unsigned int item, double newValue);
        unsigned int nItems() const { return itemCount; }

        long int nComps() const { return compCount; }

        void dump() const;

    private:
        TriHeapNode **trees;
        TriHeapNode **active;
        TriHeapNode **nodes;
        unsigned int maxNodes, maxTrees, itemCount, treeSum;
        long int compCount;

        void meld(TriHeapNode *treeList);

        static unsigned int merge(TriHeapNode **a, TriHeapNode **b);
        static void addChild(TriHeapNode *p, TriHeapNode *c);
        static void replaceChild(TriHeapNode *oldNode, TriHeapNode *newNode);

        static void dumpNodes(TriHeapNode *node, unsigned int level);
};

/*---------------------------------------------------------------------------*/
#endif



