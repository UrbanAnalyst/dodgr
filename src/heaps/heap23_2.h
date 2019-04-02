#ifndef HEAP23_2_H
#define HEAP23_2_H
/* File heap23.h - 2-3 Heap
 * ----------------------------------------------------------------------------
 *  Shane Saunders
 */

/* This implementation uses a node-pair structure. */
#include "heap.h"  /* Defines the base class for heaps. */

/*** Option to print debugging information.  Use 1 for yes, or 0 for no. ***/
#define HEAP23_DUMP 0

/* --- Heap23Node ---
 * Class for 2-3 heap node-pair objects.  This is very similar to node objects
 * in the Fibonacci heap.  The main difference is that a partner pointer is
 * used which points to a nodes partner in a (2nd, 3rd) node-pair.  We use the
 * boolean variable $extra$ to identify whether a node is the extra node in
 * the node-pair.  An extra node only uses the $partner$ and $child$ pointer,
 * and relies on the left, right and parent pointers of its partner. Also,
 * extra nodes are only pointed to by a $partner$ pointer.
 *
 * Note that while a node is an extra node its left, right, and parent pointers
 * are undefined.  When an extra node becomes the main node of the pair, the
 * correct left, right and parent pointers should be assigned to it.
 *
 * For this description, we refer to the child pointers of a node and the extra
 * node in a node pair as child1 and child2 respectively.
 * Suppose the dimension of the two nodes in a pair is i.  The left pointer of
 * child1, and the child pointers child1 and child2 point to the three
 * (2nd,3rd) node pairs of dimension i-1 in the lower work space.
 *
 * Suppose we have a pointer to one nodes in the work space.  All six
 * dimension i nodes in the workspace can be accessed by following the
 * parent pointer to the dimension i+1 node pair, then using the three pointers
 * `left', child1, and child2.
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
        Heap23Node *left, *right;
        Heap23Node *child;
        Heap23Node *partner;

        unsigned int extra;    
        unsigned int dim;

        double key;
        unsigned int item;
};

/* --- Heap23 ---
 * 2-3 Heap class.
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
 *         By maintaining this it is easy to keep track of the maximum rank
 *         tree in the heap.
 * compCount - can be used for experimental purposes when counting the number
 *             of key comparisons.
 */
class Heap23 : public Heap {
    public:
        Heap23(unsigned int maxNodes);
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
        long int compCount;

        void meld(Heap23Node *treeList);
        void removeNode(Heap23Node *cutNode);

        static unsigned int merge(Heap23Node **a, Heap23Node **b);
        static void removeChild(Heap23Node *c, Heap23Node *p);
        static void addChild(Heap23Node *p, Heap23Node *c);
        static void replaceChild(Heap23Node *oldNode, Heap23Node *newNode);

        static void dumpNodes(Heap23Node *node, unsigned int level);
};

/*---------------------------------------------------------------------------*/
#endif
