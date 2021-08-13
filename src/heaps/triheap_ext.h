#ifndef TRIHEAP_EXT_H
#define TRIHEAP_EXT_H
/* File triheap_ext.h - Extended Trinomial Heap
 * ----------------------------------------------------------------------------
 * Mark Padgham, adapted from code by Shane Saunders
 */
#include "heap.h"  /* Defines the universal Heap class. */

/* This implementation uses a node-pair structure. */



/*** Option to print debugging information.  Use 1 for yes, or 0 for no. ***/
#define SHOW_trih_ext 0

/* --- Forward declarations --- */
class ActiveItem;

/* --- TriHeapExtNode ---
 * Trinomial heap node pairs are very similar to Fibonacci heap node pairs.
 * The main difference is that a partner pointer is used which points to a
 * nodes partner in a (2nd, 3rd) node-pair, from the nodes (head, 2nd, 3rd) on
 * a trunk.  Another naming convention is to call the nodes on a trunk
 * (parent, 1st child, 2nd child), so the 1st and 2nd child nodes are paired.
 *
 * We use the boolean variable 'extra' to identify whether a node is the extra
 * node (2nd child) in the node-pair.  An extra node only uses the 'partner'
 * and 'child' pointer, and relies on the left, right and parent pointers of
 * its partner (1st child).  Also, extra nodes are only pointed to by a
 * 'partner' pointer of a 1st child node.
 *
 * Note that while a node is an extra node its left, right, and parent pointers
 * are undefined.  When an extra node becomes the 1st child node of a pair,
 * the correct left, right and parent pointers should be assigned to it.
 *
 * The pointer field activeEntry identifies whether or node the node is
 * active.  If the node is active this field will point to the corresponding
 * ActiveItem object in the linked list of active items.  If the node is not
 * active this field will be NULL.
 *
 * The remaining structure fields are:
 * dim        - the nodes dimension.
 * key        - the nodes key.
 * item       - the number of the graph vertex that the node corresponds to.
 *               Vertex numbering in the graph should be:
 *               1, 2, 3, ... maxVertex.
 *
 * In this implementation, dimensioning of nodes begins at zero, so the
 * dimension of a single node with no children is zero.
 */
class TriHeapExtNode {
    public:
        TriHeapExtNode *parent;
        TriHeapExtNode *left, *right;
        TriHeapExtNode *child;
        TriHeapExtNode *partner;

        size_t extra;
        ActiveItem *activeEntry;
        size_t dim;

        double key;
        size_t item;
};

/* --- CandidateItem ---
 * Class for linked list objects identifying the dimensions for which candidate
 * nodes are available.
 */
class CandidateItem {
    public:
        size_t dim;
        CandidateItem *next, *prev;
};

/* --- ActivePtr ---
 * Class for linked list objects identifying active nodes.
 */
class ActiveItem {
    public:
        TriHeapExtNode *node;
        size_t position;

        ActiveItem *next, *prev;
};

/* --- TriHeapExt ---
 *
 * trees - An array of pointers to trees at root level in the heap.  Entry i
 *         in the array points to the root node of a tree that has nodes of
 *         dimension i on the main trunk.
 * nodes - An array of pointers to nodes in the heap.  Nodes are indexed
 *         according to their vertex number.  This array can then be used to
 *         look up the node corresponding to a vertex number, and is useful
 *         when freeing space taken up by the heap.
 *
 * maxNodes - The maximum number of nodes allowed in the heap.
 * maxTrees - The maximum number of trees allowed in the heap (calculated from
 *             maxNodes).
 * itemCount - The current number of nodes in the heap.
 * treeSum - The binary value represented by trees in the heap.
 *         By maintaining this it is easy to keep track of the maximum rank
 *         tree in the heap.
 * activeLimit - The number of active nodes required to, trigger cleanup.
 * activeCount - The current number of active nodes in the heap.
 *
 * activeNodes - An array of pointers to active nodes in the heap.  No specific
 *          ordering.
 * activeQueues - An array of pointers to queues of active nodes for each
 *               dimension.
 * candidateItems - An array of pointers to entries in the candidate queue.
 * firstCandidate - A pointer to the head of the candidate queue.
 *
 * compCount - can be used for experimental purposes when counting the number
 *             of key comparisons.
 *
 * Note that the queues of candidate and active items are maintained using
 * circularly doubly linked lists.
 */
class TriHeapExt : public Heap {
    public:
        TriHeapExt(size_t n);
        ~TriHeapExt();

        void insert(size_t item, double k);
        size_t deleteMin();
        void decreaseKey(size_t item, double newValue);
        size_t nItems() const { return itemCount; }

        long int nComps() const { return compCount; }    

        void dump() const;

        double getmin() {
            return 0.0; // MP: dummy value not implemented yet
        }

    private:
        TriHeapExtNode **trees;
        TriHeapExtNode **activeNodes;
        TriHeapExtNode **nodes;

        ActiveItem **activeQueues;
        CandidateItem **candidateItems;
        CandidateItem *candQueueHead;

        size_t maxNodes, maxTrees, activeLimit, itemCount, activeCount, treeSum;
        long int compCount;

        void meld(TriHeapExtNode *treeList);
        void activate(TriHeapExtNode *n);
        void deactivate(TriHeapExtNode *n);
        void replaceActive(TriHeapExtNode *oldNode, TriHeapExtNode *newNode);

        static void dumpNodes(TriHeapExtNode *node, size_t level);

        static size_t merge(TriHeapExtNode **a, TriHeapExtNode **b);
        static void addChild(TriHeapExtNode *p, TriHeapExtNode *c);
        static void replaceChild(TriHeapExtNode *oldNode, TriHeapExtNode *newNode);
};

/*---------------------------------------------------------------------------*/
#endif
