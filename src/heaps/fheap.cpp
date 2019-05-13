/* File fheap.c - Fibonacci Heap
 * ----------------------------------------------------------------------------
 * Mark Padgham, adapted from code by Shane Saunders
 */
#include <cstdlib>
#include <cmath>
#if defined(FHEAP_DUMP) && FHEAP_DUMP > 0
#include <cstdio>
#endif
#include "fheap.h"

/*--- FHeap (public methods) ------------------------------------------------*/

/* --- Constructor ---
 * Creates an FHeap object capable of holding up to $n$ items.
 */
FHeap::FHeap(unsigned int n)
{
#if FHEAP_DUMP
    Rcpp::Rcout << "new, ";
#endif
    maxTrees = 1 + static_cast <unsigned int>(1.44 *
            log(static_cast <double> (n)) / log (2.0));
    maxNodes = n;

    trees = new FHeapNode *[maxTrees];
    for(unsigned int i = 0; i < maxTrees; i++) trees[i] = std::nullptr_t ();

    nodes = new FHeapNode *[n];
    for(unsigned int i = 0; i < n; i++) nodes[i] =  std::nullptr_t ();

    itemCount = 0;

    /* The $treeSum$ of the heap helps to keep track of the maximum rank while
     * nodes are inserted or deleted.
     */
    treeSum = 0;

    /* For experimental purposes, we keep a count of the number of key
     * comparisons.
     */
    compCount = 0;

#if FHEAP_DUMP
    Rcpp::Rcout << "new-exited, ";
#endif
}

/* --- Destructor ---
*/
FHeap::~FHeap()
{
#if FHEAP_DUMP
    Rcpp::Rcout << "delete, ";
#endif

    for(unsigned int i = 0; i < maxNodes; i++) delete nodes[i];
    delete [] nodes;
    delete [] trees;

#if FHEAP_DUMP
    Rcpp::Rcout << "delete-exited, ";
#endif
}

/* --- insert() ---
 * Inserts an item $item$ with associated key $k$ into the heap.
 */
void FHeap::insert(unsigned int item, double k)
{
    FHeapNode *newNode;

#if FHEAP_DUMP
    Rcpp::Rcout << "insert, ";
#endif

    /* create an initialise the new node */
    newNode = new FHeapNode;
    newNode->child = std::nullptr_t ();
    newNode->left = newNode->right = newNode;
    newNode->rank = 0;
    newNode->item = item;
    newNode->key = k;

    /* maintain a pointer to $item$'s new node in the heap */
    nodes[item] = newNode;

    /* meld the new node into the heap */
    meld(newNode);

    /* update the heaps node count */
    itemCount++;

#if FHEAP_DUMP
    Rcpp::Rcout << "insert-exited, ";
#endif
}

/* --- deleteMin() ---
 * Deletes and returns the minimum item from the heap.
 */
unsigned int FHeap::deleteMin()
{
    FHeapNode *minNode, *child, *next;
    double k, k2;
    unsigned int r, v, item;

#if FHEAP_DUMP
    Rcpp::Rcout << "deleteMin, ";
#endif

    /* First we determine the maximum rank in the heap. */
    v = treeSum;
    // MP Note: r was an int formerly intialised at -1, but it's used as a
    // direct array index below, so really should be unsigned! The r-- at end
    // fixes that.
    r = 0;
    while(v) {
        v = v >> 1;
        r++;
    };
    r--;

    /* Now determine which root node is the minimum. */
    minNode = trees[r];
    k = minNode->key;
    while(r > 0) {
        r--;
        next = trees[r];
        if(next) {
            if((k2 = next->key) < k) {
                k = k2;
                minNode = next;
            }
            compCount++;
        }
    }

    /* We remove the minimum node from the heap but keep a pointer to it. */
    r = minNode->rank;
    trees[r] = std::nullptr_t ();
    treeSum -= (1 << r);

    child = minNode->child;
    if(child) meld(child);

    /* Record the vertex no of the old minimum node before deleting it. */
    item = minNode->item;
    nodes[item] = std::nullptr_t ();
    delete minNode;
    itemCount--;

#if FHEAP_DUMP
    Rcpp::Rcout << "deleteMin-exited, ";
#endif

    return item;
}

/* --- decreaseKey() ---
 * Decreases the key used for item $item$ to the value newValue.  It is left
 * for the user to ensure that newValue is in-fact less than the current value
 */
void FHeap::decreaseKey(unsigned int item, double newValue)
{
    FHeapNode *cutNode, *parent, *newRoots, *r, *l;
    unsigned int prevRank;

#if FHEAP_DUMP
    Rcpp::Rcout << "decreaseKey on vn = " << item << ", ";
#endif

    /* Obtain a pointer to the decreased node and its parent then decrease the
     * nodes key.
     */
    cutNode = nodes[item];
    parent = cutNode->parent;
    cutNode->key = newValue;

    /* No reinsertion occurs if the node changed was a root. */
    if(!parent) {
#if FHEAP_DUMP
        Rcpp::Rcout << "decreaseKey-exited, ";
#endif
        return;
    }

    /* Update the left and right pointers of cutNode and its two neighbouring
     * nodes.
     */
    l = cutNode->left;
    r = cutNode->right;
    l->right = r;
    r->left = l;
    cutNode->left = cutNode->right = cutNode;

    /* Initially the list of new roots contains only one node. */
    newRoots = cutNode;

    /* While there is a parent node that is marked a cascading cut occurs. */
    while(parent && parent->marked) {

        /* Decrease the rank of cutNode's parent and update its child pointer.
        */
        parent->rank--;
        if(parent->rank) {
            if(parent->child == cutNode) parent->child = r;
        }
        else {
            parent->child = std::nullptr_t ();
        }

        /* Update the cutNode and parent pointers to the parent. */
        cutNode = parent;
        parent = cutNode->parent;

        /* Update the left and right pointers of cutNodes two neighbouring
         * nodes.
         */
        l = cutNode->left;
        r = cutNode->right;
        l->right = r;
        r->left = l;

        /* Add cutNode to the list of nodes to be reinserted as new roots. */
        l = newRoots->left;
        newRoots->left = l->right = cutNode;
        cutNode->left = l;
        cutNode->right = newRoots;
        newRoots = cutNode;
    }

    /* If the root node is being relocated then update the trees[] array.
     * Otherwise mark the parent of the last node cut.
     */
    if(!parent) {
        prevRank = cutNode->rank + 1;
        trees[prevRank] = std::nullptr_t ();
        treeSum -= (1 << prevRank);
    }
    else {
        /* Decrease the rank of cutNode's parent an update its child pointer.
        */
        parent->rank--;
        if(parent->rank) {
            if(parent->child == cutNode) parent->child = r;
        }
        else {
            parent->child = std::nullptr_t ();
        }

        parent->marked = 1;
    }

    /* Meld the new roots into the heap. */
    meld(newRoots);

#if FHEAP_DUMP
    Rcpp::Rcout << "decreaseKey-exited, ";
#endif
}

/*--- FHeap (private methods) -----------------------------------------------*/

/* --- meld() ---
 * melds the linked list of trees pointed to by $treeList$ into the heap.
 */
void FHeap::meld(FHeapNode *treeList)
{
    FHeapNode *first, *next, *nodePtr, *newRoot, *temp, *temp2, *lc, *rc;
    unsigned int r;

#if FHEAP_DUMP
    Rcpp::Rcout << "meld: ";
#endif

    /* We meld each tree in the circularly linked list back into the root level
     * of the heap.  Each node in the linked list is the root node of a tree.
     * The circularly linked list uses the sibling pointers of nodes.  This
     *  makes melding of the child nodes from a deleteMin operation simple.
     */
    nodePtr = first = treeList;

    do {

#if FHEAP_DUMP
        Rcpp::Rcout << nodePtr->item << ", ";
#endif

        /* Keep a pointer to the next node and remove sibling and parent links
         * from the current node.  nodePtr points to the current node.
         */
        next = nodePtr->right;
        nodePtr->right = nodePtr->left = nodePtr;
        nodePtr->parent = std::nullptr_t ();

        /* We merge the current node, nodePtr, by inserting it into the
         * root level of the heap.
         */
        newRoot = nodePtr;
        r = nodePtr->rank;

        /* This loop inserts the new root into the heap, possibly restructuring
         * the heap to ensure that only one tree for each degree exists.
         */
        do {

            /* Check if there is already a tree of degree r in the heap.
             * If there is then we need to link it with newRoot so it will be
             * reinserted into a new place in the heap.
             */
            if((temp = trees[r])) {

                /* temp will be linked to newRoot and relocated so we no
                 * longer will have a tree of degree r.
                 */
                trees[r] = std::nullptr_t ();
                treeSum -= (1 << r);

                /* Swap temp and newRoot if necessary so that newRoot always
                 * points to the root node which has the smaller key of the
                 * two.
                 */
                if(temp->key < newRoot->key) {
                    temp2 = newRoot;
                    newRoot = temp;
                    temp = temp2;
                }
                compCount++;

                /* Link temp with newRoot, making sure that sibling pointers
                 * get updated if rank is greater than 0.  Also, increase r for
                 * the next pass through the loop since the rank of new has
                 * increased.
                 */
                if(r++ > 0) {
                    rc = newRoot->child;
                    lc = rc->left;
                    temp->left = lc;
                    temp->right = rc;
                    lc->right = rc->left = temp;
                }
                newRoot->child = temp;
                newRoot->rank = r;
                temp->parent = newRoot;
                temp->marked = 0;
            }
            /* Otherwise if there is not a tree of degree r in the heap we
             * allow newRoot, which possibly carries moved trees in the heap,
             * to be a tree of degree r in the heap.
             */
            else {

                trees[r] = newRoot;
                treeSum += (1 << r);;

                /* NOTE:  Because newRoot is now a root we ensure it is
                 *        marked.
                 */
                newRoot->marked = 1;
            }

            /* Note that temp will be NULL if and only if there was not a tree
             * of degree r.
             */
        } while(temp);

        nodePtr = next;

    } while(nodePtr != first);

#if FHEAP_DUMP
    Rcpp::Rcout << "meld-exited, ";
#endif
}

/*--- FHeap (debugging) -----------------------------------------------------*/

/* --- dumpNodes() ---
 * Recursively print a text representation of the node subtree starting at the
 * node poined to by $node$, using an indentationlevel of $level$.
 */
void FHeap::dumpNodes(FHeapNode *node, unsigned int level)
{
#if FHEAP_DUMP
    FHeapNode *childNode, *partner;
    unsigned int childCount;

    /* Print leading whitespace for this level. */
    for(unsigned int i = 0; i < level; i++)
        Rcpp::Rcout << "   ";

    Rcpp::Rcout << node->item << "(" << node->key << ")[" << node->rank <<
        "]" << std::endl;

    if((childNode = node->child)) {
        childNode = node->child->right;

        childCount = 0;

        do {
            dumpNodes(childNode, level+1);
            if(childNode->dim > node->dim) {
                for(unsigned int i = 0; i < level+1; i++) Rcpp::Rcout << "   ";
                throw std::runtime_error ("error(dim)");
            }
            if(childNode->parent != node) {
                for(unsigned int i = 0; i < level+1; i++) Rcpp::Rcout << "   ";
                throw std::runtime_error ("error(parent)");
            }
            childNode = childNode->right;
            childCount++;
        } while(childNode != node->child->right);

        if(childCount != node->dim) {
            for(unsigned int i = 0; i < level; i++) Rcpp::Rcout << "   ";
            throw std::runtime_error ("error(childCount)");
        }
    }
    else { 
        if(node->dim != 0) {
            for(unsigned int i = 0; i < level; i++) Rcpp::Rcout << "   ";
            throw std::runtime_error ("error(dim)");
        }
    }
#endif
}

/* --- dump() ---
 * Print a text representation of the heap to the standard output.
 */
void FHeap::dump() const
{
#if FHEAP_DUMP
    FHeapNode *node;

    Rcpp::Rcout << std::endl;
    Rcpp::Rcout << "treeSum = " << treeSum << std::endl;
    Rcpp::Rcout << "array entries 0..maxTrees =";
    for(unsigned int i = 0; i < maxTrees; i++) {
        bool tempval = trees [i] ? 1 : 0;
        Rcpp::Rcout << " " << tempval;
    }
    Rcpp::Rcout << std::endl << std::endl;
    for(unsigned int i = 0; i < maxTrees; i++) {
        if((node = trees[i])) {
            Rcpp::Rcout << "tree " << i;
            dumpNodes(node, 0);
            Rcpp::Rcout << std::endl;
        }
    }
    Rcpp::Rcout.flush ();
#endif    
}


/*---------------------------------------------------------------------------*/
