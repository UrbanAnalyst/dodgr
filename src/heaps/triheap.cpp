/* File triheap.h - Trinomial Heap
 * ----------------------------------------------------------------------------
 * Mark Padgham, adapted from code by Shane Saunders
 */

/* This version is implemented using the node-pair pointer structure; that is,
 * nodes have a partner pointer, so that nodes can be paired.  In this
 * implementation, the child pointer points to the highest dimension child
 * node.
 */

#include "triheap.h"
#include <cstdlib>
#include <cmath>
#if SHOW_trih
#include <cstdio>
#include <Rcpp.h>
#endif


#define TRUE 1
#define FALSE 0

/*--- TriHeap ---------------------------------------------------------------*/

/* --- Constructor ---
 * Allocates a trinomial heap capable of holding $n$ items.
 */
TriHeap::TriHeap(unsigned int n)
{
#if SHOW_trih
    Rcpp::Rcout << "init, ";
    Rcpp::Rcout.flush ();
#endif

    /* The maximum number of nodes and the maximum number of trees allowed. */
    maxNodes = n;
    maxTrees = 1 + static_cast <unsigned int> (log(static_cast <double> (n)) /
            log (3.0));

    /* Allocate space for an array of pointers to trees, and nodes in the heap.
     * Initialise all array entries to zero, that is, NULL pointers.
     */
    trees = new TriHeapNode *[maxTrees];
    for (unsigned int i = 0; i < maxTrees; i++)
        trees[i] = std::nullptr_t ();

    active = new TriHeapNode *[maxTrees];
    for (unsigned int i = 0; i < maxTrees; i++)
        active[i] = std::nullptr_t ();

    nodes = new TriHeapNode *[n];
    for (unsigned int i = 0; i < n; i++)
        nodes[i] = std::nullptr_t ();

    /* We begin with no nodes in the heap. */
    itemCount = 0;

    /* The value of the heap helps to keep track of the maximum dimension while
     * nodes are inserted or deleted.
     */
    treeSum = 0;

    /* For experimental purposes, we keep a count of the number of key
     * comparisons.
     */
    compCount = 0;
}

/* --- Destructor ---
*/
TriHeap::~TriHeap()
{
#if SHOW_trih
    Rcpp::Rcout << "free, ";
    Rcpp::Rcout.flush ();
#endif

    for (unsigned int i = 0; i < maxNodes; i++) {
        delete nodes[i];
    }

    delete [] nodes;
    delete [] trees;
    delete [] active;

#if SHOW_trih
    Rcpp::Rcout << "free-exited, ";
    Rcpp::Rcout.flush ();
#endif
}

/* --- insert() ---
 * Inserts a new item, $item$, with assoicated key, $k$, into the heap.
 */
void TriHeap::insert(unsigned int item, double k)
{
    TriHeapNode *newNode;

#if SHOW_trih
    Rcpp::Rcout << "insert, ";
    Rcpp::Rcout.flush ();
#endif

    /* Create an initialise the new node.  The parent pointer will be set to
     * NULL by meld().
     */
    newNode = new TriHeapNode;
    newNode->child = std::nullptr_t ();
    newNode->extra = FALSE;
    newNode->left = newNode->right = std::nullptr_t ();
    newNode->partner = std::nullptr_t ();

    newNode->dim = 0;
    newNode->item = item;
    newNode->key = k;

    /* Maintain a pointer to item's new node in the heap. */
    nodes[item] = newNode;

    /* Meld the new node into the heap. */
    meld(newNode);

    /* Update the heap's node count. */
    itemCount++;

#if SHOW_trih
    Rcpp::Rcout << "insert-exited, ";
    Rcpp::Rcout.flush ();
#endif
}

/* --- deleteMin() ---
 * Deletes and returns the minimum node from the heap.
 */
unsigned int TriHeap::deleteMin()
{
    TriHeapNode *minNode, *child, *next, *partner;
    TriHeapNode *head, *tail, *breakNode, *firstChild;
    TriHeapNode *l, *parent, *childZero, *childHigher;
    TriHeapNode *ptr, *nextPartner, *nextParent, *nextFirstChild;
    TriHeapNode *nextChildZero, *nextChildHigher;
    double k, k2;
    unsigned int d, nextDim, v, item;
    unsigned int wasExtra;

#if SHOW_trih
    Rcpp::Rcout << "deleteMin, ";
    Rcpp::Rcout.flush ();
#endif

    /* First we determine the maximum dimension tree in the heap. */
    v = treeSum;
    d = 0;
    while(v) {
        v = v >> 1;
        d++;
    };
    d--;

    /* Now locate the root node with the smallest key, scanning from the
     * maximum dimension root position, down to dimension 0 root position.
     * At the same time, scan active nodes.  Note that the maximum dimension of
     * a active node is one less than the maximum dimension main trunk since we
     * never get active nodes on a main trunk
     */
    minNode = trees[d];
    k = minNode->key;
    while(d > 0) {
        d--;
        next = trees[d];
        if(next) {
            if((k2 = next->key) < k) {
                k = k2;
                minNode = next;
            }
            compCount++;
        }
        next = active[d];
        if(next) {
            if((k2 = next->key) < k) {
                k = k2;
                minNode = next;
            }
            compCount++;
        }
    }

#if SHOW_trih
    Rcpp::Rcout << "on vertex no " << minNode->item << ", ";
    Rcpp::Rcout.flush ();
#endif

    /* The resulting break-up depends on whether or not minNode is a root node
     * or an active node.
     */

    /* An active node may have been destroyed. */
    if(minNode->parent) {
        active[minNode->dim] = std::nullptr_t ();
    }

    /* During break-up breakNode points to the node that `appears' to have
     * been removed.  Initially this is minNode.
     */
    breakNode = minNode;
    d = breakNode->dim;
    partner = breakNode->partner;
    firstChild = breakNode->extra ? partner : breakNode;
    parent = firstChild->parent;  /* parent pointer of the broken node. */

    /* A nodes child pointer always points to its highest dimension child,
     * so child->right is the smallest dimension.  For melding the linked list
     * starting at child->right, we terminate the circular link with a NULL
     * pointer.
     *
     * For the linked list, only the right pointer is used, so the value of the
     * left pointer will be left undefined.
     */
    child = minNode->child;
    if(child) {
        head = child->right;
        ptr = tail = child;

        /* Nodes in this break up may go from active to inactive. */
        do {
            ptr = ptr->right;
            if(active[ptr->dim] == ptr)
                active[ptr->dim] = std::nullptr_t ();
        } while(ptr != tail);

        if(parent) {
            tail->right = parent;
            tail = parent;
        }
    }
    else {
        head = tail = parent;
    }


    /* Propagate break-up through the tree until the main trunk is reached.
     * Nodes on the main trunk are treated like first child and second child
     * with no parent.  Active nodes never occur as a second child.
     */
    if(parent) {
        childZero = parent->child->right;
        childHigher = firstChild->right;
        while(1) {
            /* At the start of the loop links are not yet updated to account for
             * the node that appears to have been removed.
             * parent - points to the parent of the last broken node.
             * partner - points to the partner of the last broken node.
             * firstChild - is the first child on the trunk of the broken node.
             * The smaller of parent and partner is pointed to by the linked list,
             * Since the linked list is updated in advance.
             */


            /* Make the partner of breakNode `parent's partner.  Remember the
             * old value of the parents partner pointer for the next dimension.
             * Also remember the value of parent->dim.
             */
            nextDim = parent->dim;
            parent->dim = d;
            nextPartner = parent->partner;
            parent->partner = partner;
            partner->partner = parent;

            /* For the last node pair added to the linked list
             * (i.e. (parent,partner) or (partner,parent)), we ensure that the
             * correct node is labelled as extra.
             */
            wasExtra = parent->extra;
            tail->extra = FALSE;
            tail->partner->extra = TRUE;

            /* Obtain future values of pointer variables now.  This is done because
             * constructing the linked list overwrites pointers that are followed.
             */
            if(wasExtra) {
                nextFirstChild = nextPartner;
            }
            else {
                nextFirstChild = parent;
            }
            nextParent = nextFirstChild->parent;
            nextChildZero = std::nullptr_t ();
            nextChildHigher = std::nullptr_t ();
            if(nextParent) {
                nextChildZero = nextParent->child->right;
                nextChildHigher = nextFirstChild->right;
            }


            /* Add the child pairs of `parent' that have a greater dimension than
             * firstChild to the list of broken trunks.  Keep a pointer to
             * parent->child->right for use below.
             */
            if(parent->child != firstChild) {
                ptr = tail;
                tail->right = childHigher;
                tail = parent->child;

                /* Nodes in this break up may go from active to inactive.
                */
                do {
                    ptr = ptr->right;
                    if(active[ptr->dim] == ptr)
                        active[ptr->dim] = std::nullptr_t ();
                } while(ptr != tail);
            }


            /* Update the list of children of `parent' to only include those of
             * lower dimension than firstChild.
             */
            if(d) {
                /* Lower dimension children exist.  Note that tail currently points
                 * to the highest dimension child.
                 */
                l = firstChild->left;
                l->right = childZero;
                childZero->left = l;
                parent->child = l;
            }
            else {
                /* No lower dimension children. */
                parent->child = std::nullptr_t ();
            }


            /* Now continue break up at 1 dimension higher up by treating `parent'
             * as the node that has been broken.
             * Note that on ending, nextParent will be NULL, and we only require
             * `partner' and d to be updated before exiting.
             */
            partner = nextPartner;
            d = nextDim;
            if(!nextParent) break;
            breakNode = parent;
            parent = nextParent;
            firstChild = nextFirstChild;


            /* Pointers to the dimension zero child of parent, and the next highest
             * dimension child from firstChild.
             */
            childZero = nextChildZero; 
            childHigher = nextChildHigher;


            /* Update the linked list in advance to point to the next breakNode,
             * usually `parent', but possibly `partner' if `partner' is active.
             */
            if(wasExtra) {

                /* Since breakNode is not active, `partner' may be, so we may
                 * need to swap the order of `partner' and `parent' in the trunk
                 * resulting from break-up.
                 */
                if(active[d] == partner) {
                    active[d] = std::nullptr_t ();
                    if(partner->key < parent->key) {
                        /* We make the linked list point to `partner' instead of
                         * `parent', and make parent an extra node.
                         */
                        tail->right = partner;
                        tail = partner;
                        continue;  /* Back to start of loop. */
                    }
                    compCount++;
                }
            }
            else {
                if(active[d] == breakNode)
                    active[d] = std::nullptr_t ();
            }
            tail->right = parent;
            tail = parent;

        }} /* if-while */


    /* Break up always propagates up to the main trunk level.  After break up
     * the length the main trunk decreases by one.  The current tree position
     * will become empty unless breakNode has a partner node.
     */
    if(partner) {
        partner->partner = std::nullptr_t ();

        if(partner->extra) {
            partner->extra = FALSE;
            partner->parent = std::nullptr_t ();
            partner->left = partner->right = partner;
            trees[d] = partner;
        }
    }
    else {
        trees[d] = std::nullptr_t ();
        treeSum -= (1 << d);
    }
    itemCount--;

    /* Meld the linked list of trunks resulting from break-up into the main
     * trunk level of the heap.
     */
    if(head) {
        tail->right = std::nullptr_t ();
        meld(head);
    }

    /* Record the vertex no to return. */
    item = minNode->item;

    /* Delete the old minimum node. */
    nodes[item] = std::nullptr_t ();
    delete minNode;

#if SHOW_trih
    Rcpp::Rcout << "deleteMin-exited, ";
    Rcpp::Rcout.flush ();
#endif

    return item;
}

/* --- decreaseKey() ---
 * This mthod decreases the key of the node corresponding to $item$ to
 * $newValue$.   It is up to the user to ensure that $newValue$ is in-fact less
 * than or equal to the current value.
 */
void TriHeap::decreaseKey(unsigned int item, double newValue)
{
    TriHeapNode *v, *v2, *w, *w2, *p, *above, *partner, *activeNode;
    TriHeapNode *l, *r, *lowChild, *highChild, *ptr;
    unsigned int d;

#if SHOW_trih
    Rcpp::Rcout << "decreaseKey on vn = " << item << "(" << newValue << "), ";
    Rcpp::Rcout.flush ();
#endif


    /* Pointer v points to the decreased node. */
    v = nodes[item];
    v->key = newValue;
    d = v->dim;  /* dimension */


    /* This loop allows rearrangement to propagate to higher dimensions. */
    while(1) {


        v2 = v->partner;  /* partner */
        if(v->extra) {
            p = v2->parent;  /* parent */
            above = v2;
        }
        else {
            above = p = v->parent; /* parent */
        }


        /* Determine if rearrangement is necessary */

        /* If v is a root node, then rearrangement is not necessary. */
        if(!above) return;

        if(p) {
            /* Non-main trunk. */

            activeNode = active[d];

            if(v->extra) {
                /* v is a second child.  If necessary, we swap with the first
                 * child to maintain the correct ordering.
                 */

                compCount++;
                if(v->key < above->key) {
                    /* swap */		    
                    v->extra = FALSE;
                    v2->extra = TRUE;
                    replaceChild(v2, v);

                    /* If v2 is inconsistent try promotion.  By checking if
                     * v2 was an active node we can avoid key comparison, since
                     * v2 can only be inconsistent if it is active.
                     */
                    if(activeNode == v2) {

                        compCount++;
                        if(v2->key < p->key) goto promote;  /* see below */

                        /* At this point, v2 is active but is consistent, so v
                         * will replace it as the active node, regardless of
                         * whether v is inconsistent or not.
                         */
                        active[d] = v;
                        return;
                    }

                    /* If there is some other active node for this dimension.
                     * try rearrangement with it.
                     */
                    if(activeNode) goto rearrange;  /* see below */

                    /* Otherwise we can make v an active node, regardless of
                     * whether v is inconsistent or not.
                     */
                    active[d] = v;
                    return;
                }

                /* Otherwise don't swap */

                /* We may need to promote, but only if both v and v2 are
                 * inconsistent.  By checking if v2 is already active, we
                 * can avoid key comparison since v2 won't be inconsistent
                 * if it is not active.
                 */
                if(activeNode == v2) {
                    /* If v is inconsistent, then v2 is also, so promote. */
                    compCount++;
                    if(v->key < p->key) {
                        /* Swap v and v2. */
                        v2 = v;
                        v = above;
                        goto promote;
                    }

                    /* Otherwise promotion is not necessary. */
                    return;
                }

                /* At this point, promotion is not necessary. */
                return;
            }

            /* Otherwise, v is a first child. */

            if(activeNode) {
                /* If v is already the active node leave it. */
                if(activeNode == v) return;

                /* Otherwise some other node is active, so try
                 * rearrangement.
                 */
                goto rearrange;  /* see below */
            }

            /* At this point, there is no active node, so make v the active
             * node, as it is possibly inconsistent.
             */
            active[d] = v;
            return;
        }

        /* Otherwise, v is the second node on a main trunk. */

        compCount++;
        if(v->key < above->key) {
            /* If v is smaller, we swap it with the first child
             * (i.e. the root node) to maintain the correct ordering.
             */
            v->extra = FALSE;
            v2->extra = TRUE;
            v->parent = std::nullptr_t ();
            v->left = v->right = v;
            trees[d] = v;
            return;
        }

        /* Otherwise, no rearrangement is required since heap order is still
         * satisfied.
         */
        return;

        /*------------------------------*/
rearrange:

        /* At this stage v is a first child and possibly inconsistent, but not
         * an active node.  Node v2 is consistent.  Rearrangement occurs with
         * the currently active node for this dimension, which at this point
         * must exist.
         */	


        /* The current active node and its partner are w and w2 respectively.
        */
        w = activeNode;
        w2 = w->partner;

        /* The final rearrangement always pairs v with w and v2 with w2. */
        v->partner = w;    w->partner = v;
        v2->partner = w2;  w2->partner = v2;

        /* Determine the ordering in the rearrangement. */
        compCount++;
        if(v2->key < w2->key) {
            /* Make the (1st child, 2nd child) pair (v2,w2), replacing (v,v2).
            */
            v2->extra = FALSE;
            replaceChild(v, v2);

            compCount++;
            if(v->key < w->key) {
                /* Make the pair (v,w), replacing (w,w2). */
                w->extra = TRUE;
                replaceChild(w,v);

                compCount++;
                if(w->key < v->parent->key) {
                    /* Both v and w are inconsistent, continue with
                     * promotion.  Update v2, and p;
                     */
                    v2 = w;
                    p = v->parent;
                    goto promote;  /* see below */
                }

                /* Otherwise, although w is active, it is not inconsistent.
                 * In this case we do not do promotion.  Instead, v replaces w
                 * as the active node, although v may not be inconsistent.
                 * (Avoids a key comparison to check if v is consistent)
                 */
                active[d] = v;
                return;
            }
            else {
                /* Make the pair (w,v), replacing (w,w2). */
                v->extra = TRUE;

                compCount++;
                if(v->key < w->parent->key) {
                    /* Both v and w are inconsistent, so continue with
                     * promotion.  Update v, v2, and p.
                     */
                    v2 = v;
                    v = w;
                    p = w->parent;
                    goto promote;  /* see below */
                }

                /* Only w is possibly inconsistent, so it remains the active
                 * node.
                 */
                return;
            }
        }
        else {
            /* Make the pair (w2,v2), replacing (w,w2). */
            w2->extra = FALSE;
            replaceChild(w, w2);

            compCount++;
            if(v->key < w->key) {
                /* Make the pair (v,w), replacing (v,v2). */
                w->extra = TRUE;

                compCount++;
                if(w->key < v->parent->key) {
                    /* Both v and w are inconsistent, so continue with
                     * promotion.  Update v2.
                     */
                    v2 = w;
                    goto promote;  /* see below */
                }

                /* Only v is possibly inconsistent, so it becomes the active
                 * node.
                 */
                active[d] = v;
                return;
            }
            else {
                /* Make the pair (w,v), replacing (v,v2). */
                v->extra = TRUE;
                replaceChild(v,w);

                compCount++;
                if(v->key < w->parent->key) {
                    /* Both v and w are inconsistent, so continue with
                     * promotion.  Update v, v2.
                     */
                    v2 = v;
                    v = w;
                    goto promote;
                }

                /* Only w is possibly inconsistent, so it remains the
                 * active node.
                 */
                return;
            }
        }

        /*------------------------------*/
promote:	
        /* Promotion Code.  Reached if the loop has not exited.  Uses variables
         * p, v, v2, and d.  Must ensure that on the trunk (p,v,v2), both v and
         * v2 are inconsistent nodes, so that (v,v2,p) will give heap
         * ordering.
         */

        /* First we make v2 a child node of v. */
        v2->extra = FALSE;
        addChild(v, v2);

        /* Then v replaces p.  Any child nodes of p that have a higher
         * dimension than v will become child nodes of v.  Only child nodes of
         * lower dimension than v will be left under p.
         */

        v->dim = p->dim;
        partner = v->partner = p->partner;
        highChild = p->child;

        if(d) {
            /* v has lower dimension siblings. */
            l = p->child = v->left;
            r = highChild->right;
            l->right = r;
            r->left = l;
        }
        else {
            /* v was an only child. */
            p->child = std::nullptr_t ();
        }

        if(highChild != v) {
            /* v has higher dimension siblings.  Add them to the list of v's
             * children, and update their parent pointer.  Note that v
             * currently has at least one child, since v2 was made a child
             * of v.
             */
            lowChild = v->right;
            l = v->child;
            r = v->child->right;
            l->right = lowChild;
            lowChild->left = l;
            r->left = highChild;
            highChild->right = r;

            v->child = highChild;

            ptr = v;
            do {
                ptr = ptr->right;
                ptr->parent = v;
            } while(ptr != highChild);
        }

        /* partner may be NULL if p is a root node, so don't update
         * partner->partner yet.  See update below.
         */

        /* If p was an extra node no further pointer updates are required. */

        /* Further pointer updates are only needed if p is a first child or a
         * root node.
         */
        if(!p->extra) {
            if(p->parent) {
                /* p is non-root node and a first child. */
                partner->partner = v;
                replaceChild(p, v);
            }
            else {
                /* p is a root node, so update the tree pointer. */
                if(partner) partner->partner = v;
                trees[p->dim] = v;
                v->left = v->right = v;
                v->parent = std::nullptr_t ();
            }

            /* p will become an extra node, see below. */
            p->extra = TRUE;
        }
        else {
            /* If p was an extra node then v becomes an extra node. */
            partner->partner = v;
            v->extra = TRUE;
        }

        /* Finally, make p the partner of node v2 (i.e. the 2nd child of
         * node v).
         */
        p->dim = d;
        v2->partner = p;
        p->partner = v2;

        /* The result of promotion is to release the active node. */
        active[d] = std::nullptr_t ();

        d = v->dim;  /* next dimension, moving up */

        /* If p was an active node, it no longer is, and v takes its place
         * as an active node.
         */
        if(active[d] == p) {
            active[d] = v;
            return;
        }

        /* Rearrangement continues at a higher dimension, since the position of
         * v has moved up, meaning v may be an inconsistent node.
         */
    }

#if SHOW_trih
    Rcpp::Rcout << "decrease_key-exited, ";
    Rcpp::Rcout.flush ();
#endif

}

/*--- TriHeap (private methods) ---------------------------------------------*/

/* --- meld() ---
 * Melds  the linked list of trees pointed to by *treeList into the heap.
 * This function uses the `right' sibling pointer of nodes to traverse the
 * linked list from lower dimension nodes to higher dimension nodes.
 * It expects the last nodes `right' pointer to be NULL.
 */
void TriHeap::meld(TriHeapNode *treeList)
{
    TriHeapNode *next, *addTree;
    TriHeapNode *carryTree;
    unsigned int d;

#if SHOW_trih
    Rcpp::Rcout << "meld - ";
    Rcpp::Rcout.flush ();
#endif

    /* addTree points to the tree to be merged. */
    addTree = treeList;

    carryTree = std::nullptr_t ();

    do {
        /* addTree() gets merged into the heap, and also carryTree if one
         * exists from a previous merge.
         */

        /* Keep a pointer to the next tree and remove sibling and parent links
         * from the current tree.  The dimension of the next tree is always
         * one greater than the dimension of the previous tree, so this merging
         * is like an addition of two ternary numbers.
         *
         * Note that if addTree is NULL and the loop has not exited, then
         * there is only a carryTree to be merged, so treat it like addTree.
         */
        next = std::nullptr_t ();
        if(addTree) {
            next = addTree->right;
            addTree->right = addTree->left = addTree;
            addTree->parent = std::nullptr_t ();
        }
        else {
            addTree = carryTree;
            carryTree = std::nullptr_t ();
        }

#if SHOW_trih
        Rcpp::Rcout << addTree->item << ", ";
        Rcpp::Rcout.flush ();
#endif

        /* First we merge addTree with carryTree, if there is one.  Note that
         * carryTree contains only one node in its main trunk, and addTree
         * has at most two, so the result is at most one 3-node trunk, which is
         * treated as a 1-node main trunk one dimension higher up.
         */
        if(carryTree) {
            compCount += merge(&addTree, &carryTree);
        }

        /* After the merge, if addTree is NULL, then the resulting tree
         * pointed to by carryTree carries to higher entry, so we do not need
         * to merge anything into the existing main trunk.
         * If addTree is not NULL we add it to the existing main trunk.
         */
        if(addTree) {
            d = addTree->dim;
            if(trees[d]) {
                /* Nodes already in this main trunk position, so merge. */
                compCount += merge(&trees[d], &addTree);
                if(!trees[d]) treeSum -= (1 << d);
                carryTree = addTree;
            }
            else {
                /* No nodes in this main trunk position, so use addTree. */
                trees[d] = addTree;
                treeSum += (1 << d);
            }
        }

        /* Obtain a pointer to the next tree to add. */
        addTree = next;

        /* We continue if there is still a node in the list to be merged, or
         * a carry tree remains to be merged.
         */
    } while(addTree || carryTree);


#if SHOW_trih
    Rcpp::Rcout << "meld-exited, ";
    Rcpp::Rcout.flush ();
#endif

}

/*--- TriHeap (node manipulation methods) -----------------------------------*/

/* --- merge() ---
 * Merges the two trunks pointed to by *a and *b, returning the
 * sum trunk through `a' and any carry tree through `b'.
 * When this function is used, both parameters `a' and `b' refer to either
 * a 1-node or 2-node trunk.
 *
 * Returns the number of key comparisons used.
 */
unsigned int TriHeap::merge(TriHeapNode **a, TriHeapNode **b)
{
    TriHeapNode *tree, *nextTree, *other, *nextOther;
    unsigned int c;

    /* Number of comparisons. */
    c = 0;

    /* `tree' always points to the node with the lowest key.
     * To begin with, `tree' points to the smaller head node, and `other'
     * points to the head node of the other trunk.
     */
    if((*a)->key <= (*b)->key) {
        tree = (*a);
        other = (*b);
    }
    else {
        tree = (*b);
        other = (*a);
    }
    c++;

    /* nextTree points to the next node on the trunk that `tree' is the head
     * of (if there is another node).
     * nextOther points to the next node on the trunk that `other' is the head
     * of (if there is another node).
     */
    nextTree = tree->partner;
    nextOther = other->partner;

    /* The merging depends on the existence of nodes and the values of keys. */
    if(!nextTree) {
        /* nextTree does not exist, so we simply make `other' the child of
         * `tree'. If nextOther exist the resulting 3-node trunk is a carry
         * tree.
         */

        if(nextOther) {
            addChild(tree, other);
            tree->dim++;
            *a = std::nullptr_t ();
            *b = tree;
        }
        else {
            tree->partner = other;
            other->partner = tree;
            other->extra = TRUE;

            *a = tree;
            *b = std::nullptr_t ();
        }
    }
    else if(!nextOther) {
        /* nextTree exists but nextOther does not, so the linked order of
         * nextTree and `other' in the resulting 3-node trunk depends on the
         * values of keys.  The resulting 3-node trunk becomes a carry tree.
         */

        tree->partner = std::nullptr_t ();
        other->partner = nextTree;
        nextTree->partner = other;

        if(other->key < nextTree->key) {    
            addChild(tree, other);
        }
        else {
            nextTree->extra = FALSE;
            other->extra = TRUE;	    
            addChild(tree, nextTree);
        }

        tree->dim++;

        c++;
        *a = std::nullptr_t ();
        *b = tree;
    }
    else {
        /* Otherwise, both nextTree and nextOther exist.  The result consists
         * of a 1 node trunk plus the 3-node trunk which becomes a carry tree.
         * We two trunks are made up as (tree, other, nextOther)
         * and (nextTree).  This uses no key comparisons.
         */

        tree->partner = std::nullptr_t ();
        nextTree->partner = std::nullptr_t ();
        nextTree->extra = FALSE;
        nextTree->left = nextTree->right = nextTree;
        nextTree->parent = std::nullptr_t ();

        addChild(tree, other);

        tree->dim++;

        *a = nextTree;  *b = tree;
    }

    return c;
}

/* --- addChild() ---
 * Adds the child node pointeed to by $c$ to the parent node pointed to by $p$.
 * The user should ensure that the correct dimension child is being added.
 */
void TriHeap::addChild(TriHeapNode *p, TriHeapNode *c)
{
    TriHeapNode *l, *r;

    /* If $p$ already has child nodes we must update the sibling pointers.
     * Otherwise only initialise the left and right pointers of the added
     * child.
     */
    if((l = p->child)) {
        r = l->right;
        c->left = l;
        c->right = r;
        r->left = c;
        l->right = c;
    }
    else {
        c->left = c->right = c;
    }

    p->child = c;
    c->parent = p;
}


/* --- replaceChild() ---
 * Replaces child node pointed to by $oldNode$ and its associated sub-tree with the
 * child node pointed to by $newNode$ and its associated sub-tree.
 */
void TriHeap::replaceChild(TriHeapNode *oldNode, TriHeapNode *newNode)
{
    TriHeapNode *parent, *l, *r;

    r = oldNode->right;

    /* If $oldNode$ is an only child we only need to initialise the sibling
     * pointers of the new node.  Otherwise we update sibling pointers of other
     * child nodes.
     */
    if(r == oldNode) {
        newNode->right = newNode->left = newNode;
    }
    else {
        l = oldNode->left;
        l->right = newNode;
        r->left = newNode;
        newNode->left = l;
        newNode->right = r;
    }

    /* Update parent pointer of the new node and possibly the child pointer
     * of the parent node.
     */
    parent = oldNode->parent;
    newNode->parent = parent;
    if(parent->child == oldNode) parent->child = newNode;
}

/*--- TriHeap (debugging) ---------------------------------------------------*/

/* --- dumpNodes() ---
 * Recursively print the nodes of a trinomial heap.
 */
//TriHeapNode **Active;
void TriHeap::dumpNodes(TriHeapNode *node, unsigned int level)
{
#if SHOW_trih
    TriHeapNode *childNode, *partnerNode;
    unsigned int childCount;

    /* Print leading whitespace for this level. */
    for (unsigned int i = 0; i < level; i++)
        Rcpp::Rcout << "   ";

    Rcpp::Rcout << node->item << "(" << node->key << ")" << std::endl;

    if((childNode = node->child)) {
        childNode = node->child->right;

        childCount = 0;

        do {
            dumpNodes(childNode, level+1);
            if(childNode->dim != childCount) {
                for (unsigned int i = 0; i < level+1; i++)
                    Rcpp::Rcout << "   ";
                throw std::runtime_error ("error(dim)");
            }
            if(childNode->parent != node) {
                for (unsigned int i = 0; i < level+1; i++)
                    Rcpp::Rcout << "   ";
                throw std::runtime_error ("error(parent)");
            }
            if(Active[childNode->dim] != childNode &&
                    childNode->key < node->key) {
                for (unsigned int i = 0; i < level; i++)
                    Rcpp::Rcout << "   ";
                throw std::runtime_error ("error(key)");
            }
            childNode = childNode->right;
            childCount++;
        } while(childNode != node->child->right);

        if(childCount != node->dim) {
            for (unsigned int i = 0; i < level; i++)
                Rcpp::Rcout << "   ";
            throw std::runtime_error ("error(childCount)");
        }
    }
    else { 
        if(node->dim != 0) {
            for (unsigned int i = 0; i < level; i++)
                Rcpp::Rcout << "   ";
            throw std::runtime_error ("error(dim)");
        }
    }

    if((partner=node->partner)) {
        if(node->extra==partner->extra) {
            for (unsigned int i = 0; i < level; i++)
                Rcpp::Rcout << "   ";
            Rcpp::Rcout << partner->item;
            throw std::runtime_error (" - error(extra?)");
        }
        if(partner->extra) {
            if(partner->dim != node->dim) {
                for (unsigned int i = 0; i < level; i++)
                    Rcpp::Rcout << "   ";
                Rcpp::Rcout << partner->item;
                throw std::runtime_error (" - error(dim)");
            }

            dumpNodes(partner, level);
            if(partner->key < node->key) {
                for (unsigned int i = 0; i < level; i++)
                    Rcpp::Rcout << "   ";
                throw std::runtime_error ("error(key)");
            }
        }
    }
    else if(node->parent) {
        for (unsigned int i = 0; i < level; i++)
            Rcpp::Rcout << "   ";
        throw std::runtime_error ("error(no partner)");
    }
#endif
}

/* --- dump() ---
 * Print out a trinomial heap.
 */
void TriHeap::dump() const
{
#if SHOW_trih
    TriHeapNode *node;

    Active = active;

    Rcpp::Rcout << std::endl << "value = " << treeSum << std::endl;
    Rcpp::Rcout << "array entries 0..maxTrees =";
    for (unsigned int i=0; i<maxTrees; i++) {
        bool tempval = trees [i] ? 1 : 0
            Rcpp::Rcout << " " << tempval;
    }
    Rcpp::Rcout << std::endl << "active nodes =";
    for (unsigned int i=0; i<maxTrees-1; i++) {
        if(active[i]) {
            Rcpp::Rcout << " " << active [i]->item;
            if(active[i]->dim != i) {
                throw std::runtime_error ("-error(dim)");
            }
            if(active[i]->extra) {
                throw std::runtime_error ("-error(extra)");
            }
            if(!active[i]->parent) {
                throw std::runtime_error ("-error(parent)");
            }
        }
        else {
            Rcpp::Rcout << " _";
        }
    }
    Rcpp::Rcout << std::endl << std::endl;
    for (unsigned int i=0; i<maxTrees; i++) {
        if((node = trees[i])) {
            Rcpp::Rcout << "tree " << i << std::endl << std::endl;
            dumpNodes(node, 0);
            Rcpp::Rcout << std::endl;
        }
    }
    Rcpp::Rcout.flush ();
#endif
}
/*---------------------------------------------------------------------------*/
