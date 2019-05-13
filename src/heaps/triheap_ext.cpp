/*** Trinomial Heap Implementation ***/
/*
 * Mark Padgham, adapted from code by Shane Saunders
 */

/* This version is implemented using the node-pair pointer structure; that is,
 * nodes have a partner pointer, so that nodes can be paired.  In this
 * implementation, a nodes child pointer points to its highest dimension child
 * node.
 */

#include "triheap_ext.h"
#include <cstdlib>
#include <cmath>
//#if SHOW_trih_ext
#include <cstdio>
#include <stdexcept>
//#endif


#define TRUE 1
#define FALSE 0

/*--- TriHeapExt (public methods) -------------------------------------------*/

/* --- Constructor ---
 * creates and returns a pointer to a trinomial heap.  Argument
 * maxNodes specifies the maximum number of nodes the heap can contain.
 */
TriHeapExt::TriHeapExt(unsigned int n)
{
#if SHOW_trih_ext
    Rcpp::Rcout << "init, ";
    Rcpp::Rcout.flush ();
#endif

    /* The maximum number of nodes and the maximum number of trees allowed. */
    maxNodes = n;
    maxTrees = 1 + static_cast <unsigned int> (log(static_cast <double> (n)) /
            log (3.0));

    /* The tolerance trigger of the heap.  That is, t+1; the number used for
     * detecting when we are over the tolerance and need to cleanup active
     * nodes.
     */
    activeLimit = maxTrees;

    /* Allocate space for an array of pointers to trees, and nodes in the heap.
     * Initialise all array entries to zero, that is, NULL pointers.
     */
    trees = new TriHeapExtNode *[maxTrees];
    for (unsigned int i = 0; i < maxTrees; i++)
        trees[i] = std::nullptr_t ();

    nodes = new TriHeapExtNode *[n];
    for (unsigned int i = 0; i < n; i++)
        nodes[i] = std::nullptr_t ();

    /* Allocate space for:
     *  - an unordered array of pointers to active nodes.
     *  - an array of pointers to queues of active nodes.
     *  - an array of pointers to entries in the candidates queue.
     */
    activeNodes = new TriHeapExtNode *[activeLimit];
    for (unsigned int i = 0; i < activeLimit; i++)
        activeNodes[i] = std::nullptr_t ();

    activeQueues = new ActiveItem *[activeLimit-1];
    for (unsigned int i = 0; i < activeLimit-1; i++)
        activeQueues[i] = std::nullptr_t ();

    candidateItems = new CandidateItem *[activeLimit-1];
    for (unsigned int i = 0; i < activeLimit-1; i++)
        candidateItems[i] = std::nullptr_t ();

    /* We begin with no nodes in the heap. */
    itemCount = 0;
    activeCount = 0;
    candQueueHead = std::nullptr_t ();

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
TriHeapExt::~TriHeapExt()
{
#if SHOW_trih_ext
    Rcpp::Rcout << "free, ";
    Rcpp::Rcout.flush ();
#endif

    for (unsigned int i = 0; i < maxNodes; i++) {
        delete nodes[i];
    }

    delete [] nodes;
    delete [] trees;
    delete [] activeNodes;
    delete [] activeQueues;
    delete [] candidateItems;
}

/* --- insert() ---
 * inserts item $item$ with associated key $k$ into the heap.
 */
void TriHeapExt::insert(unsigned int item, double k)
{
    TriHeapExtNode *newNode;

#if SHOW_trih_ext
    Rcpp::Rcout << "insert, ";
    Rcpp::Rcout.flush ();
#endif

    /* Create an initialise the new node.  The parent pointer will be set to
     * NULL by meld().
     */
    newNode = new TriHeapExtNode;
    newNode->child = std::nullptr_t ();
    newNode->extra = FALSE;
    newNode->left = newNode->right = std::nullptr_t ();
    newNode->partner = std::nullptr_t ();

    newNode->activeEntry = std::nullptr_t ();

    newNode->dim = 0;
    newNode->item = item;
    newNode->key = k;

    /* Maintain a pointer to item's new node in the heap. */
    nodes[item] = newNode;

    /* Meld the new node into the heap. */
    meld(newNode);

    /* Update the heap's node count. */
    itemCount++;
}

/* --- deleteMin() ---
 * Deletes and returns the minimum item from the heap.
 */
unsigned int TriHeapExt::deleteMin()
{
    static TriHeapExtNode meldListHeader;

    TriHeapExtNode *minNode, *child, *next, *partner;
    TriHeapExtNode *tail, *breakNode, *firstChild;
    TriHeapExtNode *l, *parent, *childZero, *childHigher;
    TriHeapExtNode *nextPartner, *nextParent, *nextFirstChild;
    TriHeapExtNode *node, *nodePartner;
    TriHeapExtNode *nextChildZero = std::nullptr_t (),
                   *nextChildHigher = std::nullptr_t ();
    double k, k2;
    unsigned int d, nextDim, v, item;
    unsigned int wasExtra;

#if SHOW_trih_ext
    dump();
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
     * never get active nodes on a main trunk.
     */
    minNode = trees[d];
    k = minNode->key;
    while(d > 0) {
        d--;
        next = trees[d];
        if(next) {
            compCount++;
            if((k2 = next->key) < k) {
                k = k2;
                minNode = next;
            }
        }

    }
    unsigned int i = activeCount;
    while(i > 0) {
        i--;
        next = activeNodes[i];
        compCount++;
        if((k2 = next->key) < k) {
            k = k2;
            minNode = next;
        }
    }


#if SHOW_trih_ext
    Rcpp::Rcout << "on vertex no " << minNode->item;
    Rcpp::Rcout.flush ();
#endif

    /*
     * The node-pair tree containing the minimum node must be broken up into
     * a linked list of node-pair sub-trees to be melded back into the heap.
     * This linked list first includes any children trees broken of the minimum
     * node.
     *
     */
    tail = &meldListHeader;

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
        tail->right = child->right;
        node = tail = child;

        /* Nodes in this break up may go from active to inactive. */
        do {
            node = node->right;
            nodePartner = node->partner;

            if(node->activeEntry) {
                deactivate(node);
                if(nodePartner->activeEntry) deactivate(nodePartner);
            }
            /* Note that we do not need to check a second child for activeness
             * if the first child was not active.
             */

        } while(node != tail);
    }

    /*
     * The linked list also includes higher trees that become broken up in
     * cases where a non-root (i.e. active) minimum node was removed.
     *
     */

    /* During higher break-up breakNode points to the node that `appears' to
     * have been removed.  Initially this is minNode.
     */
    d = minNode->dim;
    breakNode = minNode;
    partner = breakNode->partner;
    firstChild = breakNode->extra ? partner : breakNode;
    parent = firstChild->parent;  /* parent pointer of the broken node. */

    /* If the minimum node was active then deactivate it, and propagate
     * break-up through the tree until the main trunk is reached.
     * Nodes on the main trunk are treated like first child and second child
     * with no parent.
     */
    if(parent) {

        deactivate(minNode);

        childZero = parent->child->right;
        childHigher = firstChild->right;

        /* If `partner' is active then deactivate it and determine the order of
         * `partner' and `parent' in the trunk resulting from break-up.  Otherwise
         * parent is the first node on the trunk.
         */
        if(partner->activeEntry) {
            deactivate(partner);
            compCount++;
            if(partner->key < parent->key) {
                tail->right = partner;
                tail = partner;
            }
            else {
                tail->right = parent;
                tail = parent;
            }
        }
        else {
            tail->right = parent;
            tail = parent;
        }

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
             * Also remember the value of parent->dim.  The code has been written
             * this way because certain pointer updates overwrite pointers that
             * are needed later on.
             * Note, if we must deactivate parent, this must be done before
             * changing its dimension.
             */
            if(parent->activeEntry) deactivate(parent);
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
            if(nextParent) {
                nextChildZero = nextParent->child->right;
                nextChildHigher = nextFirstChild->right;
            }


            /* Add the child pairs of `parent' that have a greater dimension than
             * firstChild to the list of broken trunks.  Keep a pointer to
             * parent->child->right for use below.
             */
            if(parent->child != firstChild) {
                node = tail;
                tail->right = childHigher;
                tail = parent->child;

                /* Nodes in this break up may go from active to inactive.
                */
                do {
                    node = node->right;
                    nodePartner = node->partner;
                    if(node->activeEntry) {
                        deactivate(node);
                        if(nodePartner->activeEntry) deactivate(nodePartner);
                    }
                } while(node != tail);
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

            /* `partner' may be active, so we may need to swap the order of
             * `partner' and `parent' in the trunk resulting from break-up.
             */
            if(partner->activeEntry) {
                deactivate(partner);
                compCount++;
                if(partner->key < parent->key) {
                    /* We make the linked list point to `partner' instead of
                     * `parent', and make parent an extra node.
                     */
                    tail->right = partner;
                    tail = partner;
                    continue;  /* Back to start of loop. */
                }
            }

            tail->right = parent;
            tail = parent;

        }} /* if-while */

    /* Terminate the list of trees to be melded. */
    tail->right = std::nullptr_t ();

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
    if(meldListHeader.right) meld(meldListHeader.right);


    /* Record the vertex no to return. */
    item = minNode->item;

    /* Delete the old minimum node. */
    nodes[item] = std::nullptr_t ();
    delete minNode;

    return item;
}

/* --- decreaseKey() ---
 * Decrease the key of item $item$ in the heap to the value $newValue$.  It is
 * up to the user to ensure that newValue is in-fact less than or equal to the
 * current value.
 */
void TriHeapExt::decreaseKey(unsigned int item, double newValue)
{
    TriHeapExtNode *v, *v2, *w, *w2, *p, *above, *partner, *activeNode;
    TriHeapExtNode *l, *r, *lowChild, *highChild, *node;
    ActiveItem *activeEntry;
    unsigned int d;

#if SHOW_trih_ext
    dump();
    Rcpp::Rcout << "decreaseKey on vn = " << item << " (" << newValue << ")";
    Rcpp::Rcout.flush ();
#endif


    /* Pointer v points to the decreased node. */
    v = nodes[item];
    v->key = newValue;
    d = v->dim;  /* dimension */

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

    /* If v is already active then rearrangement is not necessary. */
    if(v->activeEntry) {
        /* If v is a second child and its key has became less than its partner,
         * which will also be active, then swap to maintain the correct
         * ordering.
         */
        if(v->extra && v->key < above->key) {
            v->extra = FALSE;
            v2->extra = TRUE;
            replaceChild(v2, v);
        }
        return;
    }

    if(p) {
        /* Non-main trunk. */

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
            }
            else {
                /* If v remains as a second child, then only activate it if
                 * the first child is active.
                 */
                if(v2->activeEntry) activate(v);

                if(activeCount == activeLimit) goto rearrange;  /* see below */

                return;
            }
        }

        /* At this point, v is a first child. */

        activate(v);

        /* If the number of active nodes is at tolerance level we must perform
         * some rearrangement.
         */
        if(activeCount == activeLimit) goto rearrange;  /* see below */

        /* Otherwise it is okay to exit. */
        return;

    }

    /* Otherwise, v is the second node on a main trunk. */

    compCount++;
    if(v->key < above->key) {
        /* If v is smaller, we swap it with the first child (i.e. the root
         * node) to maintain the correct ordering.
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

    /* Get a candidate for rearrangement. */
    d = candQueueHead->dim;
    activeEntry = activeQueues[d];

    activeNode = activeEntry->node;
    if(activeNode->extra) {
        v = activeNode->partner;
        v2 = activeNode;
    }
    else {
        v = activeNode;
        v2 = activeNode->partner;
    }
    p = v->parent;

    /* If we have two active nodes on the same trunk. */
    if(v2->activeEntry) {
        deactivate(v2);

        /* If the 2nd child of these is less than the parent, then
         * do promotion.
         */
        if(v2->key < v->parent->key) goto promote;

        return;
    }

    /* Try the second trunk. */
    activeNode = activeEntry->next->node;
    if(activeNode->extra) {
        w = activeNode->partner;
        w2 = activeNode;
    }
    else {
        w = activeNode;
        w2 = activeNode->partner;
    }

    /* If we have two active nodes on the same trunk. */
    if(w2->activeEntry) {
        deactivate(w2);

        /* If the 2nd child of these is less than the parent, then
         * do promotion.
         */
        if(w2->key < w->parent->key) {
            v = w;  v2 = w2;
            p = w->parent;
            goto promote;
        }

        return;
    }

    /* The final rearrangement always pairs v with w and v2 with w2. */
    v->partner = w;    w->partner = v;
    v2->partner = w2;  w2->partner = v2;

    /* Determine the ordering in the rearrangement. */
    compCount++;
    if(v2->key < w2->key) {
        /* Make the (1st child, 2nd child) pair (v2,w2), replacing (v,v2). */
        v2->extra = FALSE;
        replaceChild(v, v2);

        compCount++;
        if(v->key < w->key) {
            /* Make the pair (v,w), replacing (w,w2). */
            w->extra = TRUE;
            replaceChild(w,v);
            deactivate(w);

            compCount++;
            if(w->key < v->parent->key) {
                /* Both v and w are inconsistent, continue with promotion.
                 * Update v2, and p;
                 */
                v2 = w;
                p = v->parent;
                goto promote;  /* see below */
            }

            /* Otherwise, although w was active, it is not inconsistent.  In
             * this case we do not do promotion.  Instead, only v remains an
             * active node, although v may not be inconsistent.
             * (Avoids a key comparison to check if v is consistent)
             */
            return;
        }
        else {
            /* Make the pair (w,v), replacing (w,w2). */
            v->extra = TRUE;
            deactivate(v);

            compCount++;
            if(v->key < w->parent->key) {
                /* Both v and w are inconsistent, so continue with promotion.
                 * Update v, v2, and p.
                 */
                v2 = v;
                v = w;
                p = w->parent;
                goto promote;  /* see below */
            }

            /* Only w is possibly inconsistent, so it remains an active node.
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
            deactivate(w);

            compCount++;
            if(w->key < v->parent->key) {
                /* Both v and w are inconsistent, so continue with promotion.
                 * Update v2.
                 */
                v2 = w;
                goto promote;  /* see below */
            }

            /* Only v is possibly inconsistent, so it remains an active node.
            */
            return;
        }
        else {
            /* Make the pair (w,v), replacing (v,v2). */
            v->extra = TRUE;
            replaceChild(v,w);
            deactivate(v);

            compCount++;
            if(v->key < w->parent->key) {
                /* Both v and w are inconsistent, so continue with promotion.
                 * Update v, v2.
                 */
                v2 = v;
                v = w;
                goto promote;
            }

            /* Only w is possibly inconsistent, so it remains an active node.
            */
            return;
        }
    }

    /*------------------------------*/
promote:	
    /* Promotion Code.  Reached if the loop has not exited.  Uses variables
     * p, v, v2, and d.  Must ensure that on the trunk (p,v,v2), both v and v2
     * are inconsistent nodes, so that (v,v2,p) will give heap ordering.
     * Node v should still be active, and node v2 should have been deactivated.
     * Node v will later become active at a higher dimension.
     */
    deactivate(v);

    /* First we make v2 a child node of v. */
    v2->extra = FALSE;
    addChild(v, v2);


    /* Then v replaces p.  Any child nodes of p that have a higher dimension
     * than v will become child nodes of v.  Only child nodes of lower
     * dimension than v will be left under p.
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
         * children, and update their parent pointer.  Note that v currently
         * has at least one child, since v2 was made a child of v.
         */
        lowChild = v->right;
        l = v->child;
        r = v->child->right;
        l->right = lowChild;
        lowChild->left = l;
        r->left = highChild;
        highChild->right = r;

        v->child = highChild;

        node = v;
        do {
            node = node->right;
            node->parent = v;
        } while(node != highChild);
    }

    /* partner may be NULL if p is a root node, so don't update
     * partner->partner yet.  See update below.
     */

    /* If p was an extra node no further pointer updates are required. */

    /* Further pointer updates are only needed if p is a first child or a root
     * node.
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

    /* Finally, make p the partner of node v2 (i.e. the 2nd child of node v).
     * Note we cant update the dimension of p yet until we deactivate it.
     */
    v2->partner = p;
    p->partner = v2;

    /* If v is a second child we may need to swap it with the first child to
     * maintain correct ordering.  If v remains as a second child, we do not
     * make it active, unless the first child is also active.
     */
    if(v->extra) {
        v2 = v->partner;

        compCount++;
        if(v->key < v2->key) {
            /* swap */
            v->extra = FALSE;
            v2->extra = TRUE;

            /* Note there is a special case, where a main trunk is reached. */
            if(v2->parent) {
                /* Non-main trunk */
                replaceChild(v2, v);
            }
            else {
                /* Main trunk:  v replaces v2 as the root node. */
                v->extra = FALSE;
                v2->extra = TRUE;
                v->parent = std::nullptr_t ();
                v->left = v->right = v;
                trees[v->dim] = v;

                /* Don't make v active.
                 * If necessary, deactivate p (which v replaced).
                 */
                if(p->activeEntry) deactivate(p);
                p->dim = d;
                return;
            }
        }
        else if(!v2->activeEntry) {
            /* Don't make v active.  This is always the case for main trunks.
             * If necessary, deactivate p (which v replaced).
             */
            if(p->activeEntry) deactivate(p);
            p->dim = d;
            return;
        }
    }
    else {
        /* Don't swap, and if at main trunk level, don't make v active either.
        */
        if(!v->parent) {
            /* If necessary, deactivate p (which v replaced). */
            if(p->activeEntry) deactivate(p);
            p->dim = d;
            return;
        }
    }

    /* The result of the above code never puts active nodes on a main trunk.
     * If at main trunk level, the function will have exited above.
     */

    /* If p was an active node, it no longer is, and v takes its place as an
     * active node.
     */
    if(p->activeEntry) {
        replaceActive(p, v);
    }
    else {
        /* Otherwise make v active. */
        activate(v);
    }
    p->dim = d;
}

/*--- TriHeapExt (private methods) ------------------------------------------*/

/* --- meld() ---
 * Melds the linked list of trees pointed to by $treeList$ into the heap
 * pointed to by $h$.  This function uses the 'right' sibling pointer of nodes
 * to traverse the linked list from lower dimension nodes to higher
 * dimension nodes.  It expects the last nodes `right' pointer to be NULL.
 */
void TriHeapExt::meld(TriHeapExtNode *treeList)
{
    TriHeapExtNode *next, *addTree;
    TriHeapExtNode *carryTree;
    unsigned int d;

#if SHOW_trih_ext
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

#if SHOW_trih_ext
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
}

/* --- activate() ---
 * Add an inactive node to the queue of active nodes.  It is up to the user of
 * this function to ensure that the node was not already active.
 */
void TriHeapExt::activate(TriHeapExtNode *node)
{
    unsigned int d;
    ActiveItem *activeEntry, *first, *last;
    CandidateItem *candidate, *lastC;

    /* Note that we maintain doubly linked lists.  This allows active items
     * and candidates to be removed in O(1) time.  The circular linking allows
     * us to access the last element in the list in O(1) time.
     */

    /* Add n to the array of active nodes. */
    unsigned int i = activeCount++;
    activeNodes[i] = node;

    /* Create an entry for n in the list of pointers to active nodes. */
    activeEntry = new ActiveItem;
    activeEntry->node = node;
    activeEntry->position = i;
    node->activeEntry = activeEntry;

    /* Insertion depends on whether the list is empty or not. */
    d = node->dim;

    if((first = activeQueues[d])) {
        /* At least one item already. */
        last = first->prev;
        last->next = activeEntry;
        activeEntry->prev = last;
        activeEntry->next = first;
        first->prev = activeEntry;

        /* If there was originally one item, but now two, then insert a new
         * entry into the candidates queue.
         */
        if(first == last) {
            candidate = new CandidateItem;
            candidate->dim = d;
            candidateItems[d] = candidate;

            if(candQueueHead) {
                /* At least one candidate is already in the queue. */
                lastC = candQueueHead->prev;
                lastC->next = candidate;
                candidate->prev = lastC;
                candidate->next = candQueueHead;
                candQueueHead->prev = candidate;
            }
            else {
                /* Inserting into empty queue. */
                candQueueHead = candidate;
                candidate->next = candidate->prev = candidate;
            }
        }
    }
    else {
        /* empty */
        activeQueues[d] = activeEntry;
        activeEntry->next = activeEntry->prev = activeEntry;
    }
}

/* --- deactivate() ---
 * Remove an active node from the set of active nodes, and make it inactive.
 */
void TriHeapExt::deactivate(TriHeapExtNode *node)
{
    ActiveItem *activeEntry, *next, *prev, *first, *second;
    TriHeapExtNode *topActive;
    CandidateItem *candidate, *nextCandidate, *prevCandidate;
    unsigned int d;

    /* Obtain pointers to the corresponding entry in the active node structure
     * and remove the node from the array of active nodes.
     */

    activeEntry = node->activeEntry;

    unsigned int i = --activeCount;
    topActive = activeNodes[i];
    activeNodes[activeEntry->position] = topActive;
    topActive->activeEntry->position = activeEntry->position;
    activeNodes[i] = std::nullptr_t ();

    node->activeEntry = std::nullptr_t ();
    d = node->dim;
    first = activeQueues[d];
    second = first->next;

    /* Update the list according to the amount and position of existing list
     * items.  This is a circular doubly linked list.
     */
    if(second != first) {
        /* There are at least two items in the list. */
        prev = activeEntry->prev;
        next = activeEntry->next;

        /* May need to change pointer to first item. */
        if(activeEntry == first) {
            activeQueues[d] = second;
        }

        /* If there were only two nodes, we need to remove the candidate entry
         * from the candidates queue.
         */
        if(second->next == first) {
            /* remove candidate. */
            candidate = candidateItems[d];
            candidateItems[d] = std::nullptr_t ();
            nextCandidate = candidate->next;
            prevCandidate = candidate->prev;

            /* May need to change pointer to the first item. */
            if(nextCandidate == candidate) {
                candQueueHead = std::nullptr_t ();
            }
            else {
                if(candQueueHead == candidate) {
                    candQueueHead = nextCandidate;
                }
                prevCandidate->next = nextCandidate;
                nextCandidate->prev = prevCandidate;
            }

            delete candidate;
        }

        prev->next = next;
        next->prev = prev;
    }
    else {
        /* There is only one item in the list. */
        activeQueues[d] = std::nullptr_t ();
    }

    delete activeEntry;
}

/* --- replaceActive() ---
 * Replaces one active node with another by simply updating the entry to point
 * to the other node.  This is much quicker than performing deactivate() and
 * activate().
 */
void TriHeapExt::replaceActive(
        TriHeapExtNode *oldNode, TriHeapExtNode *newNode)
{
    ActiveItem *activeEntry;

    activeEntry = oldNode->activeEntry;
    newNode->activeEntry = activeEntry;
    oldNode->activeEntry = std::nullptr_t ();
    activeNodes[activeEntry->position] = newNode;

    activeEntry->node = newNode;
}

/*--- TriHeapExt (node manipulation methods) --------------------------------*/

/* --- merge() ---
 * Merges the two trunks pointed to by *a and *b, returning the sum trunk
 * through `a' and any carry tree through `b'.  When this function is used,
 * both parameters `a' and `b' refer to either a 1-node or 2-node trunk.
 *
 * Returns the number of key comparisons used.
 */
unsigned int TriHeapExt::merge(TriHeapExtNode **a, TriHeapExtNode **b)
{
    TriHeapExtNode *tree, *nextTree, *other, *nextOther;
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
 * Adds a new child node pointed to by $c$, to the parent node pointed to by
 * $p$.  The user should ensure that the correct dimension child is being
 * added.
 */
void TriHeapExt::addChild(TriHeapExtNode *p, TriHeapExtNode *c)
{
    TriHeapExtNode *l, *r;

    /* If p already has child nodes we must update the sibling pointers.
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
 * Replaces child node pointed to by $old$ and its sub-tree with the child node
 * pointed to by $new$ and its sub-tree.
 */
void TriHeapExt::replaceChild(
        TriHeapExtNode *oldNode, TriHeapExtNode *newNode)
{
    TriHeapExtNode *parent, *l, *r;

    r = oldNode->right;

    /* If `oldNode' is an only child we only need to initialise the sibling
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

/*--- TriHeapExt (debugging) ------------------------------------------------*/

/* --- dumpNodes() ----
 * Recursively print the nodes of a trinomial heap.
 */
void TriHeapExt::dumpNodes(TriHeapExtNode *node, unsigned int level)
{
    //#if SHOW_trih_ext
    TriHeapExtNode *childNode, *partner;
    unsigned int childCount;

    /* Print leading whitespace for this level. */
#if SHOW_trih_ext
    for (unsigned int i = 0; i < level; i++)
        Rcpp::Rcout << "   ";

    Rcpp::Rcout << node->item << "(" << node->key << ")";
    Rcpp::Rcout.flush ();
    if(node->activeEntry)
        Rcpp::Rcout << '*';
    Rcpp::Rcout << std::endl;
#endif

    if((childNode = node->child)) {
        childNode = node->child->right;

        childCount = 0;

        do {
            dumpNodes(childNode, level+1);
            if(childNode->dim != childCount) {
                throw std::runtime_error ("error(dim)");
            }
            if(childNode->parent != node) {
                throw std::runtime_error ("error(parent)");
            }
            if(childNode->activeEntry == std::nullptr_t () &&
                    childNode->key < node->key) {
                throw std::runtime_error ("error(key)");
            }
            childNode = childNode->right;
            childCount++;
        } while(childNode != node->child->right);

        if(childCount != node->dim) {
            throw std::runtime_error ("error(childCount)");
        }
    }
    else { 
        if(node->dim != 0) {
            throw std::runtime_error ("error(dim)");
        }
    }

    if((partner=node->partner)) {
        if(node->extra==partner->extra) {
            throw std::runtime_error ("error(extra?)");
        }
        if(partner->extra) {
            if(partner->dim != node->dim) {
                throw std::runtime_error ("error(dim)");
            }
            if(partner->activeEntry && ! node->activeEntry) {
                throw std::runtime_error ("error(active)");
            }

            dumpNodes(partner, level);
            if(partner->key < node->key) {
                throw std::runtime_error ("error(key)");
            }
        }
    }
    else if(node->parent) {
        throw std::runtime_error ("error(no partner)");
    }

    if(node->activeEntry) {
        if(node->activeEntry->node != node) {
            throw std::runtime_error ("error(active entry wrong)");
        }
    }
    //#endif
}

/* --- dump() ---
 * Print out a trinomial heap.
 */
void TriHeapExt::dump() const
{
#if SHOW_trih_ext
    TriHeapExtNode *node;
    TriHeapExtNode *firstChild;
    unsigned int c;

    Rcpp::Rcout << std::endl << "value = " << treeSum << std::endl <<
        "array entries 0..maxTrees =";
    for (unsigned int i=0; i<maxTrees; i++) {
        bool is1or0 = trees [i] ? 1 : 0;
        Rcpp::Rcout << " " << is1or0;
    }
    Rcpp::Rcout << std::endl << "active nodes =";
    for (unsigned int i=0; i<activeCount; i++) {
        if((node = activeNodes[i])) {
            Rcpp::Rcout << " " << node->item;
            firstChild = node->extra ? node->partner : node;
            if(!firstChild->parent) {
                throw std::runtime_error ("error(main trunk)");
            }
            if(!node->activeEntry) {
                throw std::runtime_error ("error(inactive)");
            }
            if(node->activeEntry->node != node) {
                throw std::runtime_error ("error(active entry wrong)");
            }
            if(activeNodes[node->activeEntry->position] != node) {
                throw std::runtime_error ("error(active entry index wrong)");
            }

            c = 0;
            for(unsigned int j = i+1; j < activeCount; j++) {
                if(activeNodes[j] == node) c++;
            }
            if(c) {
                throw std::runtime_error ("error(repeated)");
            }
        }
        else {
            throw std::runtime_error ("error(missing active mode)");
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
#endif
}

/*---------------------------------------------------------------------------*/
