/* File heap23.c - 2-3 Heap
 * ----------------------------------------------------------------------------
 * Mark Padgham, adapted from code by Shane Saunders
 */

/* This version is implemented using the same pointer structure as a Fibonacci
 * heap; that is, nodes have a parent pointer, and a child pointer which points
 * to a linked list of children constructed from the left and right pointers of
 * nodes.  In this implementation, the child pointer points to the highest
 * dimension child node.
 */
#include <cstdlib>
#include <cmath>
#include <cstdio>
#if defined(TTHEAP_DUMP) && TTHEAP_DUMP > 0
#include <Rcpp.h>
#endif
#include "heap23.h"

/*--- Heap23 (public methods) -----------------------------------------------*/

/* --- Constructor ---
 * Allocates a 2-3 heap capable of holding $n$ items.
 */
Heap23::Heap23(unsigned int n)
{
#if defined(TTHEAP_DUMP) && TTHEAP_DUMP > 0
    Rcpp::Rcout << "init, ";
    Rcpp::Rcout.flush ();
#endif

    /* The maximum number of nodes and the maximum number of trees allowed.
    */
    maxNodes = n; 
    maxTrees = static_cast <unsigned int>(0.5 +
            log(static_cast <double> (n) + 1.0) / log (2.0));

    /* Allocate space for an array of pointers to trees, and nodes in the heap.
     * Initialise all array entries to zero, that is, NULL pointers.
     */
    trees = new Heap23Node *[maxTrees];
    for(unsigned int i = 0; i < maxTrees; i++) trees[i] = std::nullptr_t ();

    nodes = new Heap23Node *[n];
    for(unsigned int i = 0; i < n; i++) nodes[i] = std::nullptr_t ();

    /* We begin with no nodes in the heap. */
    itemCount = 0;

    /* The value of the heap helps to keep track of the maximum rank as nodes
     * are inserted and deleted.
     */
    treeSum = 0;

    /* For experimental purposes, we keep a count of the number of key
     * comparisons.
     */
    compCount = 0;
}

/* --- Destructor ---
*/
Heap23::~Heap23()
{
#if defined(TTHEAP_DUMP) && TTHEAP_DUMP > 0
    Rcpp::Rcout << "free, ";
    Rcpp::Rcout.flush ();
#endif

    for(unsigned int i = 0; i < maxNodes; i++) delete nodes[i];

    delete [] nodes;
    delete [] trees;    
#if defined(TTHEAP_DUMP) && TTHEAP_DUMP > 0
    Rcpp::Rcout << "free-exited, ";
    Rcpp::Rcout.flush ();
#endif
}

/* --- insert() ---
 * Inserts item $item$, with associated key $k$, into the heap.
 */
void Heap23::insert(unsigned int item, double k)
{
    Heap23Node *newNode;

#if defined(TTHEAP_DUMP) && TTHEAP_DUMP > 0
    Rcpp::Rcout << "insert, ";
    dump();
    Rcpp::Rcout.flush ();
#endif

    /* Create an initialise the new node.  The parent pointer will be set to
     * NULL by meld().
     */
    newNode = new Heap23Node;
    newNode->child = std::nullptr_t ();
    newNode->left = newNode->right = std::nullptr_t ();
    newNode->dim = 0;
    newNode->item = item;
    newNode->key = k;

    /* Maintain a pointer to item's new node in the heap. */
    nodes[item] = newNode;

    /* Meld the new node into the heap. */
    meld(newNode);

    /* Update the heap's node count. */
    itemCount++;

#if defined(TTHEAP_DUMP) && TTHEAP_DUMP > 0
    Rcpp::Rcout << "insert-exited, ";
    dump(); 
    Rcpp::Rcout.flush ();
#endif
}

/* --- deleteMin() ---
 * Delete and return the minimum item from the heap.
 */
unsigned int Heap23::deleteMin()
{
    Heap23Node *minNode, *child, *next;
    double k, k2;
    unsigned int r, v, item;

#if defined(TTHEAP_DUMP) && TTHEAP_DUMP > 0
    Rcpp::Rcout << "deleteMin, ";
    Rcpp::Rcout.flush ();
#endif

    /* First we determine the maximum rank tree in the heap. */
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

    /* Now locate the root node with the smallest key, scanning from the
     * maximum rank root position, down to rank 0 root position.
     */
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
    r = minNode->dim;
    trees[r] = std::nullptr_t ();
    treeSum -= (1 << r);
    itemCount--;

    /* A nodes child pointer always points to the child with the highest rank,
     * so child->right is the smallest rank.  For melding the linked list
     * starting at child->right we terminate the circular link with a NULL
     * pointer.
     */
    child = minNode->child;
    if(child) {
        next = child->right;
        next->left = child->right = std::nullptr_t ();
        meld(next);
    }

    /* Record the vertex no to return. */
    item = minNode->item;

    /* Delete the old minimum node. */
    nodes[item] = std::nullptr_t ();
    delete minNode;

#if defined(TTHEAP_DUMP) && TTHEAP_DUMP > 0
    Rcpp::Rcout << "deleteMin-exited, ";
    Rcpp::Rcout.flush ();
#endif

    return item;
}

/* --- decreaseKey() ---
 * Decrease the key of item $item$ in the heap to the value $newValue$.  It
 * is the users reponsibility to ensure that newValue is in-fact less than or
 * equal to the current value.
 */
void Heap23::decreaseKey(unsigned int item, double newValue)
{
    Heap23Node *cutNode, *parent;

#if defined(TTHEAP_DUMP) && TTHEAP_DUMP > 0
    Rcpp::Rcout << "decreaseKey on vn = " << item;
    Rcpp::Rcout.flush ();
#endif

    /* Obtain a pointer to the decreased node and its parent and child.*/
    cutNode = nodes[item];
    parent = cutNode->parent;
    cutNode->key = newValue;

    /* No reinsertion occurs if the node changed was a root. */
    if(!parent) {
#if defined(TTHEAP_DUMP) && TTHEAP_DUMP > 0
        Rcpp::Rcout << "decreaseKey-exited, ";
        Rcpp::Rcout.flush ();
#endif
        return;
    }

    /* Now remove the node and its tree and reinsert it. */
    removeNode(cutNode);
    cutNode->right = cutNode->left = std::nullptr_t ();
    meld(cutNode);

#if defined(TTHEAP_DUMP) && TTHEAP_DUMP > 0
    Rcpp::Rcout << "decreaseKey-exited, ";
    Rcpp::Rcout.flush ();
#endif

}

/*--- Heap23 (private methods) ----------------------------------------------*/

/* --- meld() ---
 * Melds  the linked list of trees pointed to by $treeList$ into the heap
 * pointed to by $h$.  This function uses the $right$ sibling pointer
 * of nodes to traverse the linked list from lower dimension nodes to higher
 * dimension nodes.  It expects the last nodes $right$ pointer to be NULL.
 */
void Heap23::meld(Heap23Node *treeList)
{
    Heap23Node *next, *addTree;
    Heap23Node *carryTree;
    unsigned int d;

#if defined(TTHEAP_DUMP) && TTHEAP_DUMP > 0
    Rcpp::Rcout << "meld - ";
    Rcpp::Rcout.flush ();
#endif

    /* addTree points to the current tree to be merged. */
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

#if defined(TTHEAP_DUMP) && TTHEAP_DUMP > 0
        Rcpp::Rcout << addTree->item << " ";
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


#if defined(TTHEAP_DUMP) && TTHEAP_DUMP > 0
    Rcpp::Rcout << "meld-exited, ";
    Rcpp::Rcout.flush ();
#endif

}

/* --- removeNode() ---
 * Removes the node pointed to by $rNode$, and its corresponding sub-tree, from
 * the heap.  If necessary, this causes rearrangement of $rNode$'s work space.
 */
void Heap23::removeNode(Heap23Node *rNode)
{
    Heap23Node *parent, *child, *ax, *bx, *ap, *bp, *b1, *c, *p;
    unsigned int d, d1;

    parent = rNode->parent;
    child = rNode->child;
    d = rNode->dim;

    /* If this node is an extra node we simply cut the link between it and its
     * parent and update its sibling pointers.
     */
    if(d == parent->dim) {
        trimExtraNode(rNode);
    }
    /* Else if its child is an extra node then use its child to replace it. */
    else if(child && child->dim == d) {

        /* First we remove the child. */
        trimExtraNode(child);

        /* Now we put the child in rNodes position. */
        replaceNode(rNode, child);

    }
    /* Otherwise we need some rearrangement of the workspace. */
    else {

        /* Look at up to two similar nodes in the work space and determine if
         * they have an extra node under them.  Nodes relative to the node
         * being removed are pointed to by the pointers ax, ap, bx, and bp.
         * If a similar trunk lies immediately below cutNode's trunk in the
         * work space, then either ax or ap will be set to point to the node on
         * the end of that trunk.  The same applies for bx and bp, but with a
         * similar trunk immediately above in the work space. the 'x' pointers
         * are set if there is a 3rd (i.e. extra) node on the trunk.  Otherwise
         * the 'p' pointer is set to point to the 2nd node.  Pointers will be
         * set to null if a trunk does not exist or they are not used.
         */

        /* Check for nodes on a similar trunk above in the work space. */
        p = rNode->parent->left;
        if (p->dim == d) {
            c = p->child;
            if(c && c->dim == d) {
                ax = c;  ap = std::nullptr_t ();
            }
            else {
                ap = p;  ax = std::nullptr_t ();
            }
        }
        else {
            ax = ap = std::nullptr_t ();
        }

        /* Check for nodes on a similar trunk below in the work space. */
        d1 = d + 1;
        p = rNode->right;
        if (p->dim == d1) {
            p = p->child;
            if(p->dim == d1) p = p->left;

            c = p->child;
            if(c && c->dim == d) {
                bx = c;  bp = std::nullptr_t ();
            }
            else {
                bp = p;  bx = std::nullptr_t ();
            }
        }
        else {
            bx = bp = std::nullptr_t ();
        }


        if(bx) {

            /* First break `bx's parent link and sibling links. */
            trimExtraNode(bx);

            /* Then we insert bx in rNodes place. */
            replaceNode(rNode, bx);
        }
        else if(bp) {

            b1 = bp->parent;

            /* Recursively remove b1. */
            removeNode(b1);
            b1->dim = d;

            replaceNode(rNode, b1);

            /* It may improve speed by using trimExtraNode() when recursion can be
             * avoided.
             */
        }
        else if(ax) {

            /* Bend the tree to modify its shape then remove rNode. */
            swapTrunks(ax->parent, parent);
            trimExtraNode(rNode);
        }
        else if(ap) {

            /* Bend the tree, so that the node to be relocated, parent, has the
             * larger key value.
             */
            if(parent->key < ap->key) {
                swapTrunks(ap, parent);
                p = parent;
                parent = ap;
                ap = p;
            }
            compCount++;

            trimExtraNode(rNode);
            removeNode(parent);
            parent->dim = d;

            /* Make parent the child of ap. */
            addChild(ap, parent);
        }
        else {
            /* The work space only has rNode node and parent.  This only
             * occurs when parent is a root node, so after removing rNode we
             * demote parent to a lower dimension main trunk.
             */

            /* Note that parent is a root node and has dimension d + 1. */
            trees[d+1] = std::nullptr_t ();
            treeSum -= (1 << (d+1));

            parent->dim = d;
            trimExtraNode(rNode);
            parent->left = parent->right = std::nullptr_t ();

            meld(parent);
        }
    }
}

/*--- Heap23 (node manipulation methods) ------------------------------------*/

/* --- merge() ---
 * Merges the two trunks pointed to by $*a$ and $*b$, returning the sum trunk
 * through $a$ and any carry tree through $b$.  When this function is used,
 * both parameters $a$ and $b$ refer to either a 1-node or 2-node trunk.
 *
 * Returns the number of key comparisons used.
 */
unsigned int Heap23::merge(Heap23Node **a, Heap23Node **b)
{
    Heap23Node *tree, *nextTree, *other, *nextOther;
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
    nextTree = tree->child;
    if(nextTree && nextTree->dim != other->dim)
        nextTree = std::nullptr_t ();
    nextOther = other->child;
    if(nextOther && nextOther->dim != other->dim)
        nextOther = std::nullptr_t ();

    /* The merging depends on the existence of nodes and the values of keys. */
    if(!nextTree) {
        /* nextTree does not exist, so we simply make `other' the child of
         * `tree'. If nextOther exist the resulting 3-node trunk is a carry
         * tree.
         */

        addChild(tree, other);
        if(nextOther) {
            tree->dim++;
            *a = std::nullptr_t ();
            *b = tree;
        }
        else {
            *a = tree;
            *b = std::nullptr_t ();
        }
    }
    else if(!nextOther) {
        /* nextTree exists but nextOther does not, so the linked order of
         * nextTree and `other' in the resulting 3-node trunk depends on the
         * values of keys.  The resulting 3-node trunk becomes a carry tree.
         */

        if(nextTree->key <= other->key) {
            addChild(nextTree, other);
        }
        else {
            replaceNode(nextTree, other);
            addChild(other, nextTree);
        }
        c++;
        tree->dim++;
        *a = std::nullptr_t ();
        *b = tree;
    }
    else {
        /* Otherwise, both nextTree and nextOther exist.  The result consists
         * of a 1 node trunk plus the 3-node trunk which becomes a carry tree.
         * We two trunks are made up as (tree, other, nextOther)
         * and (nextTree).  This uses no key comparisons.
         */

        replaceNode(nextTree, other);
        nextTree->left = nextTree->right = nextTree;
        nextTree->parent = std::nullptr_t ();
        tree->dim++;
        *a = nextTree;  *b = tree;
    }

    return c;
}

/* --- trimExtraNode() ---
 * Trims the extra node pointed to by $x$ from the trunk to
 * which it belongs.
 */
void Heap23::trimExtraNode(Heap23Node *x)
{
    Heap23Node *l, *r;

    if(x->dim == 0) {
        /* A dimension 0 node is an only child, so cutting it leaves no
         * children.
         */

        x->parent->child = std::nullptr_t ();
    }
    else {
        /* Otherwise, sibling pointers of other child nodes must be updated. */

        l = x->left;
        r = x->right;
        l->right = r;
        r->left = l;

        x->parent->child = l;
    }
}

/* --- swapTrunks() ---
 * Where a node in an (i)th trunk, $lowNode$, and a node in an (i+1)th trunk,
 * $highNode$, share the same parent, this function is used for swapping them.
 */
void Heap23::swapTrunks(Heap23Node *lowNode, Heap23Node *highNode)
{
    unsigned int d;
    Heap23Node *parent, *l, *r;

    /* The dimensions of the two nodes are exchanged. */
    d = lowNode->dim;
    lowNode->dim = highNode->dim;
    highNode->dim = d;

    /* Obtain a pointer to the parent of both nodes. */
    parent = highNode->parent;

    /* If the left sibling of lowNode is not highNode, we need to update sibling
     * pointers.  Otherwise, the child pointer of the common parent now
     * points to lowNode.
     */
    if((l = lowNode->left) != highNode) {

        /* Update sibling pointers. */
        r = highNode->right;
        highNode->left = l;
        lowNode->right = r;
        highNode->right = lowNode;
        lowNode->left = highNode;
        l->right = highNode;
        r->left = lowNode;

        /* Determine if the child pointer of the common parent will need to be
         * updated.
         */
        if(parent->child == highNode) {
            parent->child = lowNode;
        }
    }
    else {
        parent->child = lowNode;
    }
}

/* --- addChild() ---
 * Makes node $c$ and its tree a child of node $p$.
 */
void Heap23::addChild(Heap23Node *p, Heap23Node *c)
{
    Heap23Node *l, *r;

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

/* --- replaceNode() ---
 * Replaces node $old$ and its sub-tree with node $new$ and its sub-tree.
 */
void Heap23::replaceNode(Heap23Node *oldNode, Heap23Node *newNode)
{
    Heap23Node *parent, *l, *r;

    l = oldNode->left;
    r = oldNode->right;

    /* If `oldNode' is an only child we only need to initialise the sibling
     * pointers of the new node.  Otherwise we update sibling pointers of other
     * child nodes.
     */
    if(r == oldNode) {
        newNode->right = newNode->left = newNode;
    }
    else {
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

/*--- Heap23 (debugging) ----------------------------------------------------*/

/* --- dumpNopdes() ---
 * Recursively print the nodes of a 2-3 heap.
 */
void Heap23::dumpNodes(Heap23Node *node, unsigned int level)
{
#if defined(TTHEAP_DUMP) && TTHEAP_DUMP > 0
    Heap23Node *childNode, *partner;
    unsigned int childCount;

    /* Print leading whitespace for this level. */
    for(unsigned int i = 0; i < level; i++)
        Rcpp::Rcout << "   ";

    Rcpp::Rcout << node->item << "(" << node->key << ")" << std::endl;

    if((childNode = node->child)) {
        childNode = node->child->right;

        childCount = 0;

        do {
            dumpNodes(childNode, level+1);
            if(childNode->dim != childCount) {
                for(unsigned int i = 0; i < level+1; i++)
                    Rcpp::Rcout << "   ";
                throw std::error ("error(dim)");
            }
            if(childNode->parent != node) {
                for(unsigned int i = 0; i < level+1; i++)
                    Rcpp::Rcout << "   ";
                throw std::error ("error(parent)");
            }
            if(childNode->key < node->key) {
                for(unsigned int i = 0; i < level+1; i++)
                    Rcpp::Rcout << "   ";
                throw std::error ("error(key)");
            }
            childNode = childNode->right;
            childCount++;
        } while(childNode != node->child->right);

        if(childCount != node->dim && childCount != node->dim + 1) {
            for(unsigned int i = 0; i < level; i++)
                Rcpp::Rcout << "   ";
            throw std::error ("error(childCount)");
        }
    }
    else { 
        if(node->dim != 0) {
            for(unsigned int i = 0; i < level; i++)
                Rcpp::Rcout << "   ";
            throw std::error ("error(dim)");
        }
    }
#endif
}

/* --- dump() ---
 * Print a text representation of the heap to the standard output.
 */
void Heap23::dump() const
{
#if defined(TTHEAP_DUMP) && TTHEAP_DUMP > 0
    Heap23Node *node;

    Rcpp::Rcout << std::endl << "value = " << treeSum << std::endl;
    Rcpp::Rcout << "array entries 0..maxTrees =";
    for(unsigned int i=0; i<maxTrees; i++) {
        bool tempval = trees [i] ? 1 : 0;
        Rcpp::Rcout << " " << tempval;
    }
    Rcpp::Rcout << std::endl << std::endl;
    for(unsigned int i=0; i<maxTrees; i++) {
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
