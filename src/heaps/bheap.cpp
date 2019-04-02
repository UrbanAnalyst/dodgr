/* File bheap.c - Binary Heap
 * ----------------------------------------------------------------------------
 *  Shane Saunders
 */
#include <cstdlib>
#include "bheap.h"

/* This implementation stores the binary heap in a 1 dimensional array. */


/*--- BHeap (public methods) ------------------------------------------------*/

/* --- Constructor ---
 * Allocates and initialises a binary heap capable of holding n items.
 */
BHeap::BHeap(unsigned int n)
{
    //    int i;

    /* For the purpose of indexing the binary heap, we require n+1 elements in
     * a[] since the indexing method does not use a[0].
     */
    a = new BHeapNode[n+1];
    aPos = new unsigned int[n];
    //    for(i = 0; i <= n; i++) {
    //        a[i].item = 0;
    //	a[i].key = 0;
    //    }
    //    for(i = 0; i < n; i++) aPos[i] = 0;
    itemCount = 0;
    compCount = 0;
}

/* --- Destructor ---
*/
BHeap::~BHeap()
{
    delete [] a;
    delete [] aPos;
}

/* --- min() ---
 * Returns the item with the minimum key in the heap.
 */
unsigned int BHeap::min()
{
    /* the item at the top of the binary heap has the minimum key value */
    return a[1].item;
}

double BHeap::getmin()
{
    return a[1].key;
}

/* --- insert() ---
 * Inserts an item $item$ with associated key value $key$ into the heap.
 */
void BHeap::insert(unsigned int item, double key)
{
    /* i - insertion point
     * j - parent of i
     * y - parent's entry in the heap
     */
    unsigned int i, j;
    BHeapNode y;

    /* $i$ initially indexes the new entry at the bottom of the heap */
    i = ++itemCount;

    /* stop if the insertion point reaches the top of the heap */
    while(i >= 2) {
        /* $j$ indexes the parent of $i$, and $y$ is the parent's entry */
        j = i / 2;
        y = a[j];

        /* We have the correct insertion point when the items key is >= parent
         * Otherwise we move the parent down and insertion point up.
         */
        compCount++;
        if(key >= y.key) break;

        a[i] = y;
        aPos[y.item] = i;
        i = j;
    }

    /* insert the new item at the insertion point found */
    a[i].item = item;
    a[i].key = key;
    aPos[item] = i;
}

/* --- delete() ---
 * Deletes item $item$ from the heap.
 */
void BHeap::deleteItem(unsigned int item)
{
    /* Decrease the number of entries in the heap and record the position of
     * the item to be deleted.
     */
    const unsigned int n = --itemCount;
    const unsigned int p = aPos[item];

    /* Heap needs adjusting if the position of the deleted item was not at the
     * end of the heap.
     */
    if(p <= n) {
        /* We put the item at the end of the heap in the place of the deleted
         * item and sift-up or sift-down to relocate it in the correct place in
         * the heap.
         */
        compCount++;
        if(a[p].key <= a[n+1].key) {
            a[p] = a[n + 1];
            aPos[a[p].item] = p;
            siftUp(p, n);
        }
        else {
            /* Use insert to sift-down, temporarily adjusting the size of the
             * heap for the call to insert.
             */
            itemCount = p - 1;
            insert(a[n+1].item, a[n+1].key);
            itemCount = n;
        }
    }
}

/* --- decreaseKey() ---
 * Decreases the value of $item$'s key to the value $newKey$.
 */
void BHeap::decreaseKey(unsigned int item, double newKey)
{
    const unsigned int n = itemCount;

    itemCount = aPos[item] - 1;
    insert(item, newKey);

    itemCount = n;
}

/*--- BHeap (private methods) -----------------------------------------------*/

/* --- siftUp() ---
 * Considers the sub-tree rooted at index $p$ that ends at index $q$ and moves
 * the root down, sifting up the minimum child until it is located in the
 * correct part of the binary heap.
 */
void BHeap::siftUp(unsigned int p, unsigned int q)
{
    /* y - the heap entry of the root.
     * j - the current insertion point for the root.
     * k - the child of the insertion point.
     * z - heap entry of the child of the insertion point.
     */
    unsigned int j, k;
    BHeapNode y, z;

    /* Get the value of the root and initialise the insertion point and child.
    */
    y = a[p];
    j = p;
    k = 2 * p;

    /* sift-up only if there is a child of the insertion point. */
    while(k <= q) {

        /* Choose the minimum child unless there is only one. */
        z = a[k];
        if(k < q) {
            compCount++;
            if(z.key > a[k + 1].key) z = a[++k];
        }

        /* We stop if the insertion point for the root is in the correct place.
         * Otherwise the child goes up and the root goes down.  (i.e. swap)
         */
        if(y.key <= z.key) break;
        a[j] = z;
        aPos[z.item] = j;
        j = k;
        k = 2 * j;
    }

    /* Insert the root in the correct place in the heap. */
    a[j] = y;
    aPos[y.item] = j;
}

/*--- BHeap (debugging) -----------------------------------------------------*/
void BHeap::dump() const
{

}


/*---------------------------------------------------------------------------*/
