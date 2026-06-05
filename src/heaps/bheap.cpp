/* File bheap.c - Binary Heap
 * ----------------------------------------------------------------------------
 *  Mark Padgham, adapted from code by Shane Saunders
 */
#include "bheap.h"
#include <algorithm> // std::copy, std::min

/* This implementation stores the binary heap in a 1 dimensional array. */


/*--- BHeap (public methods) ------------------------------------------------*/

/* --- Constructor ---
 * Allocates and initialises a binary heap capable of holding n items.
 */
BHeap::BHeap(size_t n)
{
    m_maxn = n;
    /* a[] is 1-indexed so needs itemCount+1 slots. Start small so the active
     * portion stays in cache; grow() doubles capacity as needed up to n+1.
     * aPos[] must remain full-size for O(1) decreaseKey position lookup. */
    m_capacity = (n < 1022) ? n + 2 : 1024;
    a = new BHeapNode[m_capacity];
    aPos = new size_t[n];
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

/* --- grow() ---
 * Doubles the capacity of a[] up to the theoretical maximum (m_maxn + 2).
 */
void BHeap::grow()
{
    size_t new_cap = std::min(m_capacity * 2, m_maxn + 2);
    BHeapNode *new_a = new BHeapNode[new_cap];
    std::copy(a, a + m_capacity, new_a);
    delete [] a;
    a = new_a;
    m_capacity = new_cap;
}

/* --- min() ---
 * Returns the item with the minimum key in the heap.
 */
size_t BHeap::min()
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
void BHeap::insert(size_t item, double key)
{
    /* i - insertion point
     * j - parent of i
     * y - parent's entry in the heap
     */
    size_t i, j;
    BHeapNode y;

    /* $i$ initially indexes the new entry at the bottom of the heap */
    i = ++itemCount;
    if (i >= m_capacity) {
        grow();
    }

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
void BHeap::deleteItem(size_t item)
{
    /* Decrease the number of entries in the heap and record the position of
     * the item to be deleted.
     */
    const size_t n = --itemCount;
    const size_t p = aPos[item];

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
void BHeap::decreaseKey(size_t item, double newKey)
{
    const size_t n = itemCount;

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
void BHeap::siftUp(size_t p, size_t q)
{
    /* y - the heap entry of the root.
     * j - the current insertion point for the root.
     * k - the child of the insertion point.
     * z - heap entry of the child of the insertion point.
     */
    size_t j, k;
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
