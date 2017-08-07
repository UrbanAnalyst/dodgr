#ifndef BHEAP_H
#define BHEAP_H
/* File bheap.h - Binary Heap
 * ----------------------------------------------------------------------------
 *  Shane Saunders
 */
#include "heap.h"  /* Defines the base class for heaps. */

/* This implementation stores the binary heap in a 1 dimensional array. */


/*--- Structure Definitions -------------------------------------------------*/

/* --- BHeapNode ---
 * A binary heap node has the following members:
 *   item - A unique integer identifying the item.  In shortest path algorithms
 *          this is the vertex number it is associated with.
 *   key  - An integer key.  In shortest path algorithms this is the tentative
 *          shortest path distance.
 */
class BHeapNode {
    public:
        int item;
        float key;
};

/* --- BHeap ---
 * Binary heap structure for frontier set in Dijkstra's algorithm.
 *   a[] - an array of binary heap nodes.
 *   p[] - stores the positions of items in the binary heap array a[].
 *   itemCount - the number of items currently in the binary heap.
 *   compCount - the number of key comparison operations
 */
class BHeap : public Heap {
    public:
        BHeap(int n);
        ~BHeap();

        void deleteItem(int item);
        unsigned int deleteMin() {
            int v;
            v = min();
            deleteItem(v);
            return v;
        }
        void insert(int item, float key);
        void decreaseKey(int item, float newKey);
        int nItems() const { return itemCount; }

        long nComps() const { return compCount; }
        void dump() const;

        /* extra functions */
        int min();

    private:
        BHeapNode *a;
        int *aPos;
        int itemCount;    
        long compCount;    

        void siftUp(int p, int q);
};

/*---------------------------------------------------------------------------*/
#endif
