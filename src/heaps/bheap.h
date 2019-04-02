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
        unsigned int item;
        double key;
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
        BHeap(unsigned int n);
        ~BHeap();

        void deleteItem(unsigned int item);
        unsigned int deleteMin() {
            unsigned int v;
            v = min();
            deleteItem(v);
            return v;
        }
        void insert(unsigned int item, double key);
        void decreaseKey(unsigned int item, double newKey);
        unsigned int nItems() const { return itemCount; }

        long int nComps() const { return compCount; }
        void dump() const;

        double getmin();

        /* extra functions */
        unsigned int min();

    private:
        BHeapNode *a;
        unsigned int *aPos;
        unsigned int itemCount;    
        long int compCount;    

        void siftUp(unsigned int p, unsigned int q);
};

/*---------------------------------------------------------------------------*/
#endif
