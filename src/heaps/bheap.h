#ifndef BHEAP_H
#define BHEAP_H
/* File bheap.h - Binary Heap
 * ----------------------------------------------------------------------------
 *  Mark Padgham, adapted from code by Shane Saunders
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
        size_t item;
        double key;
};

class BHeapNodeInt {
    public:
        size_t item;
        int key;
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
        BHeap(size_t n);
        ~BHeap();

        void deleteItem(size_t item);
        size_t deleteMin() {
            size_t v;
            v = min();
            deleteItem(v);
            return v;
        }
        void insert(size_t item, double key);
        void decreaseKey(size_t item, double newKey);
        size_t nItems() const { return itemCount; }

        long int nComps() const { return compCount; }
        void dump() const;

        double getmin();

        /* extra functions */
        size_t min();

    private:
        BHeapNode *a;
        size_t *aPos;
        size_t itemCount;    
        long int compCount;    

        void siftUp(size_t p, size_t q);
};

class BHeapInt : public HeapInt {
    public:
        BHeapInt(size_t n);
        ~BHeapInt();

        void deleteItem(size_t item);
        size_t deleteMin() {
            size_t v;
            v = min();
            deleteItem(v);
            return v;
        }
        void insert(size_t item, int key);
        void decreaseKey(size_t item, int newKey);
        size_t nItems() const { return itemCount; }

        long int nComps() const { return compCount; }
        void dump() const;

        int getmin();

        /* extra functions */
        size_t min();

    private:
        BHeapNodeInt *a;
        size_t *aPos;
        size_t itemCount;    
        long int compCount;    

        void siftUp(size_t p, size_t q);
};

/*---------------------------------------------------------------------------*/
#endif
