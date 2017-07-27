#ifndef RADIXHEAP_H
#define RADIXHEAP_H
#include "heap.h"

//#define RADIXHEAP_DEBUG


class RadixHeapNode {
    public:    
        int item, key;
        int bucket;
        RadixHeapNode *next, *prev;
};

class RadixHeap: public Heap {
    public:
        RadixHeap(int n);
        ~RadixHeap();

        unsigned int deleteMin();
        void insert(int item, float k);
        void decreaseKey(int item, float newValue);
        int nItems() const { return itemCount; }

        long nComps() const { return compCount; }
        void dump() const;

    private:
        void placeNode(int startBucket, RadixHeapNode *node);
        void insertNode(int i, RadixHeapNode *node);
        void removeNode(RadixHeapNode *node);

        static const int MaxKey = 500000;
        RadixHeapNode **nodes;
        RadixHeapNode *bucketHeaders;
        int *u;

        int nBuckets;
        int dMin;

        int itemCount;
        int compCount;
};

#endif
