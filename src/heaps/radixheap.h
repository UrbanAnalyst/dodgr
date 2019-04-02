#ifndef RADIXHEAP_H
#define RADIXHEAP_H
#include "heap.h"

//#define RADIXHEAP_DEBUG

// MP RadixHeaps need to have int item, key, bucket, rather than the unsigned
// int vals that all other heap types have. This is because they need to be
// initialised to -1 to find the min heap and work on that.

class RadixHeapNode {
    public:    
        int item, key;
        int bucket;
        RadixHeapNode *next, *prev;
};

class RadixHeap: public Heap {
    public:
        RadixHeap(unsigned int n);
        ~RadixHeap();

        unsigned int deleteMin();
        void insert(unsigned int item, double k);
        void decreaseKey(unsigned int item, double newValue);
        unsigned int nItems() const { return itemCount; }

        long int nComps() const { return compCount; }
        void dump() const;

        double getmin() {
            return 0.0; // MP: dummy value not implemented yet
        }

    private:
        void placeNode(unsigned int startBucket, RadixHeapNode *node);
        void insertNode(unsigned int i, RadixHeapNode *node);
        void removeNode(RadixHeapNode *node);

        static const unsigned int MaxKey = 500000;
        RadixHeapNode **nodes;
        RadixHeapNode *bucketHeaders;
        int *u;

        unsigned int nBuckets;
        unsigned int dMin;

        unsigned int itemCount;
        long int compCount;
};

#endif
