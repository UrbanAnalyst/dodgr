#include <iostream>
#include <cstdlib>
#include <cmath>
#include "radixheap.h"

#include <cstddef> // std::nullptr_t


RadixHeap::RadixHeap(unsigned int n)
{
    itemCount = 0;
    dMin = 0;
    nBuckets = static_cast<unsigned int>(
            ceil(log2(static_cast <double> (MaxKey) + 1.0)) + 2 );

    /* allocate node lookup array (indexed by item no) */
    nodes = new RadixHeapNode *[n];

    for(unsigned int i = 0; i < n; i++) nodes[i] = std::nullptr_t ();

    /* allocate and initialise buckets */
    RadixHeapNode blankNode;
    blankNode.next = blankNode.prev = std::nullptr_t ();
    blankNode.item = -1;
    blankNode.bucket = -1;
    blankNode.key = -1;
    bucketHeaders = new RadixHeapNode[nBuckets + 1];
    for(unsigned int i = 0; i <= nBuckets; i++) {
        bucketHeaders[i] = blankNode;
        bucketHeaders[i].next = &bucketHeaders[i];
        bucketHeaders[i].prev = &bucketHeaders[i];
    }

    /* allocate and initialse upper-limits of buckets */
    u = new int[nBuckets + 1];
    u[0] = -1;
    unsigned int l = 1;
    for(unsigned int i = 1; i <= nBuckets; i++) {
        u[i] = static_cast <int> (l) - 1;
        l *= 2;
    }
    u[nBuckets] = static_cast <int> (n*MaxKey + 1);
}

RadixHeap::~RadixHeap()
{
    delete [] nodes;
    delete [] bucketHeaders;
    delete [] u;
}

void RadixHeap::insert(unsigned int item, double k)
{
    RadixHeapNode *newNode = new RadixHeapNode;
    newNode->item = static_cast <int> (item); // MP explicit conversion for radix
    // MP radix only works for int keys, so k is rounded:
    newNode->key = static_cast <int> (round (k));
    nodes[item] = newNode;
    placeNode(nBuckets,newNode);
    itemCount++;    
#ifdef RADIXHEAP_DEBUG
    Rcpp::Rcout << "performed insert " << item << "(" << k << ")" << std::endl;
    dump();
#endif
}

void RadixHeap::decreaseKey(unsigned int item, double k)
{
    RadixHeapNode *node;
    node = nodes[item];
    removeNode(node);
    // MP radix only works for int keys, so k is rounded:
    node->key = static_cast <int> (round (k));
    placeNode (static_cast <unsigned int> (node->bucket), node);
#ifdef RADIXHEAP_DEBUG
    Rcpp::Rcout << "performed decrease-key (" << k << ") on item " << node->item << std::endl;
    dump();
#endif    
}

unsigned int RadixHeap::deleteMin()
{
    /* if bucket 1 is nonempty, return any of its nodes as the minimum */
    if(bucketHeaders[1].next != &bucketHeaders[1]) {
        RadixHeapNode *minNode = bucketHeaders[1].next;
        removeNode(minNode);
        unsigned int minItem = static_cast <unsigned int> (minNode->item);
        nodes[minItem] = std::nullptr_t ();
        delete minNode;
        itemCount--;        
        return minItem;
    }

    /* find i such that bucket i is the smallest nonempty bucket */
    unsigned int i = 2;
    while(bucketHeaders[i].next == &bucketHeaders[i]) i++;

    /* find and remove the minimum node from bucket i */
    RadixHeapNode *header = &bucketHeaders[i];
    RadixHeapNode *minNode = bucketHeaders[i].next;
    unsigned int minKey = static_cast <unsigned int> (minNode->key);
    RadixHeapNode *node = minNode->next;
    while(node != header) {
        if(static_cast <unsigned int> (node->key) < minKey) {
            minNode = node;
            minKey = static_cast <unsigned int> (node->key);
        }
        node = node->next;
    }
    removeNode(minNode);

    /* recalulate upper bounds on empty buckets */
    u[0] = static_cast <int> (minKey - 1);
    u[1] = static_cast <int> (minKey);
    unsigned int l = 1;
    int s = static_cast <int> (minKey);
    int uMax = u[i];
    for(unsigned int j = 2; j < i; j++) {
        s += l;
        u[j] = s < uMax ? s : uMax;
        l *= 2;
    }

    /* Every vertex in u[i] can now be moved to the empty lower buckets.
     * This is gauranteed since the condition u[i] = u[i-1] must hold.
     */

    /* place nodes from bucket i into lower buckets */
    RadixHeapNode *nextNode = header->next;
    while(nextNode != header) {
        node = nextNode;
        nextNode = nextNode->next;
        placeNode(i-1, node);
    }

    /* bucket i can now be marked as empty */
    bucketHeaders[i].next = bucketHeaders[i].prev = &bucketHeaders[i];

    /* delete the minimum node and return the corresponding item */
#ifdef RADIXHEAP_DEBUG
    Rcpp::Rcout << "performed delete-min " << minNode->item << "("
        << minNode->key << ")" << std::endl;
    dump();
#endif    
    unsigned int minItem = static_cast <unsigned int> (minNode->item);
    nodes[minItem] = std::nullptr_t ();
    delete minNode;
    itemCount--;
    return minItem;
}

void RadixHeap::placeNode(unsigned int startBucket, RadixHeapNode *node)
{
    /* Place the node in the bucket i <= startBucket that corresponds to its
     * key.
     */
    int key = node->key;
    unsigned int i = startBucket;
    do {
        i--;
    } while (u[i] >= key);        
    // MP: RadixHeaps do not work without the following line:
    //if (i < 0)
    //    i = 0;
    insertNode(i+1, node);
}

void RadixHeap::insertNode(unsigned int i, RadixHeapNode *node)
{
    /* link the node into bucket i */
    node->bucket = static_cast <int> (i);
    RadixHeapNode *tailNode = &bucketHeaders[i];
    RadixHeapNode *prevNode = tailNode->prev;    
    node->next = tailNode;
    tailNode->prev = node;
    node->prev = prevNode;
    prevNode->next = node;
}

void RadixHeap::removeNode(RadixHeapNode *node)
{
    /* unlink the node from its bucket */
    node->prev->next = node->next;
    node->next->prev = node->prev;
}

void RadixHeap::dump() const {    
    int i = static_cast <int> (nBuckets);
    while (i > 0 && bucketHeaders[i].next == &bucketHeaders[i])
        i--; 

    do {
        //Rcpp::Rcout << "bucket " << i << "[" << u[i] << "]:  ";
        RadixHeapNode *header = &bucketHeaders[i];
        RadixHeapNode *node = header->next;
        while(node != header) {
            //Rcpp::Rcout << node->item << "(" << node->key << "), ";
            if(node->key > u[i] || node->key <= u[i-1]) {
                throw std::runtime_error ("node in wrong bucket");
            }
            node = node->next;
        }
        //Rcpp::Rcout << std::endl;
        i--;        
    } while(i >= 0);
}
