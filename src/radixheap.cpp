#include <iostream>
#include <cstdlib>
#include <cmath>
#include "radixheap.h"

#include <Rcpp.h>

RadixHeap::RadixHeap(int n)
{
    itemCount = 0;
    dMin = 0;
    nBuckets = static_cast<int>( ceil(log2(MaxKey + 1)) + 2 );

    /* allocate node lookup array (indexed by item no) */
    nodes = new RadixHeapNode *[n];
    for(int i = 0; i < n; i++) nodes[i] = 0;

    /* allocate and initialise buckets */
    RadixHeapNode blankNode;
    blankNode.next = blankNode.prev = 0;
    blankNode.item = -1;
    blankNode.bucket = -1;
    blankNode.key = -1;
    bucketHeaders = new RadixHeapNode[nBuckets + 1];
    for(int i = 0; i <= nBuckets; i++) {
        bucketHeaders[i] = blankNode;
        bucketHeaders[i].next = &bucketHeaders[i];
        bucketHeaders[i].prev = &bucketHeaders[i];
    }

    /* allocate and initialse upper-limits of buckets */
    u = new int[nBuckets + 1];
    u[0] = -1;
    int l = 1;
    for(int i = 1; i <= nBuckets; i++) {
        u[i] = l - 1;
        l *= 2;
    }
    u[nBuckets] = n*MaxKey + 1;
}

RadixHeap::~RadixHeap()
{
    delete [] nodes;
    delete [] bucketHeaders;
    delete [] u;
}

void RadixHeap::insert(int item, float k)
{
    RadixHeapNode *newNode = new RadixHeapNode;
    newNode->item = item;
    newNode->key = k;
    nodes[item] = newNode;
    placeNode(nBuckets,newNode);
    itemCount++;    
#ifdef RADIXHEAP_DEBUG
    Rcpp::Rcout << "performed insert " << item << "(" << k << ")" << std::endl;
    dump();
#endif
}

void RadixHeap::decreaseKey(int item, float k)
{
    RadixHeapNode *node;
    node = nodes[item];
    removeNode(node);
    node->key = k;
    placeNode(node->bucket, node);
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
        int minItem = minNode->item;
        nodes[minItem] = 0;
        delete minNode;
        itemCount--;        
        return minItem;
    }

    /* find i such that bucket i is the smallest nonempty bucket */
    int i = 2;
    while(bucketHeaders[i].next == &bucketHeaders[i]) i++;

    /* find and remove the minimum node from bucket i */
    RadixHeapNode *header = &bucketHeaders[i];
    RadixHeapNode *minNode = bucketHeaders[i].next;
    int minKey = minNode->key;
    RadixHeapNode *node = minNode->next;
    while(node != header) {
        if(node->key < minKey) {
            minNode = node;
            minKey = node->key;
        }
        node = node->next;
    }
    removeNode(minNode);

    /* recalulate upper bounds on empty buckets */
    u[0] = minKey - 1;
    u[1] = minKey;
    int l = 1;
    int s = minKey;
    int uMax = u[i];
    for(int j = 2; j < i; j++) {
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
    int minItem = minNode->item;
    nodes[minItem] = 0;
    delete minNode;
    itemCount--;
    return minItem;
}

void RadixHeap::placeNode(int startBucket, RadixHeapNode *node)
{
    /* Place the node in the bucket i <= startBucket that corresponds to its
     * key.
     */
    int key = node->key;
    int i = startBucket;
    do {
        i--;
    } while (u[i] >= key);        
    // MP: RadixHeaps do not work without the following line:
    if (i < 0)
        i = 0;
    insertNode(i+1, node);
}

void RadixHeap::insertNode(int i, RadixHeapNode *node)
{
    /* link the node into bucket i */
    node->bucket = i;
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
    int i = nBuckets;
    while(i > 0 && bucketHeaders[i].next == &bucketHeaders[i]) i--; 

    do {
        Rcpp::Rcout << "bucket " << i << "[" << u[i] << "]:  ";
        RadixHeapNode *header = &bucketHeaders[i];
        RadixHeapNode *node = header->next;
        while(node != header) {
            Rcpp::Rcout << node->item << "(" << node->key << "), ";
            if(node->key > u[i] || node->key <= u[i-1]) {
                throw std::runtime_error ("node in wrong bucket");
            }
            node = node->next;
        }
        Rcpp::Rcout << std::endl;
        i--;        
    } while(i >= 0);
}
