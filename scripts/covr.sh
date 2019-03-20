#!/usr/bin/env bash

Rscript -e "covr::codecov(line_exclusions = list('src/heaps/bheap.cpp',"\
    "'src/heaps/bheap.h','src/heaps/fheap.h', 'src/heaps/fheap.cpp',"\
    "'src/heaps/heap23.h', 'src/heaps/heap23.cpp', 'src/heaps/heap23_2.h',"\
    "'src/heaps/heap.h', 'src/heaps/heap_lib.h', 'src/heaps/radixheap.h',"\
    "'src/heaps/radixheap.cpp', 'src/heaps/triheap.h', 'src/heaps/triheap.cpp',"\
    "'src/heaps/triheap_ext.h', 'src/heaps/triheap_ext.cpp', 'src/dgraph.cpp'),"\
    "function_exclusions = c('dodgr_streetnet', 'graph_from_pts'))"
