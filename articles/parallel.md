# Parallel computation

The `dodgr` package implements most calculations by default in parallel,
using the maximum number of available threads or cores. Numbers of
available threads can be determined with either of the following two
functions:

``` r
parallel::detectCores ()
```

    ## [1] 4

``` r
RcppParallel::defaultNumThreads ()
```

    ## [1] 4

Numbers of cores used in calculations can be controlled by specifying
the `numThread` parameter passed to [the `RcppParallel` function
`setThreadOptions`](https://rcppcore.github.io/RcppParallel/#threads_used).
For control over numbers of threads used, this function must be called
prior to calling any `dodgr` functions. For example, single-threaded
processing can be ensured by first making the following call:

``` r
RcppParallel::setThreadOptions (numThreads = 1L)
```
