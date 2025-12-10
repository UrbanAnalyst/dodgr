# Save a weighted streetnet to a local file

The
[weight_streetnet](https://UrbanAnalyst.github.io/dodgr/reference/weight_streetnet.md)
function returns a single `data.frame` object, the processing of which
also relies on a couple of cached lookup-tables to match edges in the
`data.frame` to objects in the original input data. It automatically
calculates and caches a contracted version of the same graph, to enable
rapid conversion between contracted and uncontracted forms. This
function saves all of these items in a single `.Rds` file, so that a the
result of a
[weight_streetnet](https://UrbanAnalyst.github.io/dodgr/reference/weight_streetnet.md)
call can be rapidly loaded into a workspace in subsequent sessions,
rather than re-calculating the entire weighted network.

## Usage

``` r
dodgr_save_streetnet(net, filename = NULL)
```

## Arguments

- net:

  `data.frame` or equivalent object representing the weighted network
  graph.

- filename:

  Name with optional full path of file in which to save the input `net`.
  The extension `.Rds` will be automatically appended, unless specified
  otherwise.

## Value

Nothing; function called for side-effect of saving network.

## Note

This may take some time if
[dodgr_cache_off](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_cache_off.md)
has been called. The contracted version of the graph is also saved, and
so must be calculated if it has not previously been automatically
cached.

## See also

Other cache:
[`clear_dodgr_cache()`](https://UrbanAnalyst.github.io/dodgr/reference/clear_dodgr_cache.md),
[`dodgr_cache_off()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_cache_off.md),
[`dodgr_cache_on()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_cache_on.md),
[`dodgr_load_streetnet()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_load_streetnet.md)

## Examples

``` r
net <- weight_streetnet (hampi)
f <- file.path (tempdir (), "streetnet.Rds")
dodgr_save_streetnet (net, f)
clear_dodgr_cache () # rm cached objects from tempdir
# at some later time, or in a new R session:
net <- dodgr_load_streetnet (f)
```
