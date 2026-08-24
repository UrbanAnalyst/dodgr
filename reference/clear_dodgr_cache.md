# Remove cached versions of `dodgr` graphs.

This function should generally *not* be needed, except if graph
structure has been directly modified other than through `dodgr`
functions. Graphs are cached based on a hash of the `edge_id`, `d`,
`d_weighted`, `time`, and `time_weighted` columns, so manual changes to
any of those will be detected automatically. Manual changes to any other
column or attribute (for example, `from`/`to` vertex identifiers) will
not be reflected in the hash, and so will not necessarily be translated
into changes in `dodgr` output unless the cached versions are cleared
using this function. See
<https://github.com/UrbanAnalyst/dodgr/wiki/Caching-of-streetnets-and-contracted-graphs>
\# nolint for details of caching process.

## Usage

``` r
clear_dodgr_cache()
```

## Value

Nothing; the function silently clears any cached objects

## See also

Other cache:
[`dodgr_cache_off()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_cache_off.md),
[`dodgr_cache_on()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_cache_on.md),
[`dodgr_load_streetnet()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_load_streetnet.md),
[`dodgr_save_streetnet()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_save_streetnet.md)

## Examples

``` r
clear_dodgr_cache ()
#> [1] TRUE TRUE TRUE TRUE TRUE TRUE
# Then call dodgr functions as usual:
graph <- weight_streetnet (hampi, wt_profile = "foot")
```
