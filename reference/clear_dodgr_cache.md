# Remove cached versions of `dodgr` graphs.

This function should generally *not* be needed, except if graph
structure has been directly modified other than through `dodgr`
functions; for example by modifying edge weights or distances. Graphs
are cached based on the vector of edge IDs, so manual changes to any
other attributes will not necessarily be translated into changes in
`dodgr` output unless the cached versions are cleared using this
function. See
<https://github.com/UrbanAnalyst/dodgr/wiki/Caching-of-streetnets-and-contracted-graphs>
for details of caching process.

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
# Then call dodgr functions as usual:
graph <- weight_streetnet (hampi, wt_profile = "foot")
```
