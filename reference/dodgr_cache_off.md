# Turn off all dodgr caching in current session.

This function is useful is speed is paramount, and if graph contraction
is not needed. Caching can be switched back on with
[dodgr_cache_on](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_cache_on.md).

## Usage

``` r
dodgr_cache_off()
```

## Value

Nothing; the function invisibly returns `TRUE` if successful.

## See also

Other cache:
[`clear_dodgr_cache()`](https://UrbanAnalyst.github.io/dodgr/reference/clear_dodgr_cache.md),
[`dodgr_cache_on()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_cache_on.md),
[`dodgr_load_streetnet()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_load_streetnet.md),
[`dodgr_save_streetnet()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_save_streetnet.md)

## Examples

``` r
dodgr_cache_off ()
# Then call dodgr functions as usual:
graph <- weight_streetnet (hampi, wt_profile = "foot")
```
