# Turn on all dodgr caching in current session.

This will only have an effect after caching has been turned off with
[dodgr_cache_off](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_cache_off.md).

## Usage

``` r
dodgr_cache_on()
```

## Value

Nothing; the function invisibly returns `TRUE` if successful.

## See also

Other cache:
[`clear_dodgr_cache()`](https://UrbanAnalyst.github.io/dodgr/reference/clear_dodgr_cache.md),
[`dodgr_cache_off()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_cache_off.md),
[`dodgr_load_streetnet()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_load_streetnet.md),
[`dodgr_save_streetnet()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_save_streetnet.md)

## Examples

``` r
dodgr_cache_on ()
# Then call dodgr functions as usual:
graph <- weight_streetnet (hampi, wt_profile = "foot")
```
