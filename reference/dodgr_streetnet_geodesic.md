# Force [weight_streetnet](https://UrbanAnalyst.github.io/dodgr/reference/weight_streetnet.md) to use geodesic distances.

Distances by default are Mapbox "cheap" distances if maximal network
distances are \< 100km, otherwise Haversine distances. Calling this
function forces all calls to
[weight_streetnet](https://UrbanAnalyst.github.io/dodgr/reference/weight_streetnet.md)
from that point on to use geodesic distances. These are more
computationally expensive to calculate, and weighting networks will
likely take more time.

## Usage

``` r
dodgr_streetnet_geodesic(unset = FALSE)
```

## Arguments

- unset:

  Calling this function with `unset = TRUE` reverts distance
  calculations to those described above, rather than geodesic.

## Value

Nothing; the function is called for its side-effect only of setting
distance calculations to geodesic.

## See also

Other extraction:
[`dodgr_streetnet()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_streetnet.md),
[`dodgr_streetnet_sc()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_streetnet_sc.md),
[`weight_railway()`](https://UrbanAnalyst.github.io/dodgr/reference/weight_railway.md),
[`weight_streetnet()`](https://UrbanAnalyst.github.io/dodgr/reference/weight_streetnet.md)

## Examples

``` r
net0 <- weight_streetnet (hampi) # Default "cheap" method
dodgr_streetnet_geodesic ()
net1 <- weight_streetnet (hampi)
cor (net0$d, net1$d) # Strongly correlated, but not perfect
#> [1] 0.9999972
max (abs (net0$d - net1$d)) # in metres
#> [1] 1.622769
```
