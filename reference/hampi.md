# Sample street network from Hampi, India.

A sample street network from the township of Hampi, Karnataka, India.

## Format

A Simple Features `sf` `data.frame` containing the street network of
Hampi.

## Note

Can be re-created with the following command, which also removes
extraneous columns to reduce size:

## See also

Other data:
[`os_roads_bristol`](https://UrbanAnalyst.github.io/dodgr/reference/os_roads_bristol.md),
[`weighting_profiles`](https://UrbanAnalyst.github.io/dodgr/reference/weighting_profiles.md)

## Examples

``` r
if (FALSE) { # \dontrun{
hampi <- dodgr_streetnet ("hampi india")
cols <- c ("osm_id", "highway", "oneway", "geometry")
hampi <- hampi [, which (names (hampi) %in% cols)]
} # }
# this 'sf data.frame' can be converted to a 'dodgr' network with
net <- weight_streetnet (hampi, wt_profile = "foot")
```
