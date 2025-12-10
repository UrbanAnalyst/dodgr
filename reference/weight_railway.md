# Weight a network for routing along railways.

Weight (or re-weight) an `sf`-formatted OSM street network for routing
along railways.

## Usage

``` r
weight_railway(
  x,
  type_col = "railway",
  id_col = "osm_id",
  keep_cols = c("maxspeed"),
  excluded = c("abandoned", "disused", "proposed", "razed")
)
```

## Arguments

- x:

  A street network represented either as `sf` `LINESTRING` objects,
  typically extracted with
  [dodgr_streetnet](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_streetnet.md).

- type_col:

  Specify column of the `sf` `data.frame` object which designates
  different types of railways to be used for weighting (default works
  with `osmdata` objects).

- id_col:

  Specify column of the sf `data.frame` object which provides unique
  identifiers for each railway (default works with `osmdata` objects).

- keep_cols:

  Vectors of columns from `sf_lines` to be kept in the resultant `dodgr`
  network; vector can be either names or indices of desired columns.

- excluded:

  Types of railways to exclude from routing.

## Value

A `data.frame` of edges representing the rail network, along with a
column of graph component numbers.

## Note

Default railway weighting is by distance. Other weighting schemes, such
as by maximum speed, can be implemented simply by modifying the
`d_weighted` column returned by this function accordingly.

## See also

Other extraction:
[`dodgr_streetnet()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_streetnet.md),
[`dodgr_streetnet_geodesic()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_streetnet_geodesic.md),
[`dodgr_streetnet_sc()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_streetnet_sc.md),
[`weight_streetnet()`](https://UrbanAnalyst.github.io/dodgr/reference/weight_streetnet.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# sample railway extraction with the 'osmdata' package
library (osmdata)
dat <- opq ("shinjuku") %>%
    add_osm_feature (key = "railway") %>%
    osmdata_sf (quiet = FALSE)
graph <- weight_railway (dat$osm_lines)
} # }
```
