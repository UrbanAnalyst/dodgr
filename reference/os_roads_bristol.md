# Sample street network from Bristol, U.K.

A sample street network for Bristol, U.K., from the Ordnance Survey.

## Format

A Simple Features `sf` `data.frame` representing motorways in Bristol,
UK.

## Note

Input data downloaded from <https://osdatahub.os.uk/downloads/open>, To
download the data from that page click on the tick box next to 'OS Open
Roads', scroll to the bottom, click 'Continue' and complete the form on
the subsequent page. This dataset is open access and can be used under
[these licensing
conditions](https://www.ordnancesurvey.co.uk/licensing), and must be
cited as follows: Contains OS data © Crown copyright and database right
(2017)

## See also

Other data:
[`hampi`](https://UrbanAnalyst.github.io/dodgr/reference/hampi.md),
[`weighting_profiles`](https://UrbanAnalyst.github.io/dodgr/reference/weighting_profiles.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library (sf)
library (dplyr)
# data must be unzipped here
# os_roads <- sf::read_sf("~/data/ST_RoadLink.shp")
# u <- paste0 (
#     "https://opendata.arcgis.com/datasets/",
#     "686603e943f948acaa13fb5d2b0f1275_4.kml"
# )
# lads <- sf::read_sf(u)
# mapview::mapview(lads)
# bristol_pol <- dplyr::filter(lads, grepl("Bristol", lad16nm))
# os_roads <- st_transform(os_roads, st_crs(lads)
# os_roads_bristol <- os_roads[bristol_pol, ] %>%
#   dplyr::filter(class == "Motorway" &
#                 roadNumber != "M32") %>%
#   st_zm(drop = TRUE)
# mapview::mapview(os_roads_bristol)
} # }
# Converting this 'sf data.frame' to a 'dodgr' network requires manual
# specification of weighting profile:
colnm <- "formOfWay" # name of column used to determine weights
wts <- data.frame (
    name = "custom",
    way = unique (os_roads_bristol [[colnm]]),
    value = c (0.1, 0.2, 0.8, 1)
)
net <- weight_streetnet (
    os_roads_bristol,
    wt_profile = wts,
    type_col = colnm, id_col = "identifier"
)
# 'id_col' tells the function which column to use to attribute IDs of ways
```
