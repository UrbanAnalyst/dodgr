# Convert a `dodgr` graph into an equivalent `sf::sfc` object.

Convert a `dodgr` graph into a `list` composed of two objects: `dat`, a
`data.frame`; and `geometry`, an `sfc` object from the (sf) package.
Works by aggregating edges into `LINESTRING` objects representing
longest sequences between all junction nodes. The resultant objects will
generally contain more `LINESTRING` objects than the original sf object,
because the former will be bisected at every junction point.

## Usage

``` r
dodgr_to_sfc(graph)
```

## Arguments

- graph:

  A `dodgr` graph

## Value

A list containing (1) A `data.frame` of data associated with the `sf`
geometries; and (ii) A Simple Features Collection (`sfc`) list of
`LINESTRING` objects.

## Note

The output of this function corresponds to the edges obtained from
`dodgr_contract_graph`. This function does not require the sf package to
be installed; the corresponding function that creates a full sf object -
[dodgr_to_sf](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_to_sf.md)
does requires sf to be installed.

## See also

Other conversion:
[`dodgr_deduplicate_graph()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_deduplicate_graph.md),
[`dodgr_to_igraph()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_to_igraph.md),
[`dodgr_to_sf()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_to_sf.md),
[`dodgr_to_tidygraph()`](https://UrbanAnalyst.github.io/dodgr/reference/dodgr_to_tidygraph.md),
[`igraph_to_dodgr()`](https://UrbanAnalyst.github.io/dodgr/reference/igraph_to_dodgr.md)

## Examples

``` r
hw <- weight_streetnet (hampi)
nrow (hw)
#> [1] 6813
xy <- dodgr_to_sfc (hw)
dim (hw) # 5.845 edges
#> [1] 6813   16
length (xy$geometry) # more linestrings aggregated from those edges
#> [1] 744
nrow (hampi) # than the 191 linestrings in original sf object
#> [1] 236
dim (xy$dat) # same number of rows as there are geometries
#> [1] 744  16
# The dodgr_to_sf function then just implements this final conversion:
# sf::st_sf (xy$dat, geometry = xy$geometry, crs = 4326)
```
