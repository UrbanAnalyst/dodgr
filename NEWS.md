# v0.1.1

Major changes:
- Crucial fix of previous typo that made all `dodgr_dist` calculations wrong
  (Earth's radius is 6371, not 3671!) - thanks to @chrijo
- `weight_streetnet` function now accepts `data.frame` objects defining
  `wt_profile`, enabling modification and direct re-submission of
  `dodgr::weighting_profiles`
- `weighting_profiles$value` modified to 0-1 scores rather than previous
  percentage values.
- `weight_streetnet` flags any highway types not present in nominated or
  submitted weighting profile.
- `dodgr_paths` now has additional `pairwise` parameter to enable paths only
  between matched pairs of `from` and `to` points (so returning `n` paths rather
  than `n^2`), thanks to @mem48.
- `dodgr_to_sf` deprecated to `dodgr_to_sfc` (#43)

Minor changes:
- Added Malcolm Morgan (@mem48; bug-finder extraordinare) as contributor 
- Bug fix with `dodgr_paths` and simple `data.frame`s, thanks to James Smith.
- Bug fix of former improper handling of one-way edges, thanks to @chrijo.
- `match_pts_to_graph` has additional `connected` parameter to allow points to
  be matched only to largest connected component.

# v0.1.0

Major changes:
- New function `dodgr_flowmap` plots maps of flows. Currently only writes .png
  files, because large networks can not be effectively plotted on graphic
  devices.
- `dodgr_flows` has option to routes flows from a set of source origins to all
  points in a network, attenuated by distance from those origins.
- `dodgr_to_sf` converts a spatially-explicit `dodgr` graph into Simple Features
  (`sf`) format.

Minor changes:
- `match_pts_to_graph` now accepts Simple Features (sf) collections of
  `sfc_POINT` objects to be matched.

# v0.0.3

Tidy C++ code that flagged errors on CRAN solaris machine. Nothing else.

# v0.0.2

Major changes:
- New function, `dodgr_paths`, for returning explicit shortest path routes.
- New function, `dodgr_flows`, for aggregting flows across a network from
  multiple origin and destination points.
- New function, `merge_directed_flows`, to reduce aggregated directional flows
  to non-directional equivalent values useful for visualisation.

Minor changes:
- `weight_streetnet` now accepts arbitrary `sf`-formatted networks via
  specification of custom weighting profiles, along with highway type and ID
  columns in data.frame.
- Distance matrices from `dodgr_dists` inherit the names of routing points
  (`from` and `to` parameters).

# v0.0.1

Initial CRAN release.
