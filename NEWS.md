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
