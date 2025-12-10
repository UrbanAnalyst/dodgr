# Weighting profiles used to route different modes of transport.

Collection of weighting profiles used to adjust the routing process to
different means of transport. Modified from data taken from the Routino
project, with additional tables for average speeds, dependence of speed
on type of surface, and waiting times in seconds at traffic lights. The
latter table (called "penalties") includes waiting times at traffic
lights (in seconds), additional time penalties for turning across
oncoming traffic ("turn"), and a binary flag indicating whether turn
restrictions should be obeyed or not.

## Format

List of `data.frame` objects with profile names, means of transport and
weights.

## References

<https://www.routino.org/xml/routino-profiles.xml>

## See also

Other data:
[`hampi`](https://UrbanAnalyst.github.io/dodgr/reference/hampi.md),
[`os_roads_bristol`](https://UrbanAnalyst.github.io/dodgr/reference/os_roads_bristol.md)
