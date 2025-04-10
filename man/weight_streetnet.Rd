% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weight-streetnet.R
\name{weight_streetnet}
\alias{weight_streetnet}
\alias{weight_streetnet.default}
\alias{weight_streetnet.sf}
\alias{weight_streetnet.sc}
\alias{weight_streetnet.SC}
\title{Weight a street network according to a specified weighting profile.}
\usage{
weight_streetnet(
  x,
  wt_profile = "bicycle",
  wt_profile_file = NULL,
  turn_penalty = FALSE,
  type_col = "highway",
  id_col = "osm_id",
  keep_cols = NULL,
  left_side = FALSE
)

\method{weight_streetnet}{default}(
  x,
  wt_profile = "bicycle",
  wt_profile_file = NULL,
  turn_penalty = FALSE,
  type_col = "highway",
  id_col = "osm_id",
  keep_cols = NULL,
  left_side = FALSE
)

\method{weight_streetnet}{sf}(
  x,
  wt_profile = "bicycle",
  wt_profile_file = NULL,
  turn_penalty = FALSE,
  type_col = "highway",
  id_col = "osm_id",
  keep_cols = NULL,
  left_side = FALSE
)

\method{weight_streetnet}{sc}(
  x,
  wt_profile = "bicycle",
  wt_profile_file = NULL,
  turn_penalty = FALSE,
  type_col = "highway",
  id_col = "osm_id",
  keep_cols = NULL,
  left_side = FALSE
)

\method{weight_streetnet}{SC}(
  x,
  wt_profile = "bicycle",
  wt_profile_file = NULL,
  turn_penalty = FALSE,
  type_col = "highway",
  id_col = "osm_id",
  keep_cols = NULL,
  left_side = FALSE
)
}
\arguments{
\item{x}{A street network represented either as \code{sf} \code{LINESTRING}
objects, typically extracted with \link{dodgr_streetnet}, or as an \code{SC}
(\code{silicate}) object typically extracted with the \link{dodgr_streetnet_sc}.}

\item{wt_profile}{Name of weighting profile, or \code{data.frame} specifying
custom values (see Details)}

\item{wt_profile_file}{Name of locally-stored, \code{.json}-formatted version of
\code{dodgr::weighting_profiles}, created with \link{write_dodgr_wt_profile}, and
modified as desired.}

\item{turn_penalty}{Including time penalty on edges for turning across
oncoming traffic at intersections (see Note).}

\item{type_col}{Specify column of the \code{sf} \code{data.frame} object
which designates different types of highways to be used for weighting
(default works with \code{osmdata} objects).}

\item{id_col}{For \code{sf}-formatted data only: Specify column of the \pkg{sf}
\code{data.frame} object which provides unique identifiers for each highway
(default works with \code{osmdata} objects).}

\item{keep_cols}{Vectors of columns from \code{x} to be kept in the resultant
\code{dodgr} network; vector can be either names, regex-patterns,  or indices of
desired columns (see notes).}

\item{left_side}{Does traffic travel on the left side of the road (\code{TRUE}) or
the right side (\code{FALSE})? - only has effect on turn angle calculations for
edge times.}
}
\value{
A \code{data.frame} of edges representing the street network, with
distances in metres and times in seconds, along with a column of graph
component numbers. Times for \pkg{sf}-formatted street networks are only
approximate, and do not take into account traffic lights, turn angles, or
elevation changes. Times for \pkg{sc}-formatted street networks take into
account all of these factors, with elevation changes automatically taken into
account for networks generated with the \pkg{osmdata} function
\code{osm_elevation()}.
}
\description{
Weight (or re-weight) an \pkg{sf} or \pkg{silicate} ("\emph{SC}") formatted OSM
street network according to a specified weighting  profile. Standard
weighting profiles may be specified by name, as one of:
\itemize{
\item foot
\item horse
\item wheelchair
\item bicycle
\item moped
\item motorcycle
\item motorcar
\item goods
\item hgv
\item psv
}
Custom weighting profiles are also possible, as explained in the Note below.
}
\details{
Distances along each edge are calculated using the \pkg{geodist} package,
defaulting to the Mapbox "cheap" metric if maximal network distances are <
100km, otherwise using Haversine distances. The
\link{dodgr_streetnet_geodesic} function can be called to force edge
distances to be calculated using more accurate yet slower geodesic
distances.

The structure of networks generated by this function is dependent
on many aspects of the input network, and in particular on specific
key-value pairs defined in the underlying OpenStreetMap (OSM) data.

Many key-value pairs influence the resultant network through being used in
specified weighting profiles. Keys used in weighting profiles are always
kept in the weighted networks, and are specified in
\link{weighting_profiles} by the "way" column in the "weighting_profiles"
item. These include:

\itemize{
\item "bridleway"
\item "cycleway"
\item "ferry"
\item "footway"
\item "living_street"
\item "motorway"
\item "motorway_link
\item "path"
\item "pedestrian"
\item "primary"
\item "primary_link"
\item "residential"
\item "secondary"
\item "secondary_link
\item "service"
\item "steps"
\item "tertiary"
\item "tertiary_link"
\item "track"
\item "trunk"
\item "trunk_link
\item "unclassified"
}

Some of these are only optionally kept, dependent on the weighting profile
chosen. For example, "cycleway" keys are only kept for bicycle weighting.
Most of the specified keys also include all possible variations on those
keys. For the example of "cycleway" again, key-value pairs are also kept for
"cycleway:left" and "cycleway:right".

The following additional keys are also automatically retained in weighted
networks:

\itemize{
\item "highway"
\item "junction"
\item "lanes"
\item "maxspeed"
\item "oneway", including with all alternative forms such as "oneway.bicycle"
\item "surface"
}

Realistic routing including factors such as access restrictions, turn
penalties, and effects of incline, can only be implemented when the objects
passed to \code{weight_streetnet} are of \pkg{sc} ("silicate") format, generated
with \link{dodgr_streetnet_sc} (and possibly enhanced through applying the
\pkg{osmdata} function \code{osm_elevation()}). Restrictions applied to ways (in
OSM terminology) may be controlled by ensuring specific columns are retained
in the \code{dodgr} network with the \code{keep_cols} argument. For example,
restrictions on access are generally specified by specifying a value for the
key of "access". Include "access" in \code{keep_cols} will ensure these values
are retained in the \code{dodgr} version, from which ways with specified values
can easily be removed or modified, as demonstrated in the examples.

Restrictions and time-penalties on turns can be implemented by setting
\code{turn_penalty = TRUE}, which will then honour turn restrictions specified in
OSM (unless the "penalties" table of \link{weighting_profiles} has
\code{restrictions = FALSE} for a specified \code{wt_profile}). Resultant graphs are
fundamentally different from the default for distance-based routing. These
graphs may be used directly in most 'dodgr' functions, but generally only if
they have been created by calling this function in the same session, or if
they have been saved and loaded with the \link{dodgr_save_streetnet} and
\link{dodgr_load_streetnet} functions. (This is because the weighted
streetnets also have accompanying data stored in a local temporary cache
directory; attempting to pass a weighted street network without these
accompanying cache files will generally error.)

Some key-value pairs may also directly influence not just the value of the
graph produced by this function, but also its size. Among these are "oneway"
flags. Without these flags, each edge will be represented in \emph{directed} form,
and so as two rows of the graph: one for A -> B, and one for B -> A. If a
way is tagged in OSM as "oneway" = "yes", and if oneway flags are respected
for a chosen weighting profile (which, for example, they are generally not
for pedestrian or "foot" weighting), then only one edge will be returned
representing travel in the direction permitted within the OSM data. Thus
weighting a network which includes "oneway" flags, and using a weighting
profile which respects these, will generate a graph with fewer rows than a
graph produced by ignoring those "oneway" flags.
}
\note{
Names for the \code{wt_profile} parameter are taken from
\link{weighting_profiles}, which is a list including a \code{data.frame} also
called \code{weighting_profiles} of weights for different modes of transport.
Values for \code{wt_profile} are taken from current modes included there, which
are "bicycle", "foot", "goods", "hgv", "horse", "moped", "motorcar",
"motorcycle", "psv", and "wheelchair". Railway routing can be implemented
with the separate function \link{weight_railway}. Alternatively, the entire
\code{weighting_profile} structures can be written to a local \code{.json}-formatted
file with \link{write_dodgr_wt_profile}, the values edited as desired, and
the name of this file passed as the \code{wt_profile_file} parameter.

The resultant graph includes only those edges for which the given
weighting profile specifies finite edge weights. Any edges of types not
present in a given weighting profile are automatically removed from the
weighted streetnet.

If the resultant graph is to be contracted via
\link{dodgr_contract_graph}, \strong{and} if the columns of the graph have been,
or will be, modified, then automatic caching must be switched off with
\link{dodgr_cache_off}. If not, the \link{dodgr_contract_graph} function will
return the automatically cached version, which is the contracted version of
the full graph prior to any modification of columns.
}
\examples{
# hampi is included with package as an 'osmdata' sf-formatted street network
net <- weight_streetnet (hampi)
class (net) # data.frame
dim (net) # 6096  11; 6096 streets
# os_roads_bristol is also included as an sf data.frame, but in a different
# format requiring identification of columns and specification of custom
# weighting scheme.
colnm <- "formOfWay"
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
dim (net) # 406 11; 406 streets

# An example for a generic (non-OSM) highway, represented as the
# `routes_fast` object of the \pkg{stplanr} package, which is a
# SpatialLinesDataFrame.
\dontrun{
library (stplanr)
# merge all of the 'routes_fast' lines into a single network
r <- overline (routes_fast, attrib = "length", buff_dist = 1)
r <- sf::st_as_sf (r, crs = 4326)
# We need to specify both a `type` and `id` column for the
# \link{weight_streetnet} function.
r$type <- 1
r$id <- seq (nrow (r))
graph <- weight_streetnet (
    r,
    type_col = "type",
    id_col = "id",
    wt_profile = 1
)
}
}
\seealso{
\link{write_dodgr_wt_profile}, \link{dodgr_times}

Other extraction: 
\code{\link{dodgr_streetnet}()},
\code{\link{dodgr_streetnet_geodesic}()},
\code{\link{dodgr_streetnet_sc}()},
\code{\link{weight_railway}()}

Other extraction: 
\code{\link{dodgr_streetnet}()},
\code{\link{dodgr_streetnet_geodesic}()},
\code{\link{dodgr_streetnet_sc}()},
\code{\link{weight_railway}()}

Other extraction: 
\code{\link{dodgr_streetnet}()},
\code{\link{dodgr_streetnet_geodesic}()},
\code{\link{dodgr_streetnet_sc}()},
\code{\link{weight_railway}()}

Other extraction: 
\code{\link{dodgr_streetnet}()},
\code{\link{dodgr_streetnet_geodesic}()},
\code{\link{dodgr_streetnet_sc}()},
\code{\link{weight_railway}()}

Other extraction: 
\code{\link{dodgr_streetnet}()},
\code{\link{dodgr_streetnet_geodesic}()},
\code{\link{dodgr_streetnet_sc}()},
\code{\link{weight_railway}()}
}
\concept{extraction}
