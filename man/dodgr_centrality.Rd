% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/centrality.R
\name{dodgr_centrality}
\alias{dodgr_centrality}
\title{Calculate betweenness centrality for a 'dodgr' network.}
\usage{
dodgr_centrality(
  graph,
  contract = TRUE,
  edges = TRUE,
  column = "d_weighted",
  vert_wts = NULL,
  dist_threshold = NULL,
  heap = "BHeap",
  check_graph = TRUE
)
}
\arguments{
\item{graph}{'data.frame' or equivalent object representing the network
graph (see Details)}

\item{contract}{If 'TRUE', centrality is calculated on contracted graph
before mapping back on to the original full graph. Note that for street
networks, in particular those obtained from the \pkg{osmdata} package, vertex
placement is effectively arbitrary except at junctions; centrality for such
graphs should only be calculated between the latter points, and thus
'contract' should always be 'TRUE'.}

\item{edges}{If 'TRUE', centrality is calculated for graph edges, returning
the input 'graph' with an additional 'centrality' column; otherwise
centrality is calculated for vertices, returning the equivalent of
'dodgr_vertices(graph)', with an additional vertex-based 'centrality' column.}

\item{column}{Column of graph defining the edge properties used to calculate
centrality (see Note).}

\item{vert_wts}{Optional vector of length equal to number of vertices
(\code{nrow(dodgr_vertices(graph))}), to enable centrality to be calculated in
weighted form, such that centrality measured from each vertex will be
weighted by the specified amount.}

\item{dist_threshold}{If not 'NULL', only calculate centrality for each point
out to specified threshold. Setting values for this will result in
approximate estimates for centrality, yet with considerable gains in
computational efficiency. For sufficiently large values, approximations will
be accurate to within some constant multiplier. Appropriate values can be
established via the \link{estimate_centrality_threshold} function.}

\item{heap}{Type of heap to use in priority queue. Options include
Fibonacci Heap (default; 'FHeap'), Binary Heap ('BHeap'),
Trinomial Heap ('TriHeap'), Extended Trinomial Heap
('TriHeapExt', and 2-3 Heap ('Heap23').}

\item{check_graph}{If \code{TRUE}, graph is first checked for duplicate edges,
which can cause incorrect centrality calculations. If duplicate edges are
detected in an interactive session, a prompt will ask whether you want to
proceed or rectify edges first. This value may be set to \code{FALSE} to skip this
check and the interactive prompt.}
}
\value{
Modified version of graph with additional 'centrality' column added.
}
\description{
Centrality can be calculated in either vertex- or edge-based form.
}
\note{
The \code{column} parameter is by default \code{d_weighted}, meaning centrality
is calculated by routing according to weighted distances. Other possible
values for this parameter are
\itemize{
\item \code{d} for unweighted distances
\item \code{time} for unweighted time-based routing
\item \code{time_weighted} for weighted time-based routing
}

Centrality is calculated by default using parallel computation with the
maximal number of available cores or threads. This number can be reduced by
specifying a value via
\verb{RcppParallel::setThreadOptions (numThreads = <desired_number>)}.
}
\examples{
\dontrun{
graph_full <- weight_streetnet (hampi)
graph <- dodgr_contract_graph (graph_full)
graph <- dodgr_centrality (graph)
# 'graph' is then the contracted graph with an additional 'centrality' column
# Same calculation via 'igraph':
igr <- dodgr_to_igraph (graph)
library (igraph)
cent <- edge_betweenness (igr)
identical (cent, graph$centrality) # TRUE
# Values of centrality between all junctions in the contracted graph can then
# be mapped back onto the original full network by "uncontracting":
graph_full <- dodgr_uncontract_graph (graph)
# For visualisation, it is generally necessary to merge the directed edges to
# form an equivalent undirected graph. Conversion to 'sf' format via
# 'dodgr_to_sf()' is also useful for many visualisation routines.
graph_sf <- merge_directed_graph (graph_full) \%>\%
    dodgr_to_sf ()
}

\dontrun{
library (mapview)
centrality <- graph_sf$centrality / max (graph_sf$centrality)
ncols <- 30
cols <- c ("lawngreen", "red")
cols <- colorRampPalette (cols) (ncols) [ceiling (ncols * centrality)]
mapview (graph_sf, color = cols, lwd = 10 * centrality)
}

# An example of flow aggregation across a generic (non-OSM) highway,
# represented as the 'routes_fast' object of the \pkg{stplanr} package,
# which is a SpatialLinesDataFrame containing commuter densities along
# components of a street network.
\dontrun{
library (stplanr)
# merge all of the 'routes_fast' lines into a single network
r <- overline (routes_fast, attrib = "length", buff_dist = 1)
r <- sf::st_as_sf (r)
# Convert to a 'dodgr' network, for which we need to specify both a 'type'
# and 'id' column.
r$type <- 1
r$id <- seq (nrow (r))
graph_full <- weight_streetnet (
    r,
    type_col = "type",
    id_col = "id",
    wt_profile = 1
)
# convert to contracted form, retaining junction vertices only, and append
# 'centrality' column
graph <- dodgr_contract_graph (graph_full) \%>\%
    dodgr_centrality ()
#' expand back to full graph; merge directed flows; and convert result to
# 'sf'-format for plotting
graph_sf <- dodgr_uncontract_graph (graph) \%>\%
    merge_directed_graph () \%>\%
    dodgr_to_sf ()
plot (graph_sf ["centrality"])
}

}
\seealso{
Other centrality: 
\code{\link{estimate_centrality_threshold}()},
\code{\link{estimate_centrality_time}()}
}
\concept{centrality}
