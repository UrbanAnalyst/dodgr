% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph-functions.R
\name{dodgr_sample}
\alias{dodgr_sample}
\title{Sample a random but connected sub-component of a graph}
\usage{
dodgr_sample(graph, nverts = 1000)
}
\arguments{
\item{graph}{A flat table of graph edges. Must contain columns labelled
\code{from} and \code{to}, or \code{start} and \code{stop}. May also contain
similarly labelled columns of spatial coordinates (for example
\code{from_x}) or \code{stop_lon}).}

\item{nverts}{Number of vertices to sample}
}
\value{
A connected sub-component of \code{graph}
}
\description{
Sample a random but connected sub-component of a graph
}
\note{
Graphs may occasionally have \code{nverts + 1} vertices, rather than
the requested \code{nverts}.
}
\examples{
graph <- weight_streetnet (hampi)
nrow (graph) # 5,742
graph <- dodgr_sample (graph, nverts = 200)
nrow (graph) # generally around 400 edges
nrow (dodgr_vertices (graph)) # 200
}
\seealso{
Other misc: 
\code{\link{compare_heaps}()},
\code{\link{dodgr_flowmap}()},
\code{\link{dodgr_full_cycles}()},
\code{\link{dodgr_fundamental_cycles}()},
\code{\link{dodgr_insert_vertex}()},
\code{\link{dodgr_sflines_to_poly}()},
\code{\link{dodgr_vertices}()},
\code{\link{merge_directed_graph}()},
\code{\link{summary.dodgr_dists_categorical}()},
\code{\link{write_dodgr_wt_profile}()}
}
\concept{misc}
