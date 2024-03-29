% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph-functions.R
\name{dodgr_components}
\alias{dodgr_components}
\title{Identify connected components of graph.}
\usage{
dodgr_components(graph)
}
\arguments{
\item{graph}{A \code{data.frame} of edges}
}
\value{
Equivalent graph with additional \code{component} column,
sequentially numbered from 1 = largest component.
}
\description{
Identify connected components of graph and add corresponding \code{component}
column to \code{data.frame}.
}
\examples{
graph <- weight_streetnet (hampi)
graph <- dodgr_components (graph)
}
\seealso{
Other modification: 
\code{\link{dodgr_contract_graph}()},
\code{\link{dodgr_uncontract_graph}()}
}
\concept{modification}
