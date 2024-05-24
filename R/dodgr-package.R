#' Distances On Directed GRaphs ("dodgr")
#'
#' Distances on dual-weighted directed graphs using priority-queue shortest
#' paths. Weighted directed graphs have weights from A to B which may differ
#' from those from B to A. Dual-weighted directed graphs have two sets of such
#' weights. A canonical example is a street network to be used for routing in
#' which routes are calculated by weighting distances according to the type of
#' way and mode of transport, yet lengths of routes must be calculated from
#' direct distances.
#'
#' @section The Main Function:
#' \itemize{
#' \item [dodgr_dists()]: Calculate pair-wise distances between
#' specified pairs of points in a graph.
#' }
#'
#' @section Functions to Obtain Graphs:
#' \itemize{
#' \item [dodgr_streetnet()]: Extract a street network in Simple
#' Features (`sf`) form.
#' \item [weight_streetnet()]: Convert an `sf`-formatted street
#' network to a `dodgr` graph through applying specified weights to all
#' edges.
#' }
#'
#' @section Functions to Modify Graphs:
#' \itemize{
#' \item [dodgr_components()]: Number all graph edges according to
#' their presence in distinct connected components.
#' \item [dodgr_contract_graph()]: Contract a graph by removing
#' redundant edges.
#' }
#'
#' @section Miscellaneous Functions:
#' \itemize{
#' \item [dodgr_sample()]: Randomly sample a graph, returning a single
#' connected component of a defined number of vertices.
#' \item [dodgr_vertices()]: Extract all vertices of a graph.
#' \item [compare_heaps()]: Compare the performance of different
#' priority queue heap structures for a given type of graph.
#' }
#'
#' @name dodgr
#' @docType package
#' @family package
#' @importFrom graphics plot
#' @importFrom grDevices colorRampPalette
#' @importFrom magrittr %>%
#' @importFrom methods is
#' @importFrom osmdata add_osm_feature getbb opq osmdata_sf osm_poly2line
#' @importFrom osmdata trim_osmdata
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
#' @useDynLib dodgr, .registration = TRUE
"_PACKAGE"

#' Weighting profiles used to route different modes of transport.
#'
#' Collection of weighting profiles used to adjust the routing process to
#' different means of transport. Modified from data taken from the Routino
#' project, with additional tables for average speeds, dependence of speed on
#' type of surface, and waiting times in seconds at traffic lights. The latter
#' table (called "penalties") includes waiting times at traffic lights (in
#' seconds), additional time penalties for turning across oncoming traffic
#' ("turn"), and a binary flag indicating whether turn restrictions should be
#' obeyed or not.
#'
#' @name weighting_profiles
#' @docType data
#' @family data
#' @keywords datasets
#' @format List of `data.frame` objects with profile names, means of transport
#' and weights.
#' @references <https://www.routino.org/xml/routino-profiles.xml>
NULL

#' Sample street network from Hampi, India.
#'
#' A sample street network from the township of Hampi, Karnataka, India.
#'
#' @name hampi
#' @docType data
#' @family data
#' @keywords datasets
#' @format A Simple Features `sf` `data.frame` containing the street
#' network of Hampi.
#'
#' @note Can be re-created with the following command, which also removes
#' extraneous columns to reduce size:
#' @examples \dontrun{
#' hampi <- dodgr_streetnet ("hampi india")
#' cols <- c ("osm_id", "highway", "oneway", "geometry")
#' hampi <- hampi [, which (names (hampi) %in% cols)]
#' }
#' # this 'sf data.frame' can be converted to a 'dodgr' network with
#' net <- weight_streetnet (hampi, wt_profile = "foot")
NULL

#' Sample street network from Bristol, U.K.
#'
#' A sample street network for Bristol, U.K., from the Ordnance Survey.
#'
#' @name os_roads_bristol
#' @docType data
#' @keywords datasets
#' @format A Simple Features `sf` `data.frame` representing
#' motorways in Bristol, UK.
#'
#' @note Input data downloaded from
#' \url{https://osdatahub.os.uk/downloads/open},
#' To download the data from that page click on the tick box next to
#' 'OS Open Roads', scroll to the bottom, click 'Continue' and complete
#' the form on the subsequent page.
#' This dataset is open access and can be used under
#' \href{https://www.ordnancesurvey.co.uk/licensing}{these licensing
#' conditions},
#' and must be cited as follows:
#' Contains OS data © Crown copyright and database right (2017)
#'
#' @family data
#' @examples \dontrun{
#' library (sf)
#' library (dplyr)
#' # data must be unzipped here
#' # os_roads <- sf::read_sf("~/data/ST_RoadLink.shp")
#' # u <- paste0 (
#' #     "https://opendata.arcgis.com/datasets/",
#' #     "686603e943f948acaa13fb5d2b0f1275_4.kml"
#' # )
#' # lads <- sf::read_sf(u)
#' # mapview::mapview(lads)
#' # bristol_pol <- dplyr::filter(lads, grepl("Bristol", lad16nm))
#' # os_roads <- st_transform(os_roads, st_crs(lads)
#' # os_roads_bristol <- os_roads[bristol_pol, ] %>%
#' #   dplyr::filter(class == "Motorway" &
#' #                 roadNumber != "M32") %>%
#' #   st_zm(drop = TRUE)
#' # mapview::mapview(os_roads_bristol)
#' }
#' # Converting this 'sf data.frame' to a 'dodgr' network requires manual
#' # specification of weighting profile:
#' colnm <- "formOfWay" # name of column used to determine weights
#' wts <- data.frame (
#'     name = "custom",
#'     way = unique (os_roads_bristol [[colnm]]),
#'     value = c (0.1, 0.2, 0.8, 1)
#' )
#' net <- weight_streetnet (
#'     os_roads_bristol,
#'     wt_profile = wts,
#'     type_col = colnm, id_col = "identifier"
#' )
#' # 'id_col' tells the function which column to use to attribute IDs of ways
NULL

#' Pipe operator
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL
