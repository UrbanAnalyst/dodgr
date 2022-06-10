get_hash <- function (graph, verts = NULL, hash = TRUE) {

    if (hash) {
        hash <- attr (graph, "hash")
        if (is.null (hash)) {
            gr_cols <- dodgr_graph_cols (graph)
            hash <- digest::digest (graph [[gr_cols$edge_id]])
        }
    } else {
        hash <- attr (graph, "hashc")
        if (is.null (hash)) {
            gr_cols <- dodgr_graph_cols (graph)
            hash <- digest::digest (list (graph [[gr_cols$edge_id]], verts))
        }
    }
    return (hash)
}

get_edge_map <- function (graph) {

    hashc <- get_hash (graph, hash = FALSE)
    if (is.null (hashc)) {
        stop ("something went wrong extracting the edge map")
    } # nocov
    fname_c <- fs::path (
        fs::path_temp (),
        paste0 ("dodgr_edge_map_", hashc, ".Rds")
    )
    if (!fs::file_exists (fname_c)) {
        stop ("something went wrong extracting the edge map")
    } # nocov
    readRDS (fname_c)
}

# cache on initial construction with weight_streetnet. This pre-calculates and
# caches the contracted graph *with no additional intermediate vertices* (that
# is, the result of `dodgr_contract_graph (graph, verts = NULL)`). Later calls
# with explicit additional vertices will generate different hashes and so will
# be re-contracted and cached directly in `dodgr_contract_graph`.
#
# A copy of the original (full) graph is also copied to a file named with the
# hash of the edge map. This is needed for graph uncontraction, so that just the
# contracted graph and edge map can be submitted, the original graph re-loaded,
# and the uncontracted version returned.
cache_graph <- function (graph, edge_col) {

    td <- fs::path_temp ()

    f <- function (graph, edge_col, td) {

        # the following line does not generate a coverage symbol because it is
        # cached, so # nocov:
        verts <- dodgr::dodgr_vertices (graph) # nocov
        hash <- attr (graph, "hash")
        fname_v <- fs::path (td, paste0 ("dodgr_verts_", hash, ".Rds"))
        if (!fs::file_exists (fname_v)) {
            saveRDS (verts, fname_v)
        }

        # save original graph to enable subsequent re-loading from the
        # contracted version
        fname <- fs::path (td, paste0 ("dodgr_graph_", hash, ".Rds"))
        saveRDS (graph, fname)

        # The hash for the contracted graph is generated from the edge IDs of
        # the full graph plus default NULL vertices:
        hashc <- digest::digest (list (graph [[edge_col]], NULL))

        graphc <- dodgr::dodgr_contract_graph (graph)
        fname_c <- fs::path (td, paste0 ("dodgr_graphc_", hashc, ".Rds"))
        saveRDS (graphc, fname_c)

        hashe <- attr (graphc, "hashe")
        verts <- dodgr::dodgr_vertices (graphc)
        fname_v <- fs::path (td, paste0 ("dodgr_verts_", hashe, ".Rds"))
        saveRDS (verts, fname_v)

        fname_e <- paste0 ("dodgr_edge_map_", hashc, ".Rds")
        fname_e_fr <- fs::path (fs::path_temp (), fname_e)
        fname_e_to <- fs::path (td, fname_e)
        if (fs::file_exists (fname_e_fr)) { # should always be
            fs::file_copy (fname_e_fr, fname_e_to, overwrite = TRUE)
        }

        fname_j <- paste0 ("dodgr_junctions_", hashc, ".Rds")
        fname_j_fr <- fs::path (fs::path_temp (), fname_j)
        fname_j_to <- fs::path (td, fname_j)
        if (fs::file_exists (fname_j_fr)) { # should always be
            fs::file_copy (fname_j_fr, fname_j_to, overwrite = TRUE)
        }
    }

    sink (file = fs::path (fs::path_temp (), "Rout.txt"))
    res <- callr::r_bg (f, list (graph, edge_col, td))
    sink ()

    return (res) # R6 processx object
}

#' clear_dodgr_cache
#'
#' Remove cached versions of `dodgr` graphs. This function should generally
#' \emph{not} be needed, except if graph structure has been directly modified
#' other than through `dodgr` functions; for example by modifying edge weights
#' or distances. Graphs are cached based on the vector of edge IDs, so manual
#' changes to any other attributes will not necessarily be translated into
#' changes in `dodgr` output unless the cached versions are cleared using this
#' function. See
#' \url{https://github.com/ATFutures/dodgr/wiki/Caching-of-streetnets-and-contracted-graphs}
#' for details of caching process.
#'
#' @return Nothing; the function silently clears any cached objects
#' @family cache
#' @export
clear_dodgr_cache <- function () {

    lf <- list.files (fs::path_temp (), full.names = TRUE, pattern = "^dodgr_")
    if (length (lf) > 0) {
        tryCatch (
            chk <- file.remove (lf),
            error = function (e) NULL
        )
    }
}

#' dodgr_cache_off
#'
#' Turn off all dodgr caching in current session. This is useful is speed is
#' paramount, and if graph contraction is not needed. Caching can be switched
#' back on with \link{dodgr_cache_on}.
#' @return Nothing; the function invisibly returns `TRUE` if successful.
#' @family cache
#' @export
dodgr_cache_off <- function () {
    Sys.setenv ("DODGR_CACHE" = "FALSE")
}

#' dodgr_cache_on
#'
#' Turn on all dodgr caching in current session. This will only have an effect
#' after caching has been turned off with \link{dodgr_cache_off}.
#' @return Nothing; the function invisibly returns `TRUE` if successful.
#' @family cache
#' @export
dodgr_cache_on <- function () {
    Sys.setenv ("DODGR_CACHE" = "TRUE")
}



is_dodgr_cache_on <- function () {
    as.logical (Sys.getenv ("DODGR_CACHE", unset = "TRUE"))
}
