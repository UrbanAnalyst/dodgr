null_to_na <- function (x)
{
    if (length (x) == 0)
        x <- NA
    return (x)
}

get_hash <- function (graph)
{
    hash_name <- ifelse (methods::is (graph, "dodgr_contracted"),
                         "hashc", "hash")
    hash <- attr (graph, hash_name)
    if (is.null (hash))
        hash <- digest::digest (graph)
    return (hash)
}

get_hashc <- function (graph, verts = NULL)
{
    hash <- attr (graph, "hash")
    if (is.null (hash))
        hash <- digest::digest (graph)

    if (is.null (verts))
        hashc <- hash
    else
        hashc <- digest::digest (list (graph, verts))
    return (hashc)
}

get_edge_map <- function (graph)
{
    hashc <- attr (graph, "hashc")
    if (is.null (hashc))
        stop ("something went wrong extracting the edge map")
    fname_e <- file.path (tempdir (), paste0 ("edge_map_", hashc, ".Rds"))
    if (!file.exists (fname_e))
        stop ("something went wrong extracting the edge map")
    readRDS (fname_e)
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
cache_graph <- function (graph, hash)
{
    td <- tempdir ()
    f <- function (graph, hash, td)
    {
        fname <- file.path (td, paste0 ("graph_", hash, ".Rds"))
        saveRDS (graph, fname)
        graphc <- dodgr::dodgr_contract_graph (graph)
        # same hash -> contracted with **NO** additional vertices
        fname_c <- file.path (td, paste0 ("graphc_", hash, ".Rds"))
        saveRDS (graphc, fname_c)

        fname_e <- paste0 ("edge_map_", hash, ".Rds")
        fname_e_fr <- file.path (tempdir (), fname_e)
        fname_e_to <- file.path (td, fname_e)
        if (file.exists (fname_e_fr)) # should always be
            file.copy (fname_e_fr, fname_e_to, overwrite = TRUE)

        fname_j <- paste0 ("junctions_", hash, ".Rds")
        fname_j_fr <- file.path (tempdir (), fname_j)
        fname_j_to <- file.path (td, fname_j)
        if (file.exists (fname_j_fr)) # should always be
            file.copy (fname_j_fr, fname_j_to, overwrite = TRUE)
    }

    sink (file = file.path (tempdir (), "Rout.txt"))
    res <- callr::r_bg (f, list (graph, hash, td))
    sink ()

    return (res) # R6 processx object
}

