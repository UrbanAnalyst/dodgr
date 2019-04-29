#' write_dodgr_wt_profile
#' 
#' Write the `dodgr` street network weighting profiles to a local
#' `.json`-formatted file for manual editing and subsequent re-reading.
#'
#' @param file Full name (including path) of file to which to write. The `.json`
#' suffix will be automatically appended.
#' @return TRUE if writing succussful.
#' @export
write_dodgr_wt_profile <- function (file = NULL)
{
    if (is.null (file))
        stop ("file name must be given")

    file <- paste0 (tools::file_path_sans_ext (file), ".json")
    con <- file (file, open = "wt")

    sc <- summary (con)
    d <- dirname (sc$description)
    if (!sc [["can write"]] == "yes")
        stop ("Unable to write to connection ", sc$description) # nocov

    wpj <- jsonlite::toJSON (dodgr::weighting_profiles, pretty = TRUE)
    writeLines (wpj, con)
    close (con)
}

read_dodgr_wt_profile <- function (file = NULL)
{
    file <- paste0 (tools::file_path_sans_ext (file), ".json")
    if (!file.exists (file))
        stop ("file [", file, "] does not exist") # nocov

    res <- jsonlite::fromJSON (file)
    # jsonlite interprets these as "integer":
    storage.mode (res$surface_speeds$max_speed) <- "numeric"
    storage.mode (res$penalties$traffic_lights) <- "numeric"
    return (res)
}
