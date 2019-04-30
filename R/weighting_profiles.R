#' write_dodgr_wt_profile
#' 
#' Write the `dodgr` street network weighting profiles to a local
#' `.json`-formatted file for manual editing and subsequent re-reading.
#'
#' @param file Full name (including path) of file to which to write. The `.json`
#' suffix will be automatically appended.
#' @return TRUE if writing succussful.
#' @seealso \link{weight_streetnet}
#' @export
write_dodgr_wt_profile <- function (file = NULL)
{
    requireNamespace ("jsonlite")

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


get_profiles <- function (wt_profile, file = NULL)
{
    if (is.null (file))
        wp <- dodgr::weighting_profiles
    else
        wp <- read_dodgr_wt_profile (file)
    return (wp)
}

get_profile <- function (wt_profile, file = NULL)
{
    wp <- get_profiles (wt_profile, file)

    profiles <- wp$weighting_profiles
    prf_names <- unique (profiles$name)
    # foot, horse, wheelchair, bicycle, moped, 
    # motorcycle, motorcar, goods, hgv, psv
    wt_profile <- match.arg (tolower (wt_profile), prf_names)
    profiles [profiles$name == wt_profile, ]
}

get_surface_speeds <- function (wt_profile, file = NULL)
{
    s <- get_profiles (wt_profile, file)$surface_speeds
    s [s$name == wt_profile, ]
}

get_turn_penalties <- function (wt_profile, file = NULL)
{
    tp <- get_profiles (wt_profile, file)$penalties
    tp [tp$name == wt_profile, ]
}
