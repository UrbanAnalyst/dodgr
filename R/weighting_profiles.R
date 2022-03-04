#' write_dodgr_wt_profile
#'
#' Write the `dodgr` street network weighting profiles to a local
#' `.json`-formatted file for manual editing and subsequent re-reading.
#'
#' @param file Full name (including path) of file to which to write. The `.json`
#' suffix will be automatically appended.
#' @return TRUE if writing successful.
#' @seealso \link{weight_streetnet}
#' @family misc
#' @export
write_dodgr_wt_profile <- function (file = NULL) {

    requireNamespace ("jsonlite")

    if (is.null (file))
        stop ("file name must be given")

    file <- paste0 (tools::file_path_sans_ext (file), ".json")
    con <- file (file, open = "wt")

    sc <- summary (con)
    if (!sc [["can write"]] == "yes")
        stop ("Unable to write to connection ", sc$description) # nocov

    wpj <- jsonlite::toJSON (dodgr::weighting_profiles, pretty = TRUE)
    writeLines (wpj, con)
    close (con)
}

read_dodgr_wt_profile <- function (file = NULL) {

    file <- paste0 (tools::file_path_sans_ext (file), ".json")
    if (!fs::file_exists (file))
        stop ("file [", file, "] does not exist") # nocov

    res <- jsonlite::fromJSON (file)
    # jsonlite interprets these as "integer":
    storage.mode (res$surface_speeds$max_speed) <- "numeric"
    storage.mode (res$penalties$traffic_lights) <- "numeric"
    return (res)
}


get_profiles <- function (file = NULL) {

    if (is.null (file))
        wp <- dodgr::weighting_profiles
    else
        wp <- read_dodgr_wt_profile (file)
    return (wp)
}

get_profile <- function (wt_profile, file = NULL) {

    profiles <- get_profiles (file)$weighting_profiles
    prf_names <- unique (profiles$name)
    if (is.numeric (wt_profile)) {
        # nocov start
        # this function is actually only called for character args
        wp <- profiles [profiles$name == "foot", ]
        wp$name <- "custom"
        wp$value <- wt_profile
        wp$max_speed <- 10
        # nocov end
    } else {
        # foot, horse, wheelchair, bicycle, moped,
        # motorcycle, motorcar, goods, hgv, psv
        wt_profile <- match.arg (tolower (wt_profile), prf_names)
        wp <- profiles [profiles$name == wt_profile, ]
    }
    return (wp)
}

get_surface_speeds <- function (wt_profile, file = NULL) {

    s <- get_profiles (file)$surface_speeds
    s [s$name == wt_profile, ]
}

get_turn_penalties <- function (wt_profile, file = NULL) {

    tp <- get_profiles (file)$penalties
    tp [tp$name == wt_profile, ]
}

are_turns_restricted <- function (wt_profile, file = NULL) {

    tp <- get_profiles (file)$penalties
    tp$restrictions [tp$name == wt_profile]
}
