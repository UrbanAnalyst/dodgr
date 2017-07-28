#' get_streetnet
#'
#' @param pts List of points presumably containing spatial coordinates
#' @param expand Relative factor by which street network should extend beyond
#' limits defined by pts
#' @param txt Used for warning messages ("from" or "to")
#' @return square matrix of distances between nodes
#'
#' @export
get_streetnet <- function (pts, expand = 0.05, txt = "from")
{
    colx <- which (grepl ("x", names (pts), ignore.case = TRUE) |
                   grepl ("lon", names (pts), ignore.case = TRUE))
    coly <- which (grepl ("y", names (pts), ignore.case = TRUE) |
                   grepl ("lat", names (pts), ignore.case = TRUE))
    if (! (length (colx) == 1 | length (coly) == 1))
        stop (paste0 ("Can not unambiguous determine coordinates in ", txt))

    x <- range (pts [, colx])
    x <- x + c (-expand, expand) * diff (x)
    y <- range (pts [, coly])
    y <- y + c (-expand, expand) * diff (y)

    dat <- osmdata::opq (c (x [1], y [1], x [2], y [2])) %>%
        osmdata::add_feature (key = "highway") %>%
        osmdata::osmdata_sf ()

    return (dat)
}
