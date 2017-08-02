#' generate_random_points
#'
#' Generate a numeric \code{matrix} of coordinates that lie within a given city.
#' The city can be specified by name.
#'
#' @param city Name of the city to be analysed
#' @param n Number of coordinates to be generated
#' @return \code{matrix} with two columns of length n with random coordinates
#'
#' @export
generate_random_points <- function (city, n)
{
    if (missing (city))
        stop ("city is missing with no default.")
    if (missing (n))
        stop ("n is missing with no default.")
    if (n < 1)
        stop ("n must be >= 1")

    poly_city <- osmdata::getbb (city, format_out = "polygon", limit = 1)
    if (any (is.na (poly_city)))
        stop (paste0 ("Could not find geometry for city '", city, "'."))

    xpoly <- poly_city [, 1]
    ypoly <- poly_city [, 2]
    xmin <- xpoly %>% min
    ymin <- ypoly %>% min
    xmax <- xpoly %>% max
    ymax <- ypoly %>% max

    n_in <- 0
    pts_out <- matrix (ncol = 2, nrow = n)
    while (n_in < n)
    {
        rnd <- random_points_in_polygon (xmin, ymin, xmax, ymax, n)
        keep <- sp::point.in.polygon (rnd [, 1], rnd [, 2], xpoly, ypoly) != 0
        i1 <- max (n_in, 0) + 1
        nkeep <- length (which (keep))
        n_in <- n_in + nkeep
        i2 <- min (n_in, n)
        rem_space <- abs (i2 - i1 + 1)
        if (nkeep > rem_space)
            nkeep <- rem_space
        pts_out [i1:i2, ] <- rnd [keep, ] [1:nkeep, ]
    }
    pts_out
}

#' random_points_in_polygon
#'
#' Generate a numeric \code{matrix} containing random coordinates within a given
#' boundary.
#'
#' @param xmin Minimum x coordinate
#' @param ymin Minimum y coordinate
#' @param xmax Maximum x coordinate
#' @param ymax Maximum y coordinate
#' @param n Number of coordinates to be generated
#' @return \code{matrix} with two columns of length n with random coordinates
#'
#' @noRd
random_points_in_polygon <- function (xmin, ymin, xmax, ymax, n)
{
    xdiff <- abs (xmax - xmin)
    ydiff <- abs (ymax - ymin)
    xrand <- xmin + (xdiff * runif (n))
    yrand <- ymin + (ydiff * runif (n))
    matrix (c (xrand, yrand), ncol = 2)
}
