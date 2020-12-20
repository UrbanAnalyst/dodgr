# nocov start
.onLoad <- function (libname, pkgname) { # nolint

    Sys.setenv ("RCPP_PARALLEL_BACKEND" = "tinythread")
    invisible ()
}
# nocov end
