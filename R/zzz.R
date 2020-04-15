# nocov start
.onLoad <- function (libname, pkgname)
{
    Sys.setenv ("RCPP_PARALLEL_BACKEND" = "tinythread")
    invisible ()
}
# nocov end

