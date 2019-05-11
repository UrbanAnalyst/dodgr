# CRAN notes for dodgr_0.1.4 submission

## Notes

This submissions generates the following NOTES on some systems:

* One possibly invalid URL for a vignette reference, due to jstor preventing
  automated queries. The URL is nevertheless valid.
* "GNU make is a SystemRequirements", which is unavoidable because of the need
  to remove compiled object files in src sub-directories.

The submission also fails on win-builder R 3.5.3 (old-rel)with a fail on
install because "Error : package 'osmdata' was installed by an R version with
different internals; it needs to be reinstalled for use with this R version".

## Test environments

Other than the above, this submission generates no additional notes, and no
warnings on:
* Ubuntu 14.04 (on `travis-ci`): R-release, R-devel
* Windows Visual Studio 2015 (on `appveyor`; `x64`): R-release, R-devel
* win-builder (R-release, R-devel, R-oldrelease)

## valgrind memory leak

Testing with "valgrind --tool=memcheck --leak-check=full" **still** reveals one
potential memory leak of around 2,000 bytes. This is due to RcppParallel and not
code within this submission; see
https://github.com/RcppCore/RcppParallel/issues/81
