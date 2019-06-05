# CRAN notes for dodgr_0.2.0 submission

## Notes

This submissions generates the following NOTES on some systems:

* One possibly invalid URL for a vignette reference, due to jstor preventing
  automated queries. The URL is nevertheless valid.
* "GNU make is a SystemRequirements", which is unavoidable because of the need
  to remove compiled object files in src sub-directories.

## Test environments

Other than the above, this submission generates no additional notes, and no
warnings on:
* Ubuntu 16.04 (on `travis-ci`): R-release, R-devel
* Windows Visual Studio 2015 (on `appveyor`; `x64`): R-release, R-devel
* win-builder (R-release, R-devel, R-oldrelease)

## valgrind memory leak

Testing with "valgrind --tool=memcheck --leak-check=full" **still** reveals one
potential memory leak of around 2,000 bytes, which must at this stage be
considered an enduring issue. This is due to RcppParallel and not code within
this submission; see
https://github.com/RcppCore/RcppParallel/issues/81
