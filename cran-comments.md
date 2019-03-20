# CRAN notes for dodgr_0.1.3 submission

The description now includes a DOI to the published manuscript accompanying this
package.

## Notes

This submissions generates the following NOTES on some systems:

* One possibly invalid URL for a vignette reference, due to jstor preventing
  automated queries. The URL is nevertheless valid.
* "GNU make is a SystemRequirements", which is unavoidable because of the need
  to remove compiled object files in src sub-directories.

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
