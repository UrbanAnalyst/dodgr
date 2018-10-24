# CRAN notes for dodgr_0.1.1 submission

We are currently unable to provide a reference in the package Description, but
hope to have a published manuscript soon, and will update asap thereafter.

## Notes

This submissions generates the following NOTES on some systems:

* One possibly invalid URL for a vignette reference, due to jstor preventing
  automated queries. The URL is nevertheless valid.
* Package size slightly only 5MB, because of large doc and libs directories. I
  will endeavour to reduce this in future versions.
* "GNU make is a SystemRequirements", which is unavoidable because of the need
  to remove compiled object files in src sub-directories.

## Test environments

Other than the above, this submission generates no additional notes, and no
warnings on:
* Ubuntu 14.04 (on `travis-ci`): R-release, R-devel
* Windows Visual Studio 2015 (on `appveyor`; `x64`): R-release, R-devel
* win-builder (R-release, R-devel, R-oldrelease)

## valgrind memory leak

Testing with "valgrind --tool=memcheck --leak-check=full" reveals one potential
memory leak of around 2,000 bytes. This is due to RcppParallel and not code
within this submission; see https://github.com/RcppCore/RcppParallel/issues/81
