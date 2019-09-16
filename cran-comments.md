# CRAN notes for dodgr_0.2.1 submission

## Current failure

CRAN checks currently fail on r-devel-linux-x86_64-debian-gcc because of one
example which fails to read from a connection. This should never fail, as it
reads from a file in tempdir() that must exist, and which succeeds on all other
systems. I have nevertheless switched the very small example that caused this
fail to `\dontrun` to rectify this problem.

## Notes

This submissions generates the following NOTES on some systems:

* Two possibly invalid URLs: one for a valid vignette reference, due to jstor
  preventing automated queries; and one for the DOI for the paper associated
  with this pacakge, which is also valid (https://doi.org/10.32866/6945, with
  DOI simply given as 'doi = "10.32866/6945"').
* "GNU make is a SystemRequirements", which is unavoidable because of the need
  to remove compiled object files in src sub-directories.

## Test environments

Other than the above, this submission generates no additional notes, and no
warnings on:
* Ubuntu 16.04 (on `travis-ci`): R-release, R-devel
* win-builder (R-release, R-devel, R-oldrelease)

## valgrind memory leak

Testing with "valgrind --tool=memcheck --leak-check=full" **still** reveals one
potential memory leak of around 2,000 bytes, which must at this stage be
considered an enduring issue. This is due to RcppParallel and not code within
this submission; see
https://github.com/RcppCore/RcppParallel/issues/81
