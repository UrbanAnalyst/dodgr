# CRAN notes for dodgr_0.2.3 submission

## Current ASAN/UBSAN notes

This submission *should* rectify previous notes from address & undefined
behaviour sanitizers. I have nevertheless been unable to confirm this as I have
not been able to reproduce the result. I tried the rocker/r-devel-ubsan-clang
container, but that failed to install package due to the valgrind memory leak
from TBB via RcppParallel mentioned below. I also checked on r-hub's sanitizer
system, but that failed for same reason, and not for reasons related to current
CRAN failures. I could also not reproduce locally. One failing test I did
manage to reproduce was the valgrind test, which this submission definitely
fixes, reducing the possibly lost byte count from >200kB back to the "usual"
2kB due to RcppParallel / TBB. I confidently presume that the fix also
successfully addresses the ASAN/UBSAN issues. Sorry for any inconvenience.

I have also rectified one previous, intermittently failing test.

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
