# CRAN notes for dodgr_0.2.9 submission

## Previous errors

This submission fixes previous errors on some systems, due to test failures reflecting recent updated to the 'sf' package.

## Notes

This submission generates NOTES on some systems regarding "possibly invalid URLs", for which I can only assure you that all URLs are valid, and all use appropriate "https", and not "http", protocol. Note in particular that https://srtm.csi.cgiar.org/srtmdata commonly times out and errors, but the https is valid, and is the sole link to an official NASA data source. Please note that all URLs identified in 'R CMD check' NOTEs only arise in documentation, and are not directly called by any functions within the package. Thus any issues with these URLs do not affect package functionality in any way.

Beyond that, this submission generates the following NOTES on some systems:

* "GNU make is a SystemRequirements", which is unavoidable because of the need to remove compiled object files in src sub-directories.

## Test environments

Other than the above, this submission generates no additional notes, and no warnings on:
* Ubuntu 20.04: R-release, R-devel
* win-builder (R-release, R-devel, R-oldrelease)
- clang UBSAN on R-devel
