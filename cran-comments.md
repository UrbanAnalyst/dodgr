# CRAN notes for dodgr_0.2.11 submission


## Notes

This submission generates the following NOTES on some systems:

* "GNU make is a SystemRequirements", which is unavoidable because of the need to remove compiled object files in src sub-directories.
- Size Notes, which are also unavoidable as they are largely due to size of compiled "libs".

## Test environments

Other than the above, this submission generates no additional notes, and no warnings on:

* Ubuntu 20.04: R-release, R-devel
* win-builder (R-release, R-devel, R-oldrelease)
- clang UBSAN on R-devel
