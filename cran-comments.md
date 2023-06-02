# CRAN notes for dodgr_0.2.21 submission

This submission removes 'CXX_STD=CXX11' statements from all Makevars files, as requested.

This submission generates the following NOTES on some systems:

* "GNU make is a SystemRequirements", which is unavoidable because of the need to remove compiled object files in src sub-directories.
- Size Notes, which are also unavoidable as they are largely due to size of compiled "libs".

Other than these, this submission generates no additional notes, and no warnings on:

* Ubuntu 22.04: R-release, R-devel
* win-builder (R-release, R-devel, R-oldrelease)
- clang UBSAN on R-devel
