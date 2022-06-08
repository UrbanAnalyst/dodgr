# CRAN notes for dodgr_0.2.14 submission

This submission generates the following NOTES on some systems:

* "GNU make is a SystemRequirements", which is unavoidable because of the need to remove compiled object files in src sub-directories.
- Size Notes, which are also unavoidable as they are largely due to size of compiled "libs".
- One "(possibly) invalid URL" which returns Status: 403. This URL is from 'jstor', and is definitely valid.

Other than these, this submission generates no additional notes, and no warnings on:

* Ubuntu 20.04: R-release, R-devel
* win-builder (R-release, R-devel, R-oldrelease)
- clang UBSAN on R-devel
