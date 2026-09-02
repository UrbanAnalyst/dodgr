# CRAN notes for dodgr_0.5.0 submission

This submission fixes the compilation failure seen with current CRAN version on one machine.

Beyond that, this submission generates the following NOTES on some systems:

- Possibly invalid URLs, all of which are GitHub because of "Too Many Requests," and not because of the URLs themselves.
* "GNU make is a SystemRequirements", which is unavoidable because of the need to remove compiled object files in src sub-directories.
- Size Notes, which are also unavoidable as they are largely due to size of compiled "libs".

Other than these, this submission generates no additional notes, and no warnings on:

* Ubuntu: R-release, R-devel
* win-builder (R-release, R-devel, R-oldrelease)
- clang UBSAN on R-devel
