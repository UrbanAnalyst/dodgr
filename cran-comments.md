# CRAN notes for dodgr_0.2.20 submission

This submission rectifies a previous submission which generated a note
regarding possibly invalid URL. Although that URL was valid, this submission
replaces it with another which should avoid that note.

This submission generates the following NOTES on some systems:

* "GNU make is a SystemRequirements", which is unavoidable because of the need to remove compiled object files in src sub-directories.
- Size Notes, which are also unavoidable as they are largely due to size of compiled "libs".

Other than these, this submission generates no additional notes, and no warnings on:

* Ubuntu 22.04: R-release, R-devel
* win-builder (R-release, R-devel, R-oldrelease)
- clang UBSAN on R-devel
