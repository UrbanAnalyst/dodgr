# CRAN notes for dodgr_0.0.3 submission

This submission should rectify all previous errors on CRAN solaris machine.
These were from internally bundled code and were only flagged on that machine.
Clang with all warning options on did not raise these errors, making it
difficult to identify and remove them. I am confident that this submission has
now achieved that.

## Test environments

This submission generates NO notes on:
* Linux (via Travis-ci): R-release, R-devel
* OSX (via Travis-ci): R-release
* win-builder: R-oldrelease, R-release, R-devel

Package also checked using both local memory sanitzer and `rocker/r-devel-san`
with clean results. 
