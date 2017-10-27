# CRAN notes for dodgr_0.0.2 submission

* Fixed previous failure on solaris (a C++ implicit type conversion)
* Fixed previous failure on Windows oldrel (it was just a test failure).

## Test environments

This submission generates NO notes on:
* Linux (via Travis-ci): R-release, R-devel
* OSX (via Travis-ci): R-release
* win-builder: R-oldrelease, R-release, R-devel

Package also checked using both local memory sanitzer and `rocker/r-devel-san`
with clean results. 
