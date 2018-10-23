# CRAN notes for dodgr_0.1.1 submission

The only NOTE generated is "GNU make is a SystemRequirements", which is
unavoidable because of the need to remove compiled object files in src
sub-directories.

## Test environments

This submission generates NO notes on:
* Ubuntu 14.04 (on `travis-ci`): R-release, R-devel
* Windows Visual Studio 2015 (on `appveyor`; `x64`): R-release, R-devel
* win-builder (R-release, R-devel, R-oldrelease)
* Package also checked using `rocker/r-devel-san` with clean results.
* local valgrind --memcheck gives clean results

Compiling with clang++ -Weverything generates only three types of warnings about
(i) non-compliance with C++98; (ii) byte padding to align objects/boundaries;
and (iii) out-of-line virtual method definitions that are raised by class
destructors only and reflect desirable behaviour here.
