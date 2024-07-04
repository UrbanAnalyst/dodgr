# CRAN notes for dodgr_0.4.1 submission

This submission fixes previous warnings from g++ 14.x which arose on CRAN Fedora checks.

It also rectifies one NOTE on namespace imports not used.

This submission generates the following NOTES on some systems:

- Possibly invalid URLs, all of which are GitHub because of "Too Many Requests," and not because of the URLs themselves.
* "GNU make is a SystemRequirements", which is unavoidable because of the need to remove compiled object files in src sub-directories.
- Size Notes, which are also unavoidable as they are largely due to size of compiled "libs".

Other than these, this submission generates no additional notes, and no warnings on:

* Ubuntu 22.04: R-release, R-devel
* win-builder (R-release, R-devel, R-oldrelease)
- clang UBSAN on R-devel
