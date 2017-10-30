# CRAN notes for dodgr_0.0.3 submission

This submission should rectify all previous errors on CRAN solaris machine.
These were from internally bundled code which had a mess of implicit
type conversions. This submission should recitfy all of these, and compiles
clean on Clang++ using "-Wconversion".

## Test environments

This submission generates NO notes on:
* Linux (via Travis-ci): R-release, R-devel
* OSX (via Travis-ci): R-release
* win-builder: R-oldrelease, R-release, R-devel

Package also checked using both memory sanitzer (valgrind) and undefined behaviour
sanitizer (llvm UBSan), with clean results in both cases.
