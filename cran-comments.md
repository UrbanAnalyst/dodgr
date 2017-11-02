# CRAN notes for dodgr_0.0.3 submission

This submission should rectify all previous errors on CRAN solaris machine.
These were from internally bundled code which had a mess of implicit
type conversions. Please accept my apologies for not discovering this myself,
and convey my gratitude to Brian Ripley for uncovering one of these. This
submission should recitfy all implicit type conversions. Compiling with clang++
-Weverything generates only three types of warnings about (i) non-compliance with
C++98; (ii) byte padding to align objects/boundaries; and (iii) out-of-line
virtual method definitions that are raised by class destructors only and reflect
desirable behaviour here.

## Test environments

This submission generates NO notes on:
* Linux (via Travis-ci): R-release, R-devel
* OSX (via Travis-ci): R-release
* win-builder: R-oldrelease, R-release, R-devel

Package also checked using both memory sanitzer (valgrind) and undefined behaviour
sanitizer (llvm UBSan), with clean results in both cases.
