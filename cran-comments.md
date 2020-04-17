# CRAN notes for dodgr_0.2.6 submission

This submission fixes the current warnings that arose due to the recent tibble upgrade. The submission should also finally fix the previously reported valgrind memory leaks, due to Intel's TBB library bundled with RcppParallel. The authors of that package were unable to find a solution, and so recommend switching to a different library for multi-thread computation, which this submission now does.

## Notes

This submission generates the following NOTES on some systems:

* Two possibly invalid URLs: one for a valid vignette reference, due to jstor preventing automated queries; and one for the DOI for the paper associated with this pacakge, which is also valid (https://doi.org/10.32866/6945, with DOI simply given as 'doi = "10.32866/6945"').
* "GNU make is a SystemRequirements", which is unavoidable because of the need to remove compiled object files in src sub-directories.

## Test environments

Other than the above, this submission generates no additional notes, and no warnings on:
* Ubuntu 16.04 (on `travis-ci`): R-release, R-devel
* win-builder (R-release, R-devel, R-oldrelease)
