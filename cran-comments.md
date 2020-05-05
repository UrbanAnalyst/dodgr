# CRAN notes for dodgr_0.2.7 submission

This submission fixes the single UBSAN error raised by previous submissions. I was able to reproduce that error locally, and can confirm that it has now been rectified.

## Notes

This submission generates the following NOTES on some systems:

* Two possibly invalid URLs: one for a valid vignette reference, due to jstor preventing automated queries; and one for the DOI for the paper associated with this pacakge, which is also valid (https://doi.org/10.32866/6945, with DOI simply given as 'doi = "10.32866/6945"').
* "GNU make is a SystemRequirements", which is unavoidable because of the need to remove compiled object files in src sub-directories.

## Test environments

Other than the above, this submission generates no additional notes, and no warnings on:
* Ubuntu 16.04 (on `travis-ci`): R-release, R-devel
* win-builder (R-release, R-devel, R-oldrelease)
- clang UBSAN on R-devel
