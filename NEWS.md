## Version TBD Date TBD 

TODO: SUMMARY OF CHANGES
**Changes to existing functionality:**

   * `get.se()` now has the option to add a small diagonal matrix with entries `EPS` in cases where the matrix inversion fails. 

   * `gbm.sc()` OPTIONS FOR TYPE OF SVD 

**Bug fixes:**

   * Use `scores` instead of `V` in `get.se()`. 
   * Use `loadings` in `denoise.U()`. 

**Other changes:**

   * Added input checks for `gbm.sc`. 

## Version 1.0.1 June 28, 2023

This patch contains minor bug fixes and changes the names of some objects. 

**Changes to existing functionality:**

   * `gbm.sc()` now returns `scores` instead of `V` to be more consistent with the model defined in the paper.

**Bug fixes:**

   * For `gbm.sc`, the case `M=1` now works. 

**Other changes:**

   * Improved documentation for `gbm.sc`. 

## Version 1.0.0 May 25, 2023

Initial release 
