## Version 1.2.0 July 29, 2025

This updated contains modifications in line with the final published version of the scGBM article in *Biostatistcs*.

**New functionality:**

- `gbm.sc()` now includes an argument `sigma`, that controls the prior mean of the singular values (these are assumed to follow an exponential prior).
- `gbm.sc()` also includes arguments `oos.Y` for out-of-sample likelihood calculation, `order.by.deviance` for finding the best ordering of factors, and `factor.init` for flexibility in how to initialize the `V`. 

## Version 1.1.1 February 8, 2024

This update contains minor changes to the code structure and also no longer emphasizes the use of the coarse cluster feature of the `CCI` function.

## Version 1.1.0 September 10, 2023

This update contains new functionality related to interpreting the loadings and also some updates to the `CCI` function.

**New functionality:**

-   `loadings.volcano()` makes a volcano plot to show the genes that are driving a particular latent factor.
-   Several changes to the `CCI` function, including a function that automatically combines clusters with low inter-CCI and, if `null.dist = TRUE`, a cutoff line corresponding to the CCI that would be expected under the null of no latent variability.

**Changes to existing functionality:**

-   `get.se()` now has the option to add a small diagonal matrix with entries `EPS` in cases where the matrix inversion fails.

**Bug fixes:**

-   Use `scores` instead of `V` in `get.se()`.
-   Use `loadings` in `denoise.U()`.

**Other changes:**

-   Added input checks for `gbm.sc`.

## Version 1.0.1 June 28, 2023

This patch contains minor bug fixes and changes the names of some objects.

**Changes to existing functionality:**

-   `gbm.sc()` now returns `scores` instead of `V` to be more consistent with the model defined in the paper.

**Bug fixes:**

-   For `gbm.sc`, the case `M=1` now works.

**Other changes:**

-   Improved documentation for `gbm.sc`.

## Version 1.0.0 May 25, 2023

Initial release
