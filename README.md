# SigDetect

SigDetect is an R package for signal detection and signal-proportion testing in one-dimensional data built upon the theory of `Compensator Based Signal Detection'. The package provides methods for analyzing both binned and unbinned data, with or without a background-only sample, and includes utilities for sensitivity analysis and signal-region calculation (in the absence of a background-only sample). For further methodological details refer to the papers https://arxiv.org/abs/2605.20508 (for the unbinned data regime) and https://arxiv.org/abs/2607.06939 (for the binned data regime).

## Features

- Signal proportion testing for unbinned data
- Signal proportion testing for binned data
- Analyses with a background-only sample
- Analyses without a background-only sample
- Sensitivity analysis for diffused signal alternatives
- Signal-region construction based on a cumulative distribution function of the signal distribution

## Package structure

- `R/` — core implementation files
- `man/` — function documentation generated via roxygen
- `DESCRIPTION` — package metadata and dependencies
- `NAMESPACE` — exported functions

## Installation

This package is not yet published to CRAN, so install it directly from GitHub:

```r
install.packages("remotes")
remotes::install_github("baner175/SigDetect")
```

If needed, install the package dependency manually:

```r
install.packages("truncdist")
```

## Main functions

The package exports the following core functions:

- `sig_detect_unbinned_wbkg()`
- `sig_detect_unbinned_wobkg()`
- `sig_detect_binned_wbkg()`
- `sig_detect_binned_wobkg()`
- `sensitivity_analysis()`
- `signal_region()`

These are designed to support signal-detection workflows in both the settings where the background-only data is either accessible or it is not.

## Typical workflow

A typical analysis follows this pattern:

#### When a background-only data is available
1. Define the signal density. The signal distribution should not depend on unknown parameters.
2. Specify the proposal background density. This may or may not depend on some unknown parameters.
3. Define the search region.
4. Run the appropriate detection function:
   i. sig_detect_unbinned_wbkg() if the data is unbinned
   ii. sig_detect_binned_wbkg() if the data is binned
7. Interpret the test statistic, p-value, confidence interval, and estimated signal proportion.

#### When a background-only data is available
1. Define the signal density and cumulative distribution function. The signal distribution should not depend on unknown parameters.
2. Specify the proposal background density. This must depend on some unknown parameter.
3. Define the search region.
4. Run the appropriate detection function:
   i. sig_detect_unbinned_wobkg() if the data is unbinned
   ii. sig_detect_binned_wobkg() if the data is binned
7. Interpret the test statistic, p-value, confidence interval, and estimated signal proportion.

## Example usage

Below is a minimal example skeleton. You would replace the placeholders with your actual signal and background model functions.

```r
library(SigDetect)

# Example signal density
fs <- function(x) {
  # e.g. Gaussian density
  dnorm(x, mean = 0, sd = 1)
}

# Example background model: density or non-density function with x and alpha arguments
bench.bkg <- function(x, alpha) {
  # placeholder example
  dnorm(x, mean = 0, sd = alpha)
}

# Search region containing the signal
search_region <- c(-5, 5)

# Example: test using unbinned data with background-only sample
res <- sig_detect_unbinned_wbkg(
  phys.samp = rnorm(200, mean = 0.3, sd = 1),
  bkg.samp = rnorm(500, mean = 0, sd = 1),
  search_region = search_region,
  fs = fs,
  prop.bkg = bench.bkg,
  eta0 = 0,
  is.density = TRUE,
  beta_lower = 0.1,
  beta_upper = 5,
  beta_start = 1,
  conf.level = 0.95
)

print(res)
```

You can also compute a signal region:

```r
Fs <- function(x) pnorm(x, mean = 0, sd = 1)
mu_sig <- 0

signal_window <- signal_region(Fs = Fs, mu_sig = mu_sig,
                              search_region = c(-5, 5), eps = 0.001)

print(signal_window)
```

## Sensitivity analysis

The package also includes a sensitivity-analysis routine useful for studying how results change under a diffused signal bump added to the background model. This is especially useful for assessing robustness of the detection procedure to alternative signal hypotheses.

```r
sens <- sensitivity_analysis(
  bench.bkg = bench.bkg,
  search_region = c(-5, 5),
  Fs = Fs,
  mu_sig = 0,
  is.density = TRUE,
  phys.data = rnorm(200, mean = 0.3, sd = 1),
  is.binned = FALSE,
  lambda = c(0.05, 0.10, 0.20),
  sig0 = 1
)

print(sens)
```

## Notes

- The package is implemented as an R package with roxygen-generated documentation.
- Many functions expect a model specification for the background and signal distributions, along with optimization bounds and starting values.
- Some routines require additional derivative information (`d_log_bkg`, `d2_log_bkg`) when the background model has parameters.

## License

This repository currently includes package metadata indicating a typical open-source license choice; please refer to the license file or project configuration for the exact terms.

## Contributions

Contributions and feedback are welcome. If you use this package for analysis or extend its functionality, consider sharing examples, bug fixes, or documentation improvements.

## Contact

Project maintainer:
- Aritra Banerjee
- Email: baner175@umn.edu

