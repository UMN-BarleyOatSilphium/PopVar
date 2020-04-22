
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PopVar

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/PopVar)](https://cran.r-project.org/package=PopVar)

<!-- badges: end -->

The main attribute of ‘PopVar’ is the prediction of genetic variance in
bi-parental populations, from which the package derives its name.
‘PopVar’ contains a set of functions that use phenotypic and genotypic
data from a set of candidate parents to 1) predict the mean, genetic
variance, and superior progeny value of all, or a defined set of
pairwise bi-parental crosses, and 2) perform cross-validation to
estimate genome-wide prediction accuracy of multiple statistical models.
More details are available in Mohammadi, Tiede, and Smith (2015). Crop
Sci. <doi:10.2135/cropsci2015.01.0030>. A dataset ‘think\_barley.rda’ is
included for reference and examples.

## Installation

You can install the released version of PopVar from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("PopVar")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("UMN-BarleyOatSilphium/PopVar")
```

## Examples

Below are some example uses of the functions in `PopVar`:

``` r
# Load the package
library(PopVar)
```

### Original functions

### Predictions using deterministic equations

### Mulit-parent populations
