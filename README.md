
<!-- README.md is generated from README.Rmd. Please edit that file -->

# templaceICAr

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/mandymejia/templateICAr.svg?branch=master)](https://travis-ci.com/mandymejia/templateICAr)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/mandymejia/templateICAr?branch=master&svg=true)](https://ci.appveyor.com/project/mandymejia/templateICAr)
[![Coveralls test
coverage](https://coveralls.io/repos/github/mandymejia/templateICAr/badge.svg)](https://coveralls.io/github/mandymejia/templateICAr)
<!-- badges: end -->

This package contains functions implementing the template ICA model
proposed in Mejia et al. (2019) and the spatial template ICA model
proposed in proposed in Mejia et al. (2020+). For both models,
subject-level brain networks are estimated as deviations from known
population-level networks, which can be estimated using standard ICA
algorithms. Both models employ an expectation-maximization algorithm for
estimation of the latent brain networks and unknown model parameters.

## Installation

You can install the development version of `templaceICAr` from Github
with:

``` r
# install.packages("devtools")
devtools::install_github("templaceICAr")
```
