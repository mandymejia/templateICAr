
<!-- README.md is generated from README.Rmd. Please edit that file -->

# templateICAr

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/templateICAr)](https://cran.r-project.org/package=templateICAr)
[![R-CMD-check](https://github.com/mandymejia/templateICAr/workflows/R-CMD-check/badge.svg)](https://github.com/mandymejia/templateICAr/actions)
[![Codecov test
coverage](https://codecov.io/gh/mandymejia/templateICAr/branch/master/graph/badge.svg)](https://app.codecov.io/gh/mandymejia/templateICAr?branch=master)
<!-- badges: end -->

This package contains functions implementing the template ICA model
proposed in Mejia et al. (2019) and the spatial template ICA model
proposed in proposed in Mejia et al. (2020+). For both models,
subject-level brain networks are estimated as deviations from known
population-level networks, which can be estimated using standard ICA
algorithms. Both models employ an expectation-maximization algorithm for
estimation of the latent brain networks and unknown model parameters.

Template ICA consists of three steps. The main functions associated with
each step are listed below.

1.  Template estimation: `estimate_template`. Can export the results
    with `export_template`.
2.  Template ICA model estimation (single-subject): `templateICA`.
3.  Identification of areas of engagement in each IC (or deviation from
    the template mean): `activations`.

## Citation

If you use `templateICAr` please cite the following papers:

| Name                                                                  | APA Citation                                                                                                                                                                                                                                                                                 |
|-----------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [Template ICA](https://doi.org/10.1080/01621459.2019.1679638)         | Mejia, A. F., Nebel, M. B., Wang, Y., Caffo, B. S., & Guo, Y. (2020). Template Independent Component Analysis: targeted and reliable estimation of subject-level brain networks using big data population priors. Journal of the American Statistical Association, 115(531), 1151-1177.      |
| [Spatial Template ICA](https://doi.org/10.1080/10618600.2022.2104289) | Mejia, A. F., Bolin, D., Yue, Y. R., Wang, J., Caffo, B. S., & Nebel, M. B. (2022). Template Independent Component Analysis with spatial priors for accurate subject-level brain network estimation and inference. Journal of Computational and Graphical Statistics, (just-accepted), 1-35. |

You can also obtain citation information from within R like so:

``` r
citation("templateICAr")
```

## Installation

You can install the development version of `templateICAr` from Github
with:

``` r
# install.packages("devtools")
devtools::install_github("mandymejia/templateICAr")
```

## Important Notes on Dependencies:

To analyze or visualize CIFTI-format data, `templateICAr` depends on the
`ciftiTools` package, which requires an installation of Connectome
Workbench. It can be installed from the [HCP
website](https://www.humanconnectome.org/software/get-connectome-workbench).

For fitting the template ICA model with surface-based priors
(`spatial_model=TRUE` in `templateICA()`), INLA is required. Due to a
CRAN policy, INLA cannot be installed automatically. You can obtain it
by running
`install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)`.
Alternatively, `dep=FALSE` can be used along with manual installation of
dependencies as necessary to avoid installing all of the many INLA
dependencies, most of which are not actually required. Binaries for
alternative Linux builds can be added with the command
`inla.binary.install()`. Note that INLA is *not* required for standard
template ICA.

Depending on the analysis, PARDISO may reduce computation time. To
obtain a free academic license forINLA-PARDISO, run `inla.pardiso()` in
R after running `library(INLA)`. Provide an academic email address. Once
you obtain a license, point to it using
`INLA::inla.setOption(pardiso.license = "pardiso.lic")` followed by
`INLA::inla.pardiso.check()` to ensure that PARDISO is successfully
installed and running.
