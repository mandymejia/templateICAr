% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EM_diagnosticICA.R
\name{UpdateTheta_diagnosticICA}
\alias{UpdateTheta_diagnosticICA}
\alias{UpdateTheta_diagnosticICA.spatial}
\alias{UpdateTheta_diagnosticICA.independent}
\title{Parameter Estimates in EM Algorithm for Diagnostic ICA Model}
\usage{
UpdateTheta_diagnosticICA.spatial(
  template_mean,
  template_var,
  meshes,
  BOLD,
  theta,
  C_diag,
  s0_vec_list,
  D_list,
  Dinv_s0_list,
  verbose = FALSE,
  return_MAP = FALSE,
  update = c("all", "kappa", "A"),
  ignore_determinant = TRUE
)

UpdateTheta_diagnosticICA.independent(
  template_mean,
  template_var,
  template_var_max,
  BOLD,
  theta,
  C_diag,
  return_MAP = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{template_mean}{(A list of \eqn{G} matrices, each \eqn{VxL}) template
mean estimates for each group \eqn{1} to \eqn{G}, where \eqn{L} is the number of ICs, \eqn{V} is
the number of data or mesh locations.}

\item{template_var}{(A list of \eqn{G} matrices, each \eqn{VxL}) template
variance estimates for each group \eqn{1} to \eqn{G}.}

\item{meshes}{Either \code{NULL} (assume spatial independence) or a list of
objects of type \code{templateICA_mesh} created by \code{make_mesh}
(each list element corresponds to one brain structure).}

\item{BOLD}{(\eqn{VxL} matrix) dimension-reduced fMRI data.}

\item{theta}{(list) current parameter estimates}

\item{C_diag}{(\eqn{Lx1}) diagonal elements of residual covariance after dimension reduction}

\item{s0_vec_list}{List of vectorized template means (one for each group \eqn{1} to \eqn{G})}

\item{D_list}{List of sparse diagonal matrices of template standard deviations (one for each group \eqn{1} to \eqn{G})}

\item{Dinv_s0_list}{List of D^{-1} times s0_vec (one for each group \eqn{1} to \eqn{G})}

\item{verbose}{If \code{TRUE}, display progress of algorithm. Default: \code{FALSE}.}

\item{return_MAP}{If \code{TRUE}, return the posterior mean and precision of
the latent fields and group membership probabilities instead of the
parameter estimates.}

\item{update}{Which parameters to update. Either \code{"all"}, \code{"A"} or \code{"kappa"}.}

\item{ignore_determinant}{For spatial model only. If \code{TRUE} (default),
ignore the normalizing constant in \eqn{p(y\|z)} when computing posterior
probabilities of \code{z}.}

\item{template_var_max}{The maximum of the template variance across group.,}
}
\value{
An updated list of parameter estimates, theta, OR if \code{return_MAP=TRUE},
 the posterior mean and precision of the latent fields and posterior
 probabilities of group membership.
}
\description{
Parameter Estimates in EM Algorithm for Diagnostic ICA Model
}
