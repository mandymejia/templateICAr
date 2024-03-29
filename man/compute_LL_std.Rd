% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EM_templateICA.R
\name{compute_LL_std}
\alias{compute_LL_std}
\title{Compute LL for EM Template ICA}
\usage{
compute_LL_std(theta, template_mean, template_var, C_diag, BOLD, verbose)
}
\arguments{
\item{theta}{The current parameter estimates}

\item{template_mean, template_var}{The template}

\item{C_diag}{The C matrix}

\item{BOLD}{(\eqn{V \times T} matrix) preprocessed fMRI data}

\item{verbose}{Print LL components?}
}
\value{
The expected log posterior at the current values
}
\description{
Compute the expected log posterior for EM Standard Template ICA
}
\keyword{internal}
