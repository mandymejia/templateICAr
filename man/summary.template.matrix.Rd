% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_template.methods.R
\name{summary.template.matrix}
\alias{summary.template.matrix}
\alias{print.summary.template.matrix}
\alias{print.template.matrix}
\title{Summarize a \code{"template.matrix"} object}
\usage{
\method{summary}{template.matrix}(object, ...)

\method{print}{summary.template.matrix}(x, ...)

\method{print}{template.matrix}(x, ...)
}
\arguments{
\item{object}{Object of class \code{"template.matrix"}.}

\item{...}{further arguments passed to or from other methods.}

\item{x}{The template from \code{estimate_template.cifti}}
}
\value{
A list summarizing the template: data dimensions, options used for
template estimation, etc.

Nothing, invisibly.

Nothing, invisibly.
}
\description{
Summary method for class \code{"template.matrix"}
}
