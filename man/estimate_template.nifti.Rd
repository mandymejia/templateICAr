% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_template.R
\name{estimate_template.nifti}
\alias{estimate_template.nifti}
\title{Estimate NIFTI template}
\usage{
estimate_template.nifti(
  nifti_fnames,
  nifti_fnames2 = NULL,
  GICA_fname,
  mask_fname,
  inds = NULL,
  scale = TRUE,
  verbose = TRUE,
  out_fname = NULL
)
}
\arguments{
\item{nifti_fnames}{Vector of file paths of NIFTI-format fMRI timeseries for template estimation.}

\item{nifti_fnames2}{(Optional) Vector of file paths of "retest" NIFTI-format fMRI
timeseries for template estimation.  Must be from the same subjects and in the same
order as nifti_fnames.  Should only be provided if nifti_fnames provided, but not required.
If none specified, will create pseudo test-retest data from single session.}

\item{GICA_fname}{File path of NIFTI-format group ICA maps (Q IC's)}

\item{mask_fname}{File path of NIFTI-format binary brain map}

\item{inds}{Indicators of which L <= Q group ICs to include in template. If NULL,
use all Q original group ICs.}

\item{scale}{Logical indicating whether BOLD data should be scaled by the
spatial standard deviation before template estimation.}

\item{verbose}{If \code{TRUE}. display progress updates}

\item{out_fname}{The path and base name prefix of the NIFTI files to write.
Will be appended with "_mean.nii" for template mean maps and "_var.nii" for
template variance maps.}
}
\value{
List of two elements: template mean of class nifti and template variance of class nifti
}
\description{
Estimate template for Template or Diagnostic ICA based on NIFTI-format data
}
