% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/groupICA.R
\name{groupICA.cifti}
\alias{groupICA.cifti}
\title{Perform group ICA based on CIFTI data}
\usage{
groupICA.cifti(
  cifti_fnames,
  subjects = NULL,
  brainstructures = "all",
  num_PCs = 100,
  num_ICs,
  max_rows_GPCA = 10000,
  GSR = FALSE,
  scale = c("local", "global", "none"),
  scale_sm_FWHM = 2,
  TR = NULL,
  hpf = 0.01,
  verbose = TRUE,
  out_fname = NULL,
  surfL = NULL,
  surfR = NULL,
  smooth = TRUE,
  smooth_surf_FWHM = 5,
  smooth_vol_FWHM = 5,
  smooth_zeroes_as_NA = FALSE,
  smooth_subcortical_merged = FALSE
)
}
\arguments{
\item{cifti_fnames}{Vector of file paths of CIFTI-format fMRI timeseries
(*.dtseries.nii) for which to calculate group ICA}

\item{subjects}{Use this argument if some subjects have more than one scan.
This is a numeric or factor vector the same length as \code{cifti_fnames}. Scans
for the same subject should have the same value. For example, if there are four
subjects with two visits each, and the scans are ordered with the first subject's
two scans first, then the second subject's two scans next, etc., then this argument
should be \code{rep(seq(4), each=2)}. If there are three subjects and four scans,
with the last two scans belonging to the same subject, this argument should be
\code{c(1,2,3,3)}. If all subjects only have one scan, keep this argument as
\code{NULL} (default).}

\item{brainstructures}{Character vector indicating which brain structure(s)
to obtain: \code{"left"} (left cortical surface), \code{"right"} (right
cortical surface) and/or \code{"subcortical"} (subcortical and cerebellar
gray matter). Can also be \code{"all"} (obtain all three brain structures).
Default: \code{c("all")}.}

\item{num_PCs}{Maximum number of PCs to retain in initial subject-level PCA}

\item{num_ICs}{Number of independent components to identify.}

\item{max_rows_GPCA}{The maximum number of rows for the matrix on which group
level PCA will be performed.  This matrix is the result of temporal concatenation
and contains \code{length(cifti_fnames) * num_PCs} rows.}

\item{GSR, scale, scale_sm_FWHM, TR, hpf}{Center BOLD columns, scale by the
standard deviation, and detrend voxel timecourses? See
\code{\link{norm_BOLD}}. Normalization is applied separately to each scan.
Defaults: Center BOLD columns, scale by the global standard deviation, and
apply a .01 Hz HPF if the \code{TR} is provided.}

\item{verbose}{If \code{TRUE}, display progress updates}

\item{out_fname}{(Optional) If not specified, a xifti object will be returned but
the GICA maps will not be written to a CIFTI file.}

\item{surfL}{(Optional) File path to a surface GIFTI for the left cortex.
If provided, will be part of xifti result object for visualization in R.
Will also be used to perform any smoothing.}

\item{surfR}{(Optional) File path to a surface GIFTI for the right cortex.
If provided, will be part of xifti result object for visualization in R.
Will also be used to perform any smoothing.}

\item{smooth}{Smooth the CIFTI data prior to reading in each file? Default:
\code{TRUE}. Use the following arguments to control the smoothing parameters.}

\item{smooth_surf_FWHM, smooth_vol_FWHM, smooth_zeroes_as_NA, smooth_subcortical_merged}{See \code{\link[ciftiTools]{smooth_cifti}}. The defaults here are the same.
Note that \code{smooth_zeroes_as_NA} will control both of the corresponding
surface and volume arguments to \code{\link[ciftiTools]{smooth_cifti}}.
These arguments only apply if \code{smooth}.}
}
\value{
\code{out_fname} if a file was written, or the GICA as a \code{"xifti"} object
if not.
}
\description{
Perform group ICA based on CIFTI data
}
