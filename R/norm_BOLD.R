#' Normalize BOLD data
#'
#' Center the data across space and/or time, scale, and detrend, in that order.
#'  For dual regression, row centering is required and column centering is not
#'  recommended. Scaling and detrending depend on the user preference.
#'
#' @param BOLD fMRI numeric data matrix (\eqn{V \times T})
#' @param center_rows,center_cols Center BOLD data across rows (each data
#'  location's time series) or columns (each time point's image)? Default:
#'  \code{TRUE} for row centering, and \code{FALSE} for column centering.
#' @param scale \code{"global"} (default), \code{"local"}, or \code{"none"}.
#'  Global scaling will divide the entire data matrix by the image standard
#'  deviation (\code{sqrt(mean(rowVars(BOLD)))}). Local scaling will divide each
#'  data location's time series by its estimated standard deviation.
#' @param scale_sm_xifti,scale_sm_FWHM Only applies if \code{scale=="local"}. To
#'  smooth the standard deviation estimates used for local scaling, provide a
#'  \code{"xifti"} object with data locations in alignment with
#'  \code{"BOLD"} and the smoothing FWHM (default: \code{2}). If no \code{"xifti"}
#'  object is provided (default) or if \code{scale_sm_FWHM=0}, do not smooth.
#' @param detrend_DCT Detrend the data? This is the number of DCT bases to use
#'  for detrending. If \code{0} (default), do not detrend.
#'
#' @return Normalized BOLD data matrix (\eqn{V \times T})
#'
#' @export
#'
#' @details In order to ensure that all fMRI data share the same scale, the BOLD
#'  data can also be scaled by the global image standard deviation, equal to
#'  \deqn{\sqrt{\frac{1}{T}\sum_{t=1}^T \sigma^2_t}}, where \eqn{\sigma^2_t} is
#'  the standard deviation across all voxels at time point \eqn{t}. If scaling
#'  is applied to the BOLD timeseries used in template estimation, it should
#'  also be applied to the BOLD timeseries being analyzed with template ICA
#'  using the resulting templates to ensure compatibility of scale. The scale is
#'  computed after detrending.
#'
#' @importFrom ciftiTools is.xifti
#' @importFrom fMRIscrub nuisance_regression dct_bases
#'
norm_BOLD <- function(
  BOLD, center_rows=TRUE, center_cols=FALSE,
  scale=c("global", "local", "none"), scale_sm_xifti=NULL, scale_sm_FWHM=2,
  detrend_DCT=0){

  nT <- ncol(BOLD)
  nV <- nrow(BOLD)
  if (nT > nV) { warning('More time points than voxels. Are you sure?') }

  stopifnot(is.logical(center_rows) && length(center_rows)==1)
  stopifnot(is.logical(center_cols) && length(center_cols)==1)
  if (is.null(scale) || isFALSE(scale)) { scale <- "none" }
  if (isTRUE(scale)) {
    warning(
      "Setting `scale='global'`. Use `'global'` or `'local'` ",
      "instead of `TRUE`, which has been deprecated."
    )
    scale <- "global"
  }
  scale <- match.arg(scale, c("global", "local", "none"))
  if (!is.null(scale_sm_xifti)) { stopifnot(is.xifti(scale_sm_xifti)) }
  stopifnot(is.numeric(scale_sm_FWHM) && length(scale_sm_FWHM)==1)
  if (isFALSE(detrend_DCT)) { detrend_DCT <- 0 }
  stopifnot(is.numeric(detrend_DCT) && length(detrend_DCT)==1)
  stopifnot(detrend_DCT >=0 && detrend_DCT==round(detrend_DCT))

  # Center.
  voxMeans <- NULL # for detrending without centering
  if (center_rows || center_cols) {
    # `BOLD` is transposed twice.
    # Center each voxel time series (across time).
    if (center_rows) {
      BOLD <- t(BOLD - rowMeans(BOLD, na.rm=TRUE))
    } else {
      if (detrend_DCT > 0) { voxMeans <- rowMeans(BOLD, na.rm=TRUE) }
      BOLD <- t(BOLD)
    }
    # Center each image (across space).
    if (center_cols) {
      BOLD <- t(BOLD - rowMeans(BOLD, na.rm=TRUE))
    } else {
      BOLD <- t(BOLD)
    }
  } else {
    if (detrend_DCT > 0) { voxMeans <- rowMeans(BOLD, na.rm=TRUE) }
  }

  # Detrend.
  if (detrend_DCT > 0) {
    BOLD <- nuisance_regression(BOLD, cbind(1, dct_bases(nT, detrend_DCT)))
    if (!center_rows) { BOLD <- BOLD + voxMeans }
  }

  # Scale by global or local SD.
  if (scale == "global") {
    sig <- rowVars(BOLD, na.rm=TRUE)
    sig <- sqrt(mean(sig, na.rm=TRUE))
    if (sig < 1e-8) {
      warning("Estimated scale is near zero. Skipping scaling.")
    } else {
      BOLD <- BOLD / sig
    }
  } else if (scale == "local") {
    if (is.null(scale_sm_xifti) || (scale_sm_FWHM==0)) {
      BOLD <- BOLD / sqrt(rowVars(BOLD, na.rm=TRUE))
    } else {
      # Check `scale_sm_xifti` is valid.
      if (nV != nrow(scale_sm_xifti)) {
        stop("`scale_sm_xifti` not compatible with `BOLD`: different spatial dimensions.")
      }
      if (!is.null(scale_sm_xifti$meta$cifti$intent) && scale_sm_xifti$meta$cifti$intent == 3007) {
        scale_sm_xifti <- convert_xifti(scale_sm_xifti, "dscalar")
      }
      # Compute and smooth the SD.
      sig <- sqrt(rowVars(BOLD, na.rm=TRUE))
      sig <- newdata_xifti(select_xifti(scale_sm_xifti, 1), sig)
      sig <- smooth_xifti(sig, surf_FWHM=scale_sm_FWHM, vol_FWHM=scale_sm_FWHM)
      # Apply the local SD scaling.
      sig <- c(as.matrix(sig))
      BOLD <- BOLD / sig
    }
  }

  BOLD
}

#' Scale BOLD (legacy version of \code{norm_BOLD}) that centers both ways.
#'
#' @param BOLD fMRI numeric data matrix (\eqn{V \times T})
#' @param scale A logical value indicating whether the fMRI timeseries should
#'  be scaled by the image standard deviation. Default: \code{FALSE}.

scale_BOLD <- function(BOLD, scale=FALSE){
  warning(
    "`scale_BOLD` has been renamed to `norm_BOLD`. ",
    " `scale_BOLD` will be removed in a future version. ",
    "Please replace instances of `scale_BOLD` with `norm_BOLD`. ",
    "Also, refer to the documentation of `norm_BOLD` to see the changes."
  )

  if (isFALSE(scale)) { scale <- "none" }
  norm_BOLD(BOLD, center_cols=TRUE, scale=scale)
}
