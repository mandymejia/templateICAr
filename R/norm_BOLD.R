#' Normalize BOLD data
#' 
#' Center the data across space and time, scale, and detrend, in that order.
#'
#' @param BOLD fMRI numeric data matrix (\eqn{V \times T})
#' @param center_rows,center_cols Center BOLD data across rows (each data location's time series) or columns (each time point's image)? Default: \code{TRUE} for both.
#' @param scale A logical value indicating whether the fMRI timeseries should be scaled by the image standard deviation (see Details).
#' @param detrend_DCT Detrend the data? This is the number of DCT bases to use for detrending. If \code{0} (default), do not detrend.
#' 
#' @return Normalized BOLD data matrix (\eqn{V \times T})
#' 
#' @export
#'
#' @details The BOLD data is centered across time, which removes the mean image by subtracting the mean of each voxel timeseries, and across space,
#' which removes the global signal.  Both steps are necessary for dual regression, which in the first regression treats time points as observations,
#' and in the second step treats voxels as observations. Neither regression includes an intercept, so the BOLD timeseries must be centered across
#' both space and time. The group independent components used in the first regression must also be centered across space.  The mixing matrix estimated
#' in the first step of dual regression is also centered across time as a result.
#'
#' In order to ensure that all fMRI data share the same scale, the BOLD data can also be scaled by the global image standard deviation, equal to
#' \deqn{\sqrt{\frac{1}{T}\sum_{t=1}^T \sigma^2_t}},
#' where \eqn{\sigma^2_t} is the standard deviation across all voxels at time point \eqn{t}. If scaling is applied to the BOLD timeseries used in
#' template estimation, it should also be applied to the BOLD timeseries being analyzed with template ICA using the resulting templates to ensure
#' compatibility of scale.
#' 
#' @importFrom fMRIscrub nuisance_regression dct_bases
#'
norm_BOLD <- function(BOLD, center_rows=TRUE, center_cols=TRUE, scale=FALSE, detrend_DCT=0){
  nT <- ncol(BOLD)
  nV <- nrow(BOLD)
  if (nT > nV) { warning('More time points than voxels. Are you sure?') }

  stopifnot(is.logical(center_rows) && length(center_rows)==1)
  stopifnot(is.logical(center_cols) && length(center_cols)==1)
  stopifnot(is.logical(scale) && length(scale)==1)
  stopifnot(is.numeric(detrend_DCT) && length(detrend_DCT)==1)
  stopifnot(detrend_DCT >=0 && detrend_DCT==round(detrend_DCT))

  # Center. 
  voxMeans <- NULL # for detrending without centering
  if (center_rows || center_cols) {
    # `BOLD` is transposed twice.
    # Center each voxel time series (across time).
    if (center_rows) {
      BOLD <- t(BOLD - rowMeans(BOLD)) 
    } else {
      if (detrend_DCT > 0) { voxMeans <- rowMeans(BOLD) }
      BOLD <- t(BOLD)
    }
    # Center each image (across space).
    if (center_cols) {
      BOLD <- t(BOLD - rowMeans(BOLD)) 
    } else {
      BOLD <- t(BOLD)
    }
  } else {
    if (detrend_DCT > 0) { voxMeans <- rowMeans(BOLD) }
  }

  # Detrend.
  if (detrend_DCT > 0) {
    BOLD <- nuisance_regression(BOLD, cbind(1, dct_bases(nT, detrend_DCT)))
    if (!center_rows) { BOLD <- BOLD + voxMeans }
  }

  # Scale by global SD.
  if (scale) { 
    sig <- sqrt(mean(rowVars(BOLD)))
    if (sig < 1e-8) {
      warning("Estimated scale is near zero. Skipping scaling.")
    } else {
      BOLD <- BOLD / sig
    }
  } 
  
  BOLD
}
