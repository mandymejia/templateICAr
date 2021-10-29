#' Scale BOLD
#' 
#' Center BOLD data across space and time, and scale
#'
#' @param BOLD (\eqn{VxT} matrix) fMRI timeseries data matrix
#' @param center Center BOLD data across space and time?
#' @param scale A logical value indicating whether the fMRI timeseries should be scaled by the image standard deviation (see Details).
#'
#' @return Centered and scaled BOLD data (\eqn{VxT} matrix)
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
scale_BOLD <- function(BOLD, center=TRUE, scale=FALSE){
  ntime <- ncol(BOLD) #length of timeseries
  nvox <- nrow(BOLD) #number of data locations
  if (ntime > nvox) warning('More time points than voxels. Are you sure?')

  if (center) {
    BOLD <- t(BOLD - rowMeans(BOLD)) #center each voxel time series (remove mean image)
    BOLD <- t(BOLD - rowMeans(BOLD)) #center each image (centering across space)
  }
  
  if (scale) { 
    sig <- sqrt(mean(rowVars(BOLD)))
    BOLD <- BOLD/sig #standardize by global SD
  } 
  
  BOLD
}
