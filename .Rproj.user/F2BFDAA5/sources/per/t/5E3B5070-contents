#' Center BOLD data across space and time and scale by global temporal standard deviation
#'
#' @param BOLD (TxV matrix) fMRI timeseries data matrix
#'
#' @return Centered and scaled BOLD data (TxV matrix)
#' @export
#'
scale_BOLD <- function(BOLD){

  dat <- BOLD
  ntime <- nrow(dat) #length of timeseries
  nvox <- ncol(dat) #number of data locations

  if(ntime > nvox) warning('More time points than voxels. Are you sure?')

  #center timeseries data across space and time and standardize scale
  dat <- scale(dat, scale=FALSE) #center each voxel time series (remove mean image)
  dat_t <- t(dat) #transpose image matrix
  sig <- sqrt(mean(colVars(dat_t))) #variance across image, averaged across time, square root to get SD
  dat <- t(scale(dat_t, scale=FALSE)) #center each image (centering across space)
  dat <- dat/sig #standardize by global SD
  dat_ctr <- dat

  return(dat_ctr)

}
