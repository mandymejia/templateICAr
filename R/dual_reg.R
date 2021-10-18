#' Dual Regression
#'
#' @param dat Subject-level fMRI data (\eqn{VxT})
#' @param GICA Group-level independent components (\eqn{VxQ})
#' @param scale A logical value indicating whether the fMRI timeseries should
#'  be scaled by the image standard deviation.
#'
#' @importFrom matrixStats colVars
#' 
#' @return A list containing the subject-level independent components \strong{S} (\eqn{VxQ}), 
#'  and subject-level mixing matrix \strong{A} (\eqn{TxQ}).
#' 
#' @export
#'
dual_reg <- function(dat, GICA, scale=FALSE){

  ntime <- ncol(dat) #length of timeseries
  nvox <- nrow(dat) #number of data locations
  if(ntime > nvox) warning('More time points than voxels. Are you sure?')
  if(nvox != nrow(GICA)) stop('The number of voxels in dat and GICA must match')

  Q <- ncol(GICA) #number of ICs
  if(Q > nvox) warning('More ICs than voxels. Are you sure?')
  if(Q > ntime) warning('More ICs than time points. Are you sure?')

  # center timeseries data across space and time (and standardize scale if scale=TRUE)
  # transpose it
  dat_ctr <- t(scale_BOLD(dat, scale=scale))

  #center each group IC over voxels
  GICA - rep(colMeans(GICA), rep.int(nvox, Q))

	#estimate A (IC timeseries)
	A <- (dat_ctr %*% GICA) %*% chol2inv(chol(crossprod(GICA)))
	#estimate S (IC maps)
	S <- solve(a=crossprod(A), b=crossprod(A, dat_ctr))

	#fix scale of spatial maps (sd=1)
	#sd_S <- sqrt(rowVars(S))
	#A <- A %*% diag(sd_S)
	#S <- diag(1/sd_S) %*% S

	#return result
	list(S = S, A = A)
}
