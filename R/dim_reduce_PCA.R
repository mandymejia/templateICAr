#' PCA
#' 
#' Efficient PCA for a tall matrix (many more rows than columns). Uses the SVD
#'  of the covariance matrix.
#' 
#' @param X \eqn{V \times T} fMRI timeseries data matrix, centered by columns.
#' @param center Center the columns of \code{X}? Default: \code{TRUE}. Set to
#'  \code{FALSE} if already centered. 
#' @param Q Number of latent dimensions to estimate. If \code{NULL} (default), 
#'  estimated using PESEL (Sobczyka et al. 2020).
#' @param Q_max Maximal number of principal components for automatic 
#'  dimensionality selection with PESEL. Default: \code{100}.
#' @param nV Number of principal directions to obtain. Default: \code{0}. Can
#'  also be \code{"Q"} to set equal to the value of \code{Q}. Note that setting
#'  this value less than \code{Q} does not speed up computation time, but does
#'  save on memory. Note that the directions will be with respect to \code{X},
#'  not its covariance matrix.
#' 
#' @importFrom fMRItools colCenter
#' @importFrom pesel pesel
#' @export
#' 
#' @return The SVD decomposition
#' 
PCA <- function(X, center=TRUE, Q=NULL, Q_max=100, nV=0) {
  nvox <- nrow(X) #number of brain locations
  ntime <- ncol(X) #number of fMRI volumes (reduce this)
  if(ntime > nvox) warning('More time points than voxels. Are you sure?')

  if (center) {
    X <- fMRItools::colCenter(X)
  } else {
    TOL <- 1e-8
    if (max(abs(colMeans(X))) > TOL) stop('Columns of X must be centered')
  }
  # Removed check for row centering
  
  # Determine PCA dimensionality
  if(is.null(Q)){
    Q <- suppressWarnings(pesel(X, npc.max=Q_max, method='homogenous'))$nPCs
  }
  if (nV == "Q") { nV <- Q }
  if (nV > Q) { warning("nV > Q, so setting nV to Q."); nV <- Q }

  # Perform dimension reduction
  out <- svd(crossprod(X) / (nvox-1), nu=Q, nv=0)

  # Compute directions. (out$v would have the directions for XtX, not X.)
  if (nV > 0) {
    out$v <- X %*% out$u[,seq(nV)] %*% diag(1/out$d[seq(nV)])
  }

  out
}