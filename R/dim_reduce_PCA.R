#' PCA-based Dimension Reduction and Prewhitening
#'
#' Performs dimension reduction and prewhitening based on probabilistic PCA using 
#'  SVD. If dimensionality is not specified, it is estimated using the method 
#'  described in Minka (2008).
#'
#' @param X \eqn{V \times T} fMRI timeseries data matrix, centered by columns.
#' @param Q Number of latent dimensions to estimate. If \code{NULL} (default), 
#'  estimated using PESEL (Sobczyka et al. 2020).
#' @param Q_max Maximal number of principal components for automatic 
#'  dimensionality selection with PESEL. Default: \code{100}.
#' @return A list containing the dimension-reduced data (\code{data_reduced}, a 
#'  \eqn{V \times Q} matrix), prewhitening/dimension reduction matrix (\code{H}, 
#'  a \eqn{QxT} matrix) and its (pseudo-)inverse (\code{Hinv}, a \eqn{TxQ} 
#'  matrix), the noise variance (\code{sigma_sq}), the correlation matrix of the
#'  dimension-reduced data (\code{C_diag}, a \eqn{QxQ} matrix), and the 
#'  dimensionality (\code{Q}).
#' 
#' @export
#' @examples
#' nT <- 30
#' nV <- 400
#' nQ <- 7
#' X <- matrix(rnorm(nV*nQ), nrow=nV) %*% diag(seq(nQ, 1)) %*% matrix(rnorm(nQ*nT), nrow=nQ) 
#' dim_reduce(X, Q=nQ)
#' 
dim_reduce <- function(X, Q=NULL, Q_max=100){

  svd_XtX <- PCA(X, Q=Q, Q_max=Q_max)
  U <- svd_XtX$u
  D1 <- svd_XtX$d[1:Q]
  D2 <- svd_XtX$d[(Q+1):length(svd_XtX$d)]

  #residual variance
  sigma_sq <- mean(D2)

  #prewhitening matrix
  H <- diag(1/sqrt(D1 - sigma_sq)) %*% t(U)
  H_inv <- U %*% diag(sqrt(D1 - sigma_sq))

  #for residual variance after prewhitening
  C_diag <- rowSums(H^2) # diag(H %*% t(H))

  #prewhitened data (transposed)
  X_new <- X %*% t(H)

  list(
    data_reduced=X_new, 
    H=H, 
    H_inv=H_inv, 
    sigma_sq=sigma_sq, 
    C_diag=C_diag, 
    Q=Q
  )
}

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
    X <- colCenter(X)
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