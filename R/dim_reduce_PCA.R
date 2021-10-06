#' PCA-based Dimension Reduction and Prewhitening for ICA
#'
#' Performs dimension reduction and prewhitening based on probablistic PCA using 
#'  SVD. If dimensionality is not specified, it is estimated using the method 
#'  described in Minka (2008).
#'
#' @param X \eqn{VxT} fMRI timeseries data matrix, centered by columns and rows 
#'  (columns are actually all that matter, but MATLAB implementation of Minka 
#'  method also assumes rows have been centered (implicit in use of cov function))
#' @param Q Number of latent dimensions to estimate. If not specified, 
#'  estimated using PESEL (Sobczyka et al. 2020).
#' @param Q_max Maximal number of principal components for automatic 
#'  dimensionality selection with PESEL. Default: \code{100}
#'
#' @importFrom pesel pesel
#' 
#' @return A list containing the dimension-reduced data (data_reduced, a 
#'  \eqn{VxQ} matrix), prewhitening/dimension reduction matrix (H, a \eqn{QxT} 
#'  matrix) and its (pseudo-)inverse (Hinv, a \eqn{TxQ} matrix), the noise variance 
#'  (sigma_sq), the correlation matrix of the dimension-reduced data 
#'  (C_diag, a \eqn{QxQ} matrix), and the dimensionality (\eqn{Q})
#' 
#' @export
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
#' @param X \eqn{VxT} fMRI timeseries data matrix, centered by columns and rows 
#'  (columns are actually all that matter, but MATLAB implementation of Minka 
#'  method also assumes rows have been centered (implicit in use of cov function))
#' @param Q Number of latent dimensions to estimate. If not specified, 
#'  estimated using PESEL (Sobczyka et al. 2020).
#' @param Q_max Maximal number of principal components for automatic 
#'  dimensionality selection with PESEL. Default: \code{100}
#' 
#' @importFrom pesel pesel
#' 
#' @return The SVD decomposition
#' 
PCA <- function(X, Xt=FALSE, Q=NULL, Q_max=100) {
  nvox <- nrow(X) #number of brain locations
  ntime <- ncol(X) #number of fMRI volumes (reduce this)
  if(ntime > nvox) warning('More time points than voxels. Are you sure?')

  # check that X has been centered both ways
  TOL <- 1e-8
  if (max(abs(colMeans(X))) > TOL) stop('Columns of X must be centered')
  if (max(abs(rowMeans(X))) > TOL) warning('Rows of X should be centered')

  #determine PCA dimensionality
  if(is.null(Q)){
    Q <- pesel(X, npc.max=Q_max, method='homogenous')$nPCs
  }

  #perform dimension reduction
  svd(crossprod(X) / nvox, nu=Q, nv=0)
}