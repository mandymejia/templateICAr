#' PCA-based Dimension Reduction and Prewhitening
#'
#' Performs dimension reduction and prewhitening based on probabilistic PCA using 
#'  SVD. If dimensionality is not specified, it is estimated using the method 
#'  described in Minka (2008).
#'
#' @param X A numeric matrix, with each column being a centered timeseries. 
#'  For fMRI data, \code{X} should be \code{T} timepoints by \code{V} brain 
#'  locations.
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
#' @importFrom fMRItools PCA
#' @export
#' 
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