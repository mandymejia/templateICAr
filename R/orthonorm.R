#' Orthonormalizes a square, invertible matrix
#'
#' @param X A square matrix to be orthonormalized.
#'
#' @return X after orthonormalization
#' @export
#'
#' @details Y is orthonormal if $YY'=Y'Y=I$. Orthonormalization of X is given by $X (X'X)^(-.5)$.
#'
orthonorm <- function(X){

  X <- as.matrix(X)

  #check that X is invertible (required for orthonorm(X) %*% t(orthonorm(X)) = I)
  if(!is.finite(determinant(X)$modulus)) stop('X not invertible')

  #compute sqrt of (X'X)
  XtX_sqrt_inv = sqrt_XtX(X, inverse=TRUE)

  #perform orthogonalization
  result = X %*% XtX_sqrt_inv # symmetric orthogonalization
  if(!all.equal(Re(result), result)) stop('Complex-valued result')
  return(result)
}

#' Compute matrix square root of X'X
#'
#' @param X A numerical matrix
#' @param inverse if inverse=TRUE, compute inverse of square root
#'
#' @return A matrix equalling the (inverse) matrix square root of X'X
#' @export
#'
sqrt_XtX <- function(X, inverse=FALSE){

  XtX = t(X) %*% X # X'X = V D^2 V'
  e = eigen(XtX)
  Vmat = e$vectors
  d2 = e$values #diagonal elements of D^2

  if(inverse) {
    if(!is.finite(determinant(X)$modulus)) stop('X not invertible')
    result = Vmat %*% diag(sqrt(1/d2)) %*% t(Vmat)
  } else {
    result = Vmat %*% diag(sqrt(d2)) %*% t(Vmat)
  }

  return(result)
}