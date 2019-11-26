#' Orthonormalizes a square, invertible matrix
#'
#' @param X A square matrix to be orthonormalized.
#'
#' @return X after orthonormalization
#' @export
#'
#' @details Y is orthonormal if $YY'=Y'Y=I$. Orthonormalization of X is given by $X (X'X)^(-.5)$.
#'
orthonorm = function(X){

  #check that X is invertible (required for orthonorm(X) %*% t(orthonorm(X)) = I)
  if(!is.finite(determinant(X)$modulus)) stop('X not invertible')

  #compute sqrt of (X'X)
  XtX_sqrt_inv = sqrt_XtX(X, inverse=TRUE)

  #perform orthogonalization
  result = X %*% XtX_sqrt_inv # symmetric orthogonalization
  if(!all.equal(Re(result), result)) stop('Complex-valued result')
  return(result)

}

####################################################################################
####################################################################################
### sqrt_XtX() -- computes matrix square root of X'X
####################################################################################
####################################################################################

#' Compute matrix square root of X'X
#'
#' @param X A numerical matrix
#' @param inverse if inverse=TRUE, compute inverse of square root
#'
#' @return A matrix equalling the (inverse) matrix square root of X'X
#' @export
#'
sqrt_XtX = function(X, inverse=FALSE){

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


#' Computes part of log-likelihood involving kappa_q for numerical optimization
#'
#' @param logkappa_q Value of log(kappa) for which to compute log-likelihood
#' @param Fmat A diagonal matrix appearing in SPDE precision
#' @param Gmat A sparse neighborhood matrix appearing in SPDE precision
#' @param GFinvG `GFinvG = Gmat %*% (1/Fmat) %*% Gmat` (available directly through SPDE object)
#' @param bigTrace1 Equal to `Trace(Dq_inv * F * Dq_inv * E[delta_q*t(delta_q)])`
#' @param bigTrace2 Equal to `Trace(Dq_inv * G * F_inv * G * Dq_inv * E[delta_q*t(delta_q)])`
#' @param C1 For the unit variance case, $tau^2 = C1/kappa^2$.  $C1 = 1/(4*pi)$ when $alpha=2$, $nu=1$, $d=2$
#'
#' @import Matrix
#' @return Value of log-likelihood at logkappa_q
#'
Q2_kappa_q <- function(logkappa_q, Fmat, Gmat, GFinvG=NULL, bigTrace1, bigTrace2, C1 = 1/(4*pi)){

  if(is.null(GFinvG)) GFinvG <- Gmat %*% solve(Fmat) %*% Gmat

  kappa_q <- exp(logkappa_q)

  #log determinant part
  part1_mat <- kappa_q^2 * Fmat + 2 * Gmat + kappa_q^(-2) * GFinvG
  logdet <- as.numeric(determinant(part1_mat)) #on log scale
  if(logdet[2] == -1) warning('negative determinant of precision matrix')
  part1 <- logdet[1]

  #first trace term in Q2
  part2a <- C1*(kappa_q^2) * bigTrace1

  #second trace term in Q2
  part2b <- C1*(kappa_q^(-2)) * bigTrace2

  result <- part1 - part2a - part2b
  return(result)

  # part1 decreases with kappa
  # part2=part2a+part2b also decreases with kappa
  # maximum depends on which part dominates
  # return(list(part1, - part2a - part2b, part1 - part2a - part2b))

}

