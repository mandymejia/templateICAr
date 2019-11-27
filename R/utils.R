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


#' Computes part of log-likelihood involving kappa (or kappa_q) for numerical optimization
#'
#' @param logkappa Value of log(kappa) for which to compute log-likelihood
#' @param Fmat A diagonal matrix appearing in SPDE precision
#' @param Gmat A sparse neighborhood matrix appearing in SPDE precision
#' @param GFinvG (Optional) Pre-computed value of \eqn{GFinvG = Gmat * (1/Fmat)* Gmat} (available from SPDE object)
#' @param bigTrace1 Value of first large trace term appearing in optimization function (see Details)
#' @param bigTrace2 Value of second large trace term appearing in optimization function (see Details)
#' @param C1 For the unit variance case, \eqn{\tau^2 = C1/\kappa^2}, where \eqn{C1 = 1/(4\pi)} when \eqn{\alpha=2}, \eqn{\nu=1}, \eqn{d=2}
#' @param Q Equal to the number of ICs for the common smoothness model, or NULL for the IC-specific smoothness model
#'
#' @import Matrix
#' @return Value of log-likelihood at logkappa
#'
#' @details This is the function to be maximized in order to determine the MLE for \eqn{\kappa} or the \eqn{\kappa_q}'s in the M-step of the EM algorithm in spatial
#' template ICA.  In the model where \eqn{\kappa_q} can be different for each IC \eqn{q}, the optimization function factorizes over the \eqn{\kappa_q}'s.  This function computes
#' the value of the part of the optimization function pertaining to one of the \eqn{\kappa_q}'s. The optimization function includes two Trace terms involving
#' very large matrices, which are passed to this function as `bigTrace1` and `bigTrace2`.  `bigTrace1` is equal to
#'
#' \deqn{K1q = Trace(Dq_inv * F * Dq_inv * E[delta_q*t(delta_q)])}
#'
#' and `bigTrace2` is equal to
#'
#' \deqn{K2q = Trace(Dq_inv * G * F_inv * G * Dq_inv * E[delta_q*t(delta_q)]).}
#'
#' For the common smoothness model, a single kappa value is assumed over all of the IC's. The two Trace terms are then equal to \eqn{\sum_q(K1q)} and \eqn{\sum_q(K2q)}.
#'
Q2_kappa <- function(logkappa, Fmat, Gmat, GFinvG=NULL, bigTrace1, bigTrace2, C1 = 1/(4*pi), Q=NULL){

  if(is.null(GFinvG)) GFinvG <- Gmat %*% solve(Fmat) %*% Gmat

  kappa <- exp(logkappa)

  #log determinant part
  part1_mat <- kappa^2 * Fmat + 2 * Gmat + kappa^(-2) * GFinvG
  logdet <- as.numeric(determinant(part1_mat)) #on log scale
  if(logdet[2] == -1) warning('negative determinant of precision matrix')
  part1 <- logdet[1]
  if(!is.null(Q)) part1 <- Q*part1

  #first trace term in Q2
  part2a <- C1*(kappa^2) * bigTrace1

  #second trace term in Q2
  part2b <- C1*(kappa^(-2)) * bigTrace2

  result <- part1 - part2a - part2b
  return(result)

  # part1 decreases with kappa
  # part2=part2a+part2b also decreases with kappa
  # maximum depends on which part dominates
  # return(list(part1, - part2a - part2b, part1 - part2a - part2b))

}

