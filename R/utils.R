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
#' @param kappa Value of kappa for which to compute log-likelihood
#' @param spde Object of class SPDE representing spatial process
#' @param OplusW Sparse matrix containing estimated values of RHS of trace in part 2 of log-likelihood. In common smoothness model, represents the sum over q=1,...,Q.
#' @param u Vector needed for part 3 of log-likelihood
#' @param v Vector needed for part 3 of log-likelihood
#' @param C1 For the unit variance case, \eqn{\tau^2 = C1/\kappa^2}, where \eqn{C1 = 1/(4\pi)} when \eqn{\alpha=2}, \eqn{\nu=1}, \eqn{d=2}
#' @param Q Equal to the number of ICs for the common smoothness model, or NULL for the IC-specific smoothness model
#'
#' @import Matrix
#' @return Value of log-likelihood at logkappa
#'
#' @details This is the function to be maximized in order to determine the MLE for \eqn{\kappa} or the \eqn{\kappa_q}'s in the M-step of the EM algorithm in spatial
#' template ICA.  In the model where \eqn{\kappa_q} can be different for each IC \eqn{q}, the optimization function factorizes over the \eqn{\kappa_q}'s.  This function computes
#' the value of the part of the optimization function pertaining to one of the \eqn{\kappa_q}'s. (UPDATE THIS:) The optimization function includes two Trace terms involving
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
LL2_kappa <- function(kappa, spde, OplusW, u, v, C1 = 1/(4*pi), Q=NULL){

  #COMPUTE R_q_inv FOR CURRENT VALUE OF kappa

  Fmat = spde$param.inla$M0
  Gmat = 1/2*(spde$param.inla$M1 + t(spde$param.inla$M1))
  GFinvG = spde$param.inla$M2 #this equals G %*% solve(F) %*% G
  Qmat = C1*(kappa^2 * Fmat + 2 * Gmat + kappa^(-2) * GFinvG)
  Q11 = Qmat[inmesh,inmesh] # = Amat %*% Qmat %*% t(Amat)
  Q12 = Qmat[inmesh, notinmesh]
  Q21 = Qmat[notinmesh, inmesh]
  Q22 = Qmat[notinmesh,notinmesh]
  Q22_inv <- solve(Q22)
  R_q_inv = Q11 - (Q12 %*% Q22_inv %*% Q21)


  #CALCULATE LOG LIKELIHOOD FOR KAPPA IN 3 PARTS


  ### PART 1: log(det(R_q_inv))

  chol_Rinv <- chol(R_q_inv) #R_inv = chol_Rinv' * chol_Rinv
  LL2_part1 <- 2*sum(log(diag(chol_Rinv))) #log determinant of R_q_inv
  if(!is.null(Q)) LL2_part1 <- Q*LL2_part1 #for common smoothness model

  ### PART 2: Trace terms

  #OplusW = Omega_inv_qq + W_hat_qq (sum over q for common smoothness case)
  LL2_part2 <- sum(R_q_inv * OplusW) #equivalent to sum(diag(R_inv_qq %*% OplusW[q])) and FAST (note: R_inv_qq symmetric)

  ### PART 3: u_q' * R_q_inv * v_hat_q (sum over q for common smoothness case)

  if(!is.null(Q)) {
    R_inv <- bdiag(rep(list(R_q_inv), Q))
    LL2_part3 <- t(u) %*% R_inv %*% v
  } else {
    LL2_part3 <- t(u) %*% R_q_inv %*% v
  }

  result <- LL2_part1 - LL2_part2 + LL2_part3
  return(as.numeric(result))

}


# ####################################################################
# # loglik_kappa() - computes log likelihood of kappa for known deviations
# ####################################################################
#
# #ASSUME DELTA AND D_DIAG ALREADY PROJECTED TO MESH LOCATIONS
# #par = c(log-kappa, log-sigma_sq) for type='estimated' or just kappa for type='true'
# #delta is the known or estimated subject effect (deviation)
# #D_diag is the diagonal elements of the matrix D_q, which contains the sqrt template variance values
# loglik_kappa_est <- function(par, delta, D_diag, mesh, C1 = 1/(4*pi)){
#
#   require(Matrix)
#
#   kappa <- exp(par[1])
#   sigma_sq <- exp(par[2])
#
#   #SPDE matrices, needed for oracle MLE of kappa
#   spde = mesh$spde
#   n_mesh = spde$n.spde
#   Fmat = spde$param.inla$M0
#   Gmat = 1/2*(spde$param.inla$M1 + t(spde$param.inla$M1))
#   GFinvG = spde$param.inla$M2 #this equals G %*% solve(F) %*% G
#   Rinv = C1*(kappa^2 * Fmat + 2 * Gmat + kappa^(-2) * GFinvG)
#   cholR = chol(Rinv) #Rmat = cholR'cholR, log(det(Rmat)) = 2*sum(log(diag(cholR)))
#   det_Rinv <- 2*sum(log(diag(cholR)))
#
#   D_diag[D_diag == 0] <- 0.0001
#   D = Diagonal(length(D_diag), as.vector(D_diag))
#   delta = as.vector(delta)
#
#   W = Rinv + 1/sigma_sq* t(mesh$A) %*% (D^2) %*% mesh$A
#   cholW = chol(W)
#   det_W <- 2*sum(log(diag(cholW)))
#
#   #compute determinant part of log-likelihood
#   det_part <- det_Rinv - det_W - n_mesh*log(sigma_sq)
#   if(det_Rinv == Inf | det_W == Inf) {
#     stop('negative determinant of precision matrix, returning NA')
#     return(NA)
#   }
#
#   #compute exponential part of log-likelihood
#   D_delta <- t(mesh$A) %*% D %*% delta
#   Winv_D_delta <- inla.qsolve(Q = W, B=matrix(D_delta, ncol=1), method='solve')
#   exp_part1 <- 1/sigma_sq * sum(delta^2)
#   exp_part2 <- 1/(sigma_sq^2) * t(D_delta) %*% Winv_D_delta
#   exp_part <- -1+as.numeric(exp_part1) + as.numeric(exp_part2)
#
#   return(-1*(det_part + exp_part)) #return negative log-likelihood for minimization
#
#
# }

