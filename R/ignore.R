#' #' Orthonormalizes a square, invertible matrix
#' #'
#' #' @param X A square matrix to be orthonormalized.
#' #'
#' #' @return X after orthonormalization
#' #' @export
#' #'
#' #' @details Y is orthonormal if $YY'=Y'Y=I$. Orthonormalization of X is given by $X (X'X)^(-.5)$.
#' #'
#' orthonorm = function(X){
#'
#'   X <- as.matrix(X)
#'
#'   #check that X is invertible (required for orthonorm(X) %*% t(orthonorm(X)) = I)
#'   if(!is.finite(determinant(X)$modulus)) stop('X not invertible')
#'
#'   #compute sqrt of (X'X)
#'   XtX_sqrt_inv = sqrt_XtX(X, inverse=TRUE)
#'
#'   #perform orthogonalization
#'   result = X %*% XtX_sqrt_inv # symmetric orthogonalization
#'   if(!all.equal(Re(result), result)) stop('Complex-valued result')
#'   return(result)
#'
#' }

####################################################################################
####################################################################################
### sqrt_XtX() -- computes matrix square root of X'X
####################################################################################
####################################################################################

#' #' Compute matrix square root of X'X
#' #'
#' #' @param X A numerical matrix
#' #' @param inverse if inverse=TRUE, compute inverse of square root
#' #'
#' #' @return A matrix equalling the (inverse) matrix square root of X'X
#' #' @export
#' #'
#' sqrt_XtX = function(X, inverse=FALSE){
#'
#'   XtX = t(X) %*% X # X'X = V D^2 V'
#'   e = eigen(XtX)
#'   Vmat = e$vectors
#'   d2 = e$values #diagonal elements of D^2
#'
#'   if(inverse) {
#'     if(!is.finite(determinant(X)$modulus)) stop('X not invertible')
#'     result = Vmat %*% diag(sqrt(1/d2)) %*% t(Vmat)
#'   } else {
#'     result = Vmat %*% diag(sqrt(d2)) %*% t(Vmat)
#'   }
#'
#'   return(result)
#'
#' }
#'

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
