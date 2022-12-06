
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
#   Winv_D_delta <- INLA::inla.qsolve(Q = W, B=matrix(D_delta, ncol=1), method='solve')
#   exp_part1 <- 1/sigma_sq * sum(delta^2)
#   exp_part2 <- 1/(sigma_sq^2) * t(D_delta) %*% Winv_D_delta
#   exp_part <- -1+as.numeric(exp_part1) + as.numeric(exp_part2)
#
#   return(-1*(det_part + exp_part)) #return negative log-likelihood for minimization
#
#
# }
