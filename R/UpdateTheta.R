#' @name UpdateTheta
#' @rdname UpdateTheta
#'
#' @title Parameter Estimates in EM Algorithm for Template ICA Model
#'
#' @param template_mean mean maps for each IC in template
#' @param template_var between-subject variance maps for each IC in template
#' @param BOLD dimension-reduced fMRI data
#' @param spde NULL for spatial independence model, otherwise SPDE object representing spatial prior on deviations.
#' @param theta A list of current parameter estimates (mixing matrix A, noise variance nu0_sq and (for spatial model) SPDE parameters kappa)
#' @param C_diag (Qx1) diagonal elements of matrix proportional to residual variance.
#' @param common_smoothness If TRUE, use the common smoothness version of the spatial template ICA model, which assumes that all IC's have the same smoothness parameter, \eqn{\kappa}
#' @param verbose If TRUE, print progress updates for slow steps.
#' @param return_kappa_fun If TRUE, return the log likelihood as a function of kappa in a neighborhood of the MLE (common smoothness model only)
#' @param return_MAP If TRUE, returns the posterior mean and precision of the latent fields instead of the parameter estimates
#' @param dim_reduce_flag If FALSE, data is in the original resolution (no dimension reduction).
#'
#' @return An updated list of parameter estimates, theta, OR if return_MAP=TRUE, the posterior mean and precision of the latent fields
#'
NULL

#' @rdname UpdateTheta
#' @export
#' @importFrom stats optimize
#' @importFrom INLA inla.qsample inla.qsolve inla.setOption
#' @import Matrix
UpdateTheta.spatial = function(template_mean, template_var, spde, BOLD, theta, C_diag, Hinv, common_smoothness=TRUE, verbose=FALSE, return_kappa_fun=FALSE, return_MAP=FALSE, dim_reduce_flag){

  Q = nrow(template_mean)
  V = ncol(BOLD)
  ntime = nrow(BOLD)

  #initialize new parameter values
  theta_new = list(A = NULL, nu0_sq = NULL, kappa = NULL)

  #1. Compute mu_sy by solving for x in the system of equations Omega*x = m
  #   Compute m and Omega (using previous parameter estimates)
  #   Omega = Sigma_inv + P'B'(nu0C_inv)BP
  #   m = P'B'(nu0C_inv)y + Sigma_inv*s0

  #2. Compute constant term in equation (10): 1/(1+m'mu_sy)

  #3. Compute diagonal blocks size Q of y*m'
  #   - Sum along the way and divide by V at the end

  #4. Compute A-hat:
  #   - Multiply block avg by constant term from (2)
  #   - Orthogonalize
  #   - Make each column have var=1 on original dimensionality

  ### Ingredients for Omega and m

  # s0_vec = "s_0" (long vector organized by ICs)
  # y (long vector organized by locations)

  # Sigma_inv = D_inv*R_inv*D_inv, where D_inv contains template sd (big diagonal matrix) & R_inv is SPDE spatial corr matrix (big block diagonal matrix)

  # P (permutation matrix)

  # nu0C_inv = "(nu_0^2 C)^{-1}" (big diagonal matrix)
  # B = "I_v \otimes A" (big block diagonal matrix)

  if(verbose) print('Computing Posterior Moments of S')

  s0_vec = as.vector(t(template_mean))
  y_vec = as.vector(BOLD)

  #Compute SPDE matrices (F, G, GFinvG) and Sigma_inv (QVxQV), a sparse block diagonal matrix
  stuff <- compute_Sigma_inv(spde, kappa=theta$kappa, template_var, C1=1/(4*pi))
  R_inv <- stuff$R_inv
  D <- stuff$D
  Fmat <- stuff$Fmat
  Gmat <- stuff$Gmat
  GFinvG <- stuff$GFinvG

  #set up P as a sparse matrix (see OneNote for illustration of this)
  P <- make_Pmat(Q, V)

  #1. Compute mu_s
  stuff <- compute_mu_s(y_vec, s0_vec, R_inv, D, theta, P, C_diag)
  mu_s <- stuff$mu
  m_vec <- stuff$m
  Omega <- stuff$Omega
  Omega_inv_m <- stuff$Omega_inv_m

  if(return_MAP){
    return(list(mu_s = mu_s, Omega_s = Omega))
    stop()
  }

  if(verbose) print('Updating A')

  #2. Compute constant term in equation (10): c2 = 1/(1+m'mu_sy)
  #c2 <- as.numeric(1/(1 + t(m_vec) %*% mu_s))

  #3. Compute diagonal blocks size Q of y*m'P'
  #   - Sum along the way and divide by V at the end

  #Compute first matrix appearing in A-hat
  Pmu_vec <- P %*% mu_s
  yPmu <- matrix(0, nrow=ntime, ncol=Q)
  for(v in 1:V){
    inds_y <- (1:ntime) + (v-1)*ntime
    inds_Pmu <- (1:Q) + (v-1)*Q
    y_v <- y_vec[inds_y]
    Pmu_v <- Pmu_vec[inds_Pmu]
    yPmu_v <- outer(y_v, Pmu_v)
    yPmu <- yPmu + yPmu_v #sum up to later divide by V to average
  }

  #Compute second matrix appearing in A-hat
  Omega_PP <- P %*% Omega %*% t(P)
  D_PP <- P %*% D %*% t(P)

  nsamp <- 1000 #number of samples
  if(verbose) print(paste0('Drawing ',nsamp,' Monte Carlo samples for estimation of covariance terms'))
  tmp <- system.time(musamp <- inla.qsample(nsamp, Q = Omega_PP, b=rep(0, V*Q), mu=rep(0, V*Q), num.threads=4))
  if(verbose) print(tmp) #50 seconds for Q=3 and V=4214

  T_mat <- matrix(0, Q, Q)
  for(v in 1:V){
    inds <- (1:Q) + (v-1)*Q

    #T2 matrix
    Pmu_v <- Pmu_vec[inds]
    T2_vv <- outer(Pmu_v, Pmu_v)

    #T1 matrix
    Xctr_v <- scale(t(musamp[inds,]), scale=FALSE) # < 1 sec
    Omega_PP_inv_vv <- crossprod(Xctr_v)/(nsamp-1)
    T1_vv <- D_PP[inds,inds] %*% Omega_PP_inv_vv %*% D_PP[inds,inds]

    #T(v,v)
    T_vv <- T1_vv + T2_vv
    T_mat <- T_mat + T_vv
  }
  #T_mat_inv <- solve(T_mat)

  #A-hat = yPmu * T_mat^(-1)
  A_hat_t <- solve(t(T_mat), t(yPmu)) #A_hat <- yPmu %*% T_mat_inv

  # yPm <- matrix(0, nrow=ntime, ncol=Q)
  # Pm_vec <- P %*% m_vec
  # for(v in 1:V){
  #   inds_y <- (1:ntime) + (v-1)*ntime
  #   inds_Pm <- (1:Q) + (v-1)*Q
  #   y_v <- y_vec[inds_y]
  #   Pm_v <- Pm_vec[inds_Pm]
  #   yPm_v <- outer(y_v, Pm_v)
  #   yPm <- yPm + yPm_v #sum up to later divide by V to average
  # }
  # A_hat <- c2*yPm
  print(colVars(as.matrix(A_hat)))
  if(dim_reduce_flag==TRUE) A_hat = orthonorm(A_hat)
  # #fix SD of columns of A on original dimensionality?
  # sd_A <- sqrt(colVars(Hinv %*% A_hat)) #get scale of A (after reverse-prewhitening)
  # sd_A <- sqrt(colVars(A_hat)) #get scale of A
  # A_hat <- A_hat %*% diag(1/sd_A) #standardize scale of A


  #keep value of nu0sq_hat from PCA-based dimension reduction
  nu0sq_hat <- theta$nu0_sq

  theta_new$A <- A_hat
  theta_new$nu0_sq <- nu0sq_hat[1]

  ##########################################
  ### Update posterior moments of S with A_hat
  ##########################################

  if(verbose) print('Updating posterior moments of S')

  #1. Compute m and Omega, then compute mu_sy by solving for x in the system of equations Omega*x = m
  stuff <- compute_mu_s(y_vec, s0_vec, R_inv, D, theta_new, P, C_diag)
  mu_s <- stuff$mu
  m_vec <- stuff$m
  Omega <- stuff$Omega
  Omega_inv_m <- stuff$Omega_inv_m

  ##########################################
  ### E-STEP for kappa_q
  ##########################################

  if(verbose) print('Updating SPDE Parameters (kappa)')

  # Likelihood in terms of kappa_q's.  With K1=F and K2=GFinvG:
  # PartA = Tr(K1 * Omega_inv_qq) + Tr(K1 * M_q) + u_q'*K1*v_q
  # PartB = Tr(K2 * Omega_inv_qq) + Tr(K2 * M_q) + u_q'*K2*v_q
  # PartC = 1/c1 * log(det(kappa^2*Fmat + 2*Gmat + kappa^(-2)*GFinvG))

  # Tr(K1 * Omega_inv_qq) --> Use inla.qsample to draw from N(0, Omega_inv_qq), then estimate necessary elements (non-zeros in K1) of Omega_inv_qq for q=1,...,Q
  # Tr(K1 * M_q) --> M = outer(Omega_inv*m,Omega_inv*m), where Omega_inv*m is known. Just calculate necessary elements (non-zeros in K1) of M_q for q=1,...,Q

  if(verbose) print('Computing trace terms in log-likelihood of kappa')

  # SETUP FOR FIRST TERM IN PARTS A & B (Trace with Omega-inv-hat)

  nsamp <- 1000 #number of samples
  if(verbose) print(paste0('Drawing ',nsamp,' Monte Carlo samples for covariance estimation:'))
  tmp <- system.time(musamp <- inla.qsample(nsamp, Q = Omega, b=rep(0, V*Q), mu=rep(0, V*Q), num.threads=4)) #sample from N(0, Omega^(-1)) to estimate Omega^(-1) (diagonal blocks only)
  if(verbose) print(tmp)
  #30-60 sec for 1000 samples with V*Q = 4214*3 = 12,642

  #set up estimation of only terms necessary for trace computation
  nonzero_GFinvG <- as.matrix(1*(GFinvG != 0))
  diag(nonzero_GFinvG) <- 1 #make sure all diagonal elements estimated (required for Trace 1)
  nonzero_cols <- which(GFinvG != 0, arr.ind = TRUE)[,2] #column indices of non-zero locations

  X <- matrix(1, nrow=nsamp, ncol=V) #placeholder for X in Omega^(-1) = XX'
  bigX_left <- KhatriRao(diag(1, V), X) # -- 13 SEC (do one time instead of KhatriRao(diag(1, V), Xctr_q) for each q)
  bigX_right <- KhatriRao(nonzero_GFinvG, X) # -- 160 SEC (do one time instead of KhatriRao(as.matrix(1*(K_q != 0)), Xctr_q) for each q)
  inds_left_X <- which(bigX_left != 0, arr.ind=TRUE)
  inds_right_X <- which(bigX_right != 0, arr.ind=TRUE)

  # SETUP FOR SECOND TERM IN PARTS A & B (Trace with M-hat)

  x <- matrix(1, nrow=1, ncol=V) #placeholder for x in M = xx'
  bigx_left <- KhatriRao(diag(1, V), x) # -- 13 SEC (do one time instead of KhatriRao(diag(1, V), Xctr_q) for each q)
  bigx_right <- KhatriRao(nonzero_GFinvG, x) # -- 160 SEC (do one time instead of KhatriRao(as.matrix(1*(K_q != 0)), Xctr_q) for each q)
  inds_left_x <- which(bigx_left != 0, arr.ind=TRUE)
  inds_right_x <- which(bigx_right != 0, arr.ind=TRUE)

  # SETUP FOR THIRD TERM IN PARTS A & B (u'*K*v)

  u <- inla.qsolve(Q = D, B=matrix(s0_vec, ncol=1), method='solve') #Dinv_s0 = D^(-1)*s0
  v <- 2*Omega_inv_m - u

  t0 <- Sys.time()
  PartA1 <- PartB1 <- rep(NA, Q) #Tr(Kr * Omega_inv_qq), r=1,2
  PartA2 <- PartB2 <- rep(NA, Q) #Tr(Kr * M_q), r=1,2
  PartA3 <- PartB3 <- rep(NA, Q) #u_q'*Kr*v_q, r=1,2

  K1 <- Fmat
  K2 <- GFinvG

  for(q in 1:Q){
    if(verbose) print(paste('Block ',q,' of ',Q))
    inds_q <- (1:V) + (q-1)*V

    # COMPUTE OMEGA_INV[q,q] (NECESSARY ENTRIES ONLY!) -- 13 SEC

    Xctr_q <- scale(t(musamp[inds_q,]), scale=FALSE) # < 1 sec

    #left-multiplication matrix
    vals_left_X <- as.vector(Xctr_q)
    bigX_left_q <- sparseMatrix(i = inds_left_X[,1], j = inds_left_X[,2], x = vals_left_X)
    #right-multiplication matrix
    X_repcols <- Xctr_q[,nonzero_cols]
    vals_right_X <- as.vector(X_repcols)
    bigX_right_q <- sparseMatrix(i = inds_right_X[,1], j = inds_right_X[,2], x = vals_right_X)
    #multiply together
    Omega_inv_qq <- t(bigX_left_q) %*% bigX_right_q / (nsamp - 1)

    # COMPUTE M[q,q]

    # T_q = matrix that selects the qth block of size V from a matrix with V*Q rows
    e_q <- Matrix(0, nrow=1, ncol=Q, sparse=TRUE)
    e_q[1,q] <- 1
    T_q <- kronecker(e_q, Diagonal(V))

    Omega_inv_m_q <- T_q %*% Omega_inv_m

    #compute M = Omega_inv_m * Omega_inv_m' (necessary entries only)

    #left-multiplication matrix
    vals_left_x <- as.vector(Omega_inv_m_q)
    bigx_left_q <- sparseMatrix(i = inds_left_x[,1], j = inds_left_x[,2], x = vals_left_x)
    #right-multiplication matrix
    x_repcols <- t(Omega_inv_m_q)[,nonzero_cols]
    vals_right_x <- as.vector(x_repcols)
    bigx_right_q <- sparseMatrix(i = inds_right_x[,1], j = inds_right_x[,2], x = vals_right_x)
    #multiply together
    M_hat_qq <- t(bigx_left_q) %*% bigx_right_q

    ### PARTS A1 & B1

    #compute Trace for this part of the matrix
    PartA1[q] <- sum(K1 * Omega_inv_qq) #equivalent to sum(diag(K_q1 %*% cov_s_qq)) and FAST (note: K_q symmetric)
    PartB1[q] <- sum(K2 * Omega_inv_qq) #equivalent to sum(diag(K_q2 %*% cov_s_qq)) and FAST (note: K_q symmetric)

    ### PARTS A2 & B2

    PartA2[q] <- sum(K1 * M_hat_qq)
    PartB2[q] <- sum(K2 * M_hat_qq)

    ### PARTS A3 & B3

    u_q <- T_q %*% u
    v_q <- T_q %*% v
    PartA3[q] <- t(u_q) %*% K1 %*% v_q
    PartB3[q] <- t(u_q) %*% K2 %*% v_q

  }
  if(verbose) print(Sys.time() - t0)

  # NUMERICALLY ESTIMATE MLEs for kappa_q's -- 10 MIN

  print("Performing numerical optimization for kappa")

  t0 <- Sys.time()
  if(common_smoothness==FALSE){
    kappa_opt <- rep(NA, Q)
    for(q in 1:Q){
      if(verbose) print(paste('Optimization ',q,' of ',Q))
      kappa_opt_q <- optimize(Q2_kappa, lower=-10, upper=10, maximum=TRUE,
                              Fmat=Fmat, Gmat=Gmat, GFinvG=GFinvG,
                              part1 = -1*PartA1[q] - PartA2[q] + PartA3[q], #factor times kappa^2
                              part2 = -1*PartB1[q] - PartB2[q] + PartB3[q]) #factor times kappa^(-2)
      kappa_opt[q] <- (kappa_opt_q$maximum)
    }
  }

  if(common_smoothness==TRUE){
    #if(verbose) print('Optimizing kappa')
    kappa_opt <- optimize(Q2_kappa, lower=0, upper=5, maximum=TRUE,
                          Fmat=Fmat, Gmat=Gmat, GFinvG=GFinvG,
                          part1 = sum(-1*PartA1 - PartA2 + PartA3),
                          part2 = sum(-1*PartB1 - PartB2 + PartB3),
                          Q=Q) #to indicate common smoothness model to the Q2_kappa function
    kappa_opt <- rep(kappa_opt$maximum, Q)

    ### COMPUTE AND RETURN OBJECTIVE FUNCTION IN A NEIGHBORHOOD OF KAPPA-OPT
    if(return_kappa_fun){
      kappa_vals <- exp(seq(floor(log(kappa_opt[1])-1), ceiling(log(kappa_opt[1])+2), 0.05)) #pick range on log scale to avoid values of kappa below zero
      #kappa_vals <- seq(-3, -1, 0.05)
      kappa_opt_fun <- rep(NA, length(kappa_vals))
      part1 = sum(-1*PartA1 - PartA2 + PartA3) #factor times kappa^2
      part2 = sum(-1*PartB1 - PartB2 + PartB3) #factor times kappa^(-2)
      for(k in 1:length(kappa_vals)){ kappa_opt_fun[k] <- Q2_kappa(kappa_vals[k], Fmat=Fmat, Gmat=Gmat, GFinvG=GFinvG, part1, part2, Q=Q) }
    }
  }
  print(Sys.time() - t0)

  # plot(log(kappa_vals), kappa_opt_fun, type='l', ylim=c(-1e6, 6e5),
  #      main = expression(paste("MLE of ",kappa[q])), xlab = expression(log(kappa[q])), ylab = expression(f[q]))
  # abline(v=log(kappa_opt), lty=2)

  # RETURN NEW PARAMETER ESTIMATES

  theta_new$kappa <- kappa_opt
  if(return_kappa_fun) theta_new$kappa_fun <- cbind(kappa_vals, kappa_opt_fun)

  return(theta_new)

}

#' @rdname UpdateTheta
#' @export
UpdateTheta.independent = function(template_mean, template_var, BOLD, theta, C_diag){

  Q = nrow(BOLD)
  V = ncol(BOLD)

  #initialize new objects
  theta_new = list(A = matrix(NA, Q, Q), nu0_sq = NA)
  A_part1 = A_part2 = matrix(0, Q, Q) #two parts of product for A-hat (construct each looping over voxels)

  A = theta$A
  nu0_sq = theta$nu0_sq
  nu0C_inv = diag(1/(C_diag*nu0_sq)) #Sigma0_inv in matlab code
  At_nu0Cinv = t(A) %*% nu0C_inv
  At_nu0Cinv_A = At_nu0Cinv %*% A

  print('Updating A')

  #store posterior moments for M-step of nu0_sq
  miu_s = matrix(NA, nrow=Q, ncol=V)
  miu_ssT = array(NA, dim=c(Q, Q, V))

  for(v in 1:V){

    y_v = BOLD[,v]
    s0_v = template_mean[,v]

    ##########################################
    ### E-STEP FOR A AND nu0^2: POSTERIOR MOMENTS OF s_i(v)
    ##########################################

    E_v_inv = diag(1/template_var[,v])
    Sigma_s_v = solve(E_v_inv + At_nu0Cinv_A)
    miu_s_v = Sigma_s_v	%*% (At_nu0Cinv %*% y_v + E_v_inv %*% s0_v) #Qx1
    miu_ssT_v = (miu_s_v %*% t(miu_s_v)) + Sigma_s_v #QxQ
    miu_s[,v] = miu_s_v #save for M-step of nu0_sq
    miu_ssT[,,v] = miu_ssT_v #save for M-step of nu0_sq

    ##########################################
    ### M-STEP FOR A: CONSTRUCT PARAMETER ESTIMATES
    ##########################################

    A_part1 = A_part1 + y_v %*% t(miu_s_v) #QxQ
    A_part2 = A_part2 + miu_ssT_v #QxQ

  }

  A_hat = orthonorm(A_part1 %*% solve(A_part2))

  ##########################################
  ### M-STEP FOR nu0^2: CONSTRUCT PARAMETER ESTIMATES
  ##########################################

  # print('Updating Error Variance nu0_sq')
  #
  # #use A-hat or A?
  #
  # Cinv = diag(1/C_diag)
  # Cinv_A = Cinv %*% A_hat
  # At_Cinv_A = t(A_hat) %*% Cinv %*% A_hat
  # nu0sq_part1 = nu0sq_part2 = nu0sq_part3 = 0
  #
  # for(v in 1:V){
  #
  #   y_v = BOLD[,v]
  #   nu0sq_part1 = nu0sq_part1 + t(y_v) %*% Cinv %*% y_v
  #   nu0sq_part2 = nu0sq_part2 + t(y_v) %*% Cinv_A %*% miu_s[,v]
  #   nu0sq_part3 = nu0sq_part3 + sum(diag(At_Cinv_A %*% miu_ssT[,,v]))
  # }
  #
  # nu0sq_hat = 1/(Q*V)*(nu0sq_part1 - 2*nu0sq_part2 + nu0sq_part3)

  nu0sq_hat <- theta$nu0_sq


  # RETURN NEW PARAMETER ESTIMATES

  theta_new$A <- A_hat
  theta_new$nu0_sq <- nu0sq_hat[1]
  return(theta_new)
}



#' Computes posterior mean and precision matrix of s
#'
#' @param y_vec Vectorized, dimension-reduced fMRI data, grouped by locations. A vector of length QV.
#' @param s0_vec Vectorized template mean, grouped by ICs.A vector of length QV.
#' @param R_inv Estimate of inverse spatial correlation matrix (sparse)
#' @param D Diagonal matrix containing square root of template variance values along the diagonal (sparse)
#' @param theta List of current parameter estimates
#' @param P Permutation matrix for regrouping by locations (instead of by ICs.)
#' @param C_diag Diagonals of residual covariance of the first level model. A vector of length Q.
#'
#' @return A list containing the posterior mean \eqn{\mu} (mu) and precision \eqn{\Omega} (Omega) of s=(s1,...,sQ), along with the supporting vector m, where \eqn{\mu = \Omega^{-1}m}.
#'
#' @details
#'
#' @import Matrix
#' @importFrom INLA inla.qsolve
#'
compute_mu_s <- function(y_vec, s0_vec, R_inv, D, theta, P, C_diag){

  ntime <- length(C_diag)
  Q <- ncol(theta$A)
  V <- nrow(P)/Q

  A <- theta$A
  nu0_sq <- theta$nu0_sq

  #set up B, C, and products thereof
  ones = Diagonal(V)
  B = kronecker(ones, A)
  nu0C_inv = kronecker(ones, diag(1/(C_diag*nu0_sq)))
  Pt_Bt_nu0C_inv = t(P) %*% t(B) %*% nu0C_inv

  #compute D^(-1)*s0
  Dinv_s0 <- inla.qsolve(Q = D, B=matrix(s0_vec, ncol=1), method='solve')

  #compute m (using current parameter estimates in theta)
  m1_vec <- D %*% Pt_Bt_nu0C_inv %*% y_vec
  m2_vec <- R_inv %*% Dinv_s0
  m_vec <- m1_vec + m2_vec

  # Compute Omega (using current parameter estimates in theta)
  Omega <- R_inv + D %*% Pt_Bt_nu0C_inv %*% B %*% P %*% D

  # Compute mu_s|y by solving for x in the system of equations Omega*x = m
  Omega_inv_m <- inla.qsolve(Q = Omega, B=m_vec, method='solve')
  mu_s <- D %*% Omega_inv_m

  return(list(mu=mu_s, m=m_vec, Omega=Omega, Omega_inv_m=Omega_inv_m))
}


#' Creates precision matrix for S and SPDE matrices (F, G and GFinvG)
#'
#' @param spde Object of class SPDE representing spatial process
#' @param kappa Current estimates of SPDE parameter kappa for each latent field
#' @param template_var QxV matrix of template variances
#' @param C1 Constant, equal to 1/(4*pi) for a 2-dimensional field with alpha=2
#'
#' @return A list containing Sigma inverse and SPDE matrices
#'
#' @details
#'
#' @import Matrix
#'
compute_Sigma_inv <- function(spde, kappa, template_var, C1=1/(4*pi)){

  Q <- length(kappa)

  #SPDE matrices, needed to construct R_q_inv
  Fmat = spde$param.inla$M0
  Gmat = 1/2*(spde$param.inla$M1 + t(spde$param.inla$M1))
  GFinvG = spde$param.inla$M2 #this equals G %*% solve(F) %*% G

  #set up Sigma (QVxQV) as a sparse block diagonal matrix
  R_inv_list = vector('list', Q)
  D_list = vector('list', Q)
  for(q in 1:Q){
    kappa_q = kappa[q]
    R_q_inv = C1 * (kappa_q^2 * Fmat + 2 * Gmat + kappa_q^(-2) * GFinvG)
    D_q = Diagonal(V, sqrt(template_var[q,])) #sparse diagonal matrix
    R_inv_list[[q]] = R_q_inv
    D_list[[q]] = D_q
  }
  R_inv = bdiag(R_inv_list)
  D = bdiag(D_list)

  return(list(R_inv=R_inv, D=D, Fmat=Fmat, Gmat=Gmat, GFinvG=GFinvG))
}

#' Creates permutation matrix P to regroup elements of s by locations instead of by ICs
#'
#' @param Q The number of template ICs
#' @param V The number of spatial locations
#'
#' @return P Permutation matrix size QVxQV
#'
#' @details If s=(s1,...,sQ) is grouped by ICs 1,...Q, then Ps=(s(1),...,s(V)) is grouped by locations 1,...,V
#'
#' @importFrom Matrix sparseMatrix
#'
make_Pmat <- function(Q, V){
  cols = 1:(Q*V)
  rows_P1 = seq(1, (Q-1)*V+1, by=V)
  offset = rep(0:(V-1), each=Q)
  rows = rep(rows_P1, V) + offset
  P = t(sparseMatrix(i = rows, j = cols))
}







