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
#'
#' @return An updated list of parameter estimates, theta, OR if return_MAP=TRUE, the posterior mean and precision of the latent fields
#'
NULL

#' @rdname UpdateTheta
#' @export
#' @importFrom stats optimize
#' @importFrom INLA inla.qsample inla.qsolve inla.setOption
#' @import Matrix
UpdateTheta.spatial = function(template_mean, template_var, spde, BOLD, theta, C_diag, common_smoothness=TRUE, verbose=FALSE, return_kappa_fun=FALSE, return_MAP=FALSE){

  Q = nrow(template_mean)
  V = ncol(BOLD)

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
  Sigma_inv <- stuff$Sigma_inv
  Fmat <- stuff$Fmat
  Gmat <- stuff$Fmat
  GFinvG <- stuff$GFinvG


  #set up P as a sparse matrix (see OneNote for illustration of this)
  P <- make_Pmat(Q, V)

  #1. Compute m and Omega, then compute mu_sy by solving for x in the system of equations Omega*x = m
  stuff <- compute_mu_s(y_vec, s0_vec, Sigma_inv, theta, P, C_diag)
  mu_s <- stuff$mu
  m_vec <- stuff$m
  Omega <- stuff$Omega

  if(return_MAP){
    return(list(mu_s = mu_s, Omega_s = Omega))
    stop()
  }

  if(verbose) print('Updating A')

  #2. Compute constant term in equation (10): c2 = 1/(1+m'mu_sy)
  c2 <- as.numeric(1/(1 + t(m_vec) %*% mu_s))

  #3. Compute diagonal blocks size Q of y*m'
  #   - Sum along the way and divide by V at the end

  yPm <- matrix(0, nrow=Q, ncol=Q)
  Pm_vec <- P %*% m_vec
  for(v in 1:V){
    inds_v <- (1:Q) + (v-1)*Q
    y_v <- y_vec[inds_v]
    Pm_v <- Pm_vec[inds_v]
    yPm_v <- outer(y_v, Pm_v)
    yPm <- yPm + yPm_v #sum up to later divide by V to average
  }
  A_hat <- c2*yPm/V
  #A_hat = orthonorm(A_hat)
  #fix SD of columns of A on original dimensionality?

  #keep value of nu0sq_hat from PCA-based dimension reduction
  nu0sq_hat <- theta$nu0_sq

  theta_new$A <- A_hat
  theta_new$nu0_sq <- nu0sq_hat[1]

  ##########################################
  ### Update posterior moments of S with A_hat
  ##########################################

  if(verbose) print('Updating posterior moments of S')

  #1. Compute m and Omega, then compute mu_sy by solving for x in the system of equations Omega*x = m
  stuff <- compute_mu_s(y_vec, s0_vec, Sigma_inv, theta_new, P, C_diag)
  mu_s <- stuff$mu
  m_vec <- stuff$m
  Omega <- stuff$Omega

  ##########################################
  ### E-STEP for kappa_q
  ##########################################

  if(verbose) print('Updating SPDE Parameters (kappa)')

  # COMPUTE TERMS INVOLVING POSTERIOR MEAN AND VARIANCE OF s_q

  # Likelihood in terms of kappa_q's.  For r=1,2:
  # PartA = Tr(K_qr * cov_sq) = Tr(K_qr * T_q * cov_s * T_q') --> Use inla.qsample to draw from N(0, cov_s), then estimate necessary elements of cov_sq for q=1,...,Q
  # PartB = mu_sq' * K_qr * mu_sq --> We already have mu_s, just compute mu_sq = T_q * mu_s
  # PartC = s_0q'K_qr*(2*mu_sq - s_0q) --> We already have s_0, just compute s_0q = T_q * s_0 and 2*mu_sq - s_0q
  # PartD = 1/c1 * log(det(kappa^2*Fmat + 2*Gmat + kappa^(-2)*GFinvG))

  nsamp <- 1000 #number of samples

  if(verbose) print('Computing matrix terms in log-likelihood of kappa')

  if(verbose) print(paste0('Drawing ',nsamp,' Monte Carlo samples for covariance estimation:'))
  tmp <- system.time(musamp <- inla.qsample(nsamp, Q = Omega, b=rep(0, V*Q), mu=rep(0, V*Q), num.threads=4))
  if(verbose) print(tmp)
  #30 sec for 1000 samples with V*Q = 4214*3 = 12,642

  t0 <- Sys.time()
  PartA1 <- PartA2 <- rep(NA, Q) #trace part with K_qr, r=1,2
  PartB1 <- PartB2 <- rep(NA, Q) #parts with mu'*K_qr*mu, r=1,2
  PartC1 <- PartC2 <- rep(NA, Q) #parts with s0'*K_qr*(2mu - s0), r=1,2

  #set up estimation of only terms of cov_s necessary for trace computation
  bigX_left <- KhatriRao(diag(1, V), matrix(1, nrow=nsamp, ncol=V)) # -- 13 SEC (do one time instead of KhatriRao(diag(1, V), Xctr_q) for each q)
  mat_nonzero <- as.matrix(1*(GFinvG != 0))
  diag(mat_nonzero) <- 1 #make sure all diagonal elements estimated (required for Trace 1)
  bigX_right <- KhatriRao(mat_nonzero, matrix(1, nrow=nsamp, ncol=V)) # -- 160 SEC (do one time instead of KhatriRao(as.matrix(1*(K_q != 0)), Xctr_q) for each q)
  bigX_right_cols <- which(GFinvG != 0, arr.ind = TRUE)[,2] #column indices of non-zero locations
  inds_left <- which(bigX_left != 0, arr.ind=TRUE)
  inds_right <- which(bigX_right != 0, arr.ind=TRUE)
  for(q in 1:Q){
    if(verbose) print(paste('Block ',q,' of ',Q))
    inds_q <- (1:V) + (q-1)*V
    D_q_inv = Diagonal(V, 1/sqrt(template_var[q,]))
    K_q1 <- D_q_inv %*% Fmat %*% D_q_inv
    K_q2 <- D_q_inv %*% GFinvG %*% D_q_inv

    Xctr_q <- scale(t(musamp[inds_q,]), scale=FALSE) # < 1 sec

    # #compute cov_delta_qq, set unnecessary terms to zero -- 33 SEC
    # cov_delta_qq <- crossprod(Xctr_q)/(nsamp-1)
    # nonzero_q <- which(K_q != 0)
    # cov_delta_qq_zeros <- K_q
    # cov_delta_qq_zeros[nonzero_q] <- cov_delta_qq[nonzero_q]

    # #compute cov_delta_qq (necessary entries only!) -- 13 SEC

    #left-multiplication matrix
    vals_left <- as.vector(Xctr_q)
    bigX_left_q <- sparseMatrix(i = inds_left[,1], j = inds_left[,2], x = vals_left)
    #right-multiplication matrix
    X_repcols <- Xctr_q[,bigX_right_cols]
    vals_right <- as.vector(X_repcols)
    bigX_right_q <- sparseMatrix(i = inds_right[,1], j = inds_right[,2], x = vals_right)
    #multiply together
    cov_s_qq_sparse <- t(bigX_left_q) %*% bigX_right_q / (nsamp - 1)

    # T_q = matrix that selects the qth block of size V from a matrix with V*Q rows
    e_q <- Matrix(0, nrow=1, ncol=Q, sparse=TRUE)
    e_q[1,q] <- 1
    T_q <- kronecker(e_q, Diagonal(V))

    #compute mu_sq = T_q * mu_s and s0_q = T_q * s0
    mu_sq <- T_q %*% mu_s
    s0_q <- T_q %*% s0_vec

    ### PART A

    #compute Trace for this part of the matrix
    PartA1[q] <- sum(K_q1 * cov_s_qq_sparse) #equivalent to sum(diag(K_q1 %*% cov_s_qq)) and FAST (note: K_q symmetric)
    PartA2[q] <- sum(K_q2 * cov_s_qq_sparse) #equivalent to sum(diag(K_q2 %*% cov_s_qq)) and FAST (note: K_q symmetric)

    ### PART B
    PartB1[q] <- t(mu_sq) %*% K_q1 %*% mu_sq
    PartB2[q] <- t(mu_sq) %*% K_q2 %*% mu_sq

    ### PART C
    PartC1[q] <- t(s0_q) %*% K_q1 %*% (2*mu_sq - s0_q)
    PartC2[q] <- t(s0_q) %*% K_q2 %*% (2*mu_sq - s0_q)

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
                              part1 = -1*PartA1[q] - PartB1[q] + PartC1[q], #factor times kappa^2
                              part2 = -1*PartA2[q] - PartB2[q] + PartC2[q]) #factor times kappa^(-2)
      kappa_opt[q] <- (kappa_opt_q$maximum)
    }
  }

  if(common_smoothness==TRUE){
    #if(verbose) print('Optimizing kappa')
    kappa_opt <- optimize(Q2_kappa, lower=0, upper=5, maximum=TRUE,
                          Fmat=Fmat, Gmat=Gmat, GFinvG=GFinvG,
                          part1 = sum(-1*PartA1 - PartB1 + PartC1),
                          part2 = sum(-1*PartA2 - PartB2 + PartC2),
                          Q=Q) #to indicate common smoothness model to the Q2_kappa function
    kappa_opt <- rep(kappa_opt$maximum, Q)

    ### COMPUTE AND RETURN OBJECTIVE FUNCTION IN A NEIGHBORHOOD OF KAPPA-OPT
    if(return_kappa_fun){
      kappa_vals <- exp(seq(floor(log(kappa_opt[1])-1), ceiling(log(kappa_opt[1])+1), 0.05)) #pick range on log scale to avoid values of kappa below zero
      #kappa_vals <- seq(-3, -1, 0.05)
      kappa_opt_fun <- rep(NA, length(kappa_vals))
      part1 = sum(-1*PartA1 - PartB1 + PartC1) #factor times kappa^2
      part2 = sum(-1*PartA2 - PartB2 + PartC2) #factor times kappa^(-2)
      for(k in 1:length(kappa_vals)){ kappa_opt_fun[k] <- Q2_kappa(kappa_vals[k], Fmat=Fmat, Gmat=Gmat, GFinvG=GFinvG, part1, part2, Q=Q) }
    }
  }
  print(Sys.time() - t0)

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
#' @param theta Current parameter estimates
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
compute_mu_s <- function(y_vec, s0_vec, Sigma_inv, theta, P, C_diag){

  Q <- length(C_diag)
  V <- nrow(P)/Q

  A <- theta$A
  nu0_sq <- theta$nu0_sq

  #set up B, C, and products thereof
  ones = Diagonal(V)
  B = kronecker(ones, A)
  nu0C_inv = kronecker(ones, diag(1/(C_diag*nu0_sq)))
  Pt_Bt_nu0C_inv = t(P) %*% t(B) %*% nu0C_inv

  # Compute m and Omega (using previous parameter estimates in theta)
  m_vec <- Pt_Bt_nu0C_inv %*% y_vec + Sigma_inv %*% s0_vec
  Omega <- Sigma_inv + Pt_Bt_nu0C_inv %*% B %*% P

  # Compute mu_s|y by solving for x in the system of equations Omega*x = m
  mu_s <- inla.qsolve(Q = Omega, B=m_vec, method='solve')

  return(list(mu=mu_s, m=m_vec, Omega=Omega))
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

  Q <- nrow(template_var)
  V <- ncol(template_var)

  #SPDE matrices, needed to construct R_q_inv
  Fmat = spde$param.inla$M0
  Gmat = 1/2*(spde$param.inla$M1 + t(spde$param.inla$M1))
  GFinvG = spde$param.inla$M2 #this equals G %*% solve(F) %*% G

  #set up Sigma (QVxQV) as a sparse block diagonal matrix
  Sigma_inv_list = vector('list', Q)
  for(q in 1:Q){
    D_q_inv = Diagonal(V, 1/sqrt(template_var[q,])) #sparse diagonal matrix
    kappa_q = kappa[q]
    R_q_inv = C1 * (kappa_q^2 * Fmat + 2 * Gmat + kappa_q^(-2) * GFinvG)
    Sigma_q_inv = D_q_inv %*% R_q_inv %*% D_q_inv
    Sigma_inv_list[[q]] = Sigma_q_inv
  }
  Sigma_inv = bdiag(Sigma_inv_list)

  return(list(Sigma_inv=Sigma_inv, Fmat=Fmat, Gmat=Gmat, GFinvG=GFinvG))
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







