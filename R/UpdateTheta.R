#' @name UpdateTheta
#' @rdname UpdateTheta
#'
#' @title Parameter Estimates in EM Algorithm for Template ICA Model
#'
#' @param template_mean (QxV matrix) mean maps for each IC in template
#' @param template_var (QxV matrix) between-subject variance maps for each IC in template
#' @param BOLD dimension-reduced fMRI data
#' @param spde NULL for spatial independence model, otherwise SPDE object representing spatial prior on deviations.
#' @param theta A list of current parameter estimates (mixing matrix A, noise variance nu0_sq and (for spatial model) SPDE parameters kappa)
#' @param C_diag (Qx1) diagonal elements of matrix proportional to residual variance.
#' @param s0_vec
#' @param D
#' @param Dinv_s0
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
UpdateTheta.spatial = function(template_mean, template_var, mesh, BOLD, theta, C_diag, s0_vec, D, Dinv_s0, common_smoothness=TRUE, verbose=FALSE, return_kappa_fun=FALSE, return_MAP=FALSE, dim_reduce_flag){

  Q = nrow(template_mean)
  V = ncol(BOLD)
  ntime = nrow(BOLD)
  spde = mesh$spde

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

  if(verbose) cat('Computing Posterior Moments of S \n')

  y_vec = as.vector(BOLD)

  if(verbose) cat('...posterior precision \n') # less than 1 sec
  #Compute SPDE matrices (F, G, GFinvG) and Sigma_inv (QVxQV), a sparse block diagonal matrix
  stuff <- compute_Sigma_inv(mesh, kappa=theta$kappa, template_var, C1=1/(4*pi))
  R_inv <- stuff$R_inv
  Fmat <- stuff$Fmat
  Gmat <- stuff$Gmat
  GFinvG <- stuff$GFinvG

  #set up P as a sparse matrix (see OneNote for illustration of this)
  P <- make_Pmat(Q, V)

  #1. Compute mu_s
  if(verbose) cat('...posterior mean \n') #20 seconds with pardiso! (1 hour without)
  stuff <- compute_mu_s(y_vec, D, Dinv_s0, R_inv, theta, P, C_diag)
  mu_s <- stuff$mu
  m_vec <- stuff$m
  Omega <- stuff$Omega
  Omega_inv_m <- stuff$Omega_inv_m

  if(verbose) print(summary(as.vector(mu_s)))

  if(return_MAP){
    return(list(mu_s = mu_s, Omega_s = Omega))
    stop()
  }

  if(verbose) cat('Updating A \n')
  t00 <- Sys.time()

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

  if(verbose) cat(paste0('Drawing ',nsamp,' Monte Carlo samples for estimation of covariance estimation \n'))

  nsamp <- 500 #number of samples
  musamp <- inla.qsample(nsamp, Q = Omega_PP, b=rep(0, V*Q), mu=rep(0, V*Q), num.threads=8)
  #1000 samples --> 50 seconds for Q=3 and V=4214
  #1000 samples --> > 1 hour for Q=16 and V=5500
  #100 samples --> 45 minutes for Q=16 and V=5500

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

  A_hat <- t(solve(t(T_mat), t(yPmu))) #A_hat <- yPmu %*% T_mat_inv

  if(dim_reduce_flag==TRUE) A_hat = orthonorm(A_hat)

  #keep value of nu0sq_hat from PCA-based dimension reduction
  nu0sq_hat <- theta$nu0_sq

  theta_new$A <- A_hat
  theta_new$nu0_sq <- nu0sq_hat[1]

  ##########################################
  ### Update posterior moments of S with A_hat
  ##########################################

  # if(verbose) cat('Updating posterior mean of S \n')
  #
  # #1. Compute m and Omega, then compute mu_sy by solving for x in the system of equations Omega*x = m
  # system.time(stuff <- compute_mu_s(y_vec, D, Dinv_s0, R_inv, theta, P, C_diag))
  # mu_s <- stuff$mu
  # m_vec <- stuff$m
  # Omega <- stuff$Omega
  # Omega_inv_m <- stuff$Omega_inv_m
  #
  # print(summary(as.vector(mu_s)))

  ##########################################
  ### E-STEP for kappa
  ##########################################

  if(verbose) cat('Updating kappa  \n')

  # Likelihood in terms of kappa_q's.
  # LL2_part1 = log(det(R_q_inv))
  # LL2_part2 = Tr(R_q_inv * Omega_inv_q) + Tr(R_q_inv * W_hat_q)
  # LL2_part3 = u_q' * R_q_inv * v_hat_q
  #             u_q = Dinv * s0
  #             v_q = 2 Omega_inv * m - Dinv * s0

  # Tr(R_q_inv * Omega_inv_q) --> Use inla.qsample to draw from N(0, Omega_inv_q), then estimate necessary elements (non-zeros in R_q_inv) of Omega_inv_q for q=1,...,Q
  # Tr(R_q_inv * W_hat_q) --> W_hat_q = outer(Omega_inv*m,Omega_inv*m)_qq, where Omega_inv*m is known. Just calculate necessary elements (non-zeros in R_q_inv) of W_hat_q for q=1,...,Q

  if(verbose) cat('..Computing trace terms in log-likelihood of kappa \n')

  ### SET UP FOR PART 2

  #1. Generate Monte Carlo samples for estimation of Omega_hat^{-1}_qq

  if(verbose) cat(paste0('....drawing ',nsamp,' Monte Carlo samples for covariance estimation \n'))

  nsamp <- 500 #number of samples
  musamp <- inla.qsample(nsamp, Q = Omega, b=rep(0, V*Q), mu=rep(0, V*Q), num.threads=8) #sample from N(0, Omega^(-1)) to estimate Omega^(-1) (diagonal blocks only)
  #9 sec for 500 samples with V*Q = 2530*3 = 12,642

  #2. Determine non-zero terms of R_q_inv

  if(verbose) cat('....setting up helper objects for trace computation \n') #15 sec (Q=16, V=5500)

  #set up estimation of only terms necessary for trace computation
  nonzero_Rinv <- as.matrix(1*(R_inv[1:V,1:V] != 0))
  #diag(nonzero_Rinv) <- 1 #make sure all diagonal elements estimated (required for Trace 1)
  nonzero_cols <- which(R_inv[1:V,1:V] != 0, arr.ind = TRUE)[,2] #column indices of non-zero locations

  #3. Set up matrices needed for computation of only necessary elements of Omega_hat^{-1}

  X <- matrix(1, nrow=nsamp, ncol=V) #placeholder for X in Omega^(-1) = XX'
  bigX_left <- KhatriRao(diag(1, V), X) # -- 13 SEC (do one time instead of KhatriRao(diag(1, V), Xctr_q) for each q)
  bigX_right <- KhatriRao(nonzero_Rinv, X) # -- 160 SEC (do one time instead of KhatriRao(as.matrix(1*(K_q != 0)), Xctr_q) for each q)
  inds_left_X <- which(bigX_left != 0, arr.ind=TRUE)
  inds_right_X <- which(bigX_right != 0, arr.ind=TRUE)

  #4. Set up vectors needed for computation of only necessary elements of W_hat

  x <- matrix(1, nrow=1, ncol=V) #placeholder for x in M = xx'
  bigx_left <- KhatriRao(diag(1, V), x) # -- 13 SEC (do one time instead of KhatriRao(diag(1, V), Xctr_q) for each q)
  bigx_right <- KhatriRao(nonzero_Rinv, x) # -- 160 SEC (do one time instead of KhatriRao(as.matrix(1*(K_q != 0)), Xctr_q) for each q)
  inds_left_x <- which(bigx_left != 0, arr.ind=TRUE)
  inds_right_x <- which(bigx_right != 0, arr.ind=TRUE)

  if(verbose) cat('....computing necessary elements of RHS matrices in trace terms \n')
  #15 sec for Q=16, V=5500

  if(!common_smoothness) OplusW <- vector('list', length=Q) #Omega_inv_qq + W_hat_qq

  for(q in 1:Q){
    #if(verbose) cat(paste('......block ',q,' of ',Q,' \n'))
    inds_q <- (1:V) + (q-1)*V

    # COMPUTE OMEGA_INV[q,q] (NECESSARY ENTRIES ONLY)

    Xctr_q <- scale(t(musamp[inds_q,]), scale=FALSE)

    #left-multiplication matrix
    vals_left_X <- as.vector(Xctr_q)
    bigX_left_q <- sparseMatrix(i = inds_left_X[,1], j = inds_left_X[,2], x = vals_left_X)
    #right-multiplication matrix
    X_repcols <- Xctr_q[,nonzero_cols]
    vals_right_X <- as.vector(X_repcols)
    bigX_right_q <- sparseMatrix(i = inds_right_X[,1], j = inds_right_X[,2], x = vals_right_X)
    #multiply together
    Omega_inv_qq <- t(bigX_left_q) %*% bigX_right_q / (nsamp - 1)

    # COMPUTE W_hat[q,q] = Omega_inv_m[q] * Omega_inv_m[q' (NECESSARY ENTRIES ONLY)

    # T_q = matrix that selects the qth block of size V from a matrix with V*Q rows
    e_q <- Matrix(0, nrow=1, ncol=Q, sparse=TRUE)
    e_q[1,q] <- 1
    T_q <- kronecker(e_q, Diagonal(V))
    Omega_inv_m_q <- T_q %*% Omega_inv_m

    #left-multiplication matrix
    vals_left_x <- as.vector(Omega_inv_m_q)
    bigx_left_q <- sparseMatrix(i = inds_left_x[,1], j = inds_left_x[,2], x = vals_left_x)
    #right-multiplication matrix
    x_repcols <- t(Omega_inv_m_q)[,nonzero_cols]
    vals_right_x <- as.vector(x_repcols)
    bigx_right_q <- sparseMatrix(i = inds_right_x[,1], j = inds_right_x[,2], x = vals_right_x)
    #multiply together
    W_hat_qq <- t(bigx_left_q) %*% bigx_right_q

    # COMBINE OMEGA_INV[q,q] AND W_hat[q,q] --> two trace terms involving R_q_inv can be combined

    OplusW_qq <- Omega_inv_qq + W_hat_qq
    if(common_smoothness){
      if(q==1) OplusW <- OplusW_qq else OplusW <- OplusW + OplusW_qq
    } else {
      OplusW[[q]] <- OplusW_qq
    }
  }


  ### SET UP FOR PART 3

  u <- Dinv_s0
  v <- 2*Omega_inv_m - u

  # NUMERICALLY SOLVE FOR MLE OF KAPPA

  cat("..Performing numerical optimization \n")
  #~8 seconds for V=5500, Q=3 (comon smoothness)

  if(common_smoothness){
    kappa_opt <- optimize(LL2_kappa, lower=0, upper=5, maximum=TRUE,
             spde=spde, OplusW=OplusW, u=u, v=v, Q=Q)  #Q=Q to indicate common smoothness model to the LL2_kappa function
    kappa_opt <- rep(kappa_opt$maximum, Q)
  } else {
    kappa_opt <- rep(NA, Q)
    for(q in 1:Q){
      if(verbose) cat(paste('Optimization ',q,' of ',Q,' \n'))
      inds_q <- (1:V) + (q-1)*V
      kappa_opt_q <- optimize(LL2_kappa, lower=0, upper=5, maximum=TRUE,
             spde=spde, OplusW=OplusW[[q]], u=u[inds_q], v=v[inds_q], Q=NULL)
      kappa_opt[q] <- (kappa_opt_q$maximum)
    }
  }

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

  cat('Updating A \n')

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

  # cat('Updating Error Variance nu0_sq \n')
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
#' @import Matrix
#' @importFrom INLA inla.qsolve
#'
compute_mu_s <- function(y_vec, D, Dinv_s0, R_inv, theta, P, C_diag){

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

  #compute m (using current parameter estimates in theta)
  m1_vec <- D %*% Pt_Bt_nu0C_inv %*% y_vec
  m2_vec <- R_inv %*% Dinv_s0
  m_vec <- m1_vec + m2_vec

  # Compute Omega (using current parameter estimates in theta)
  Omega <- R_inv + D %*% Pt_Bt_nu0C_inv %*% B %*% P %*% D

  # Compute mu_s|y by solving for x in the system of equations Omega*x = m
  Omega_inv_m <- inla.qsolve(Q = Omega, B=m_vec, method='solve') # <-- slowest part (up to 1 hour for Q=16, V=5500), but with inla.setOption(smtp='pardiso') goes down to 20 seconds!
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
#' @import Matrix
#'
compute_Sigma_inv <- function(mesh, kappa, template_var, C1=1/(4*pi)){

  Q <- length(kappa)
  V <- ncol(template_var)

  #SPDE matrices, needed to construct R_q_inv
  spde = mesh$spde
  Fmat = spde$param.inla$M0
  Gmat = 1/2*(spde$param.inla$M1 + Matrix::t(spde$param.inla$M1))
  GFinvG = spde$param.inla$M2 #this equals G %*% solve(F) %*% G

  if(var(kappa)==0) onekappa <- TRUE
  if(length(kappa)==1) onekappa <- TRUE
  if(onekappa) kappa <- kappa[1]

  #get inmesh and notinmesh indices
  Amat = mesh$A #n_loc x n_mesh
  N = ncol(mesh$A) #number of mesh locations
  inmesh = which(colSums(Amat) > 0)
  notinmesh = setdiff(1:N, inmesh)


  #set up R^{-1} (QVxQV) as a sparse block diagonal matrix
  if(onekappa){ #just compute Rinv once
    Qmat = C1*(kappa^2 * Fmat + 2 * Gmat + kappa^(-2) * GFinvG)
    Q11 = Qmat[inmesh,inmesh] # = Amat %*% Qmat %*% t(Amat)
    Q12 = Qmat[inmesh, notinmesh]
    Q21 = Qmat[notinmesh, inmesh]
    Q22 = Qmat[notinmesh,notinmesh]
    Q22_inv <- solve(Q22)
    R_q_inv = Q11 - (Q12 %*% Q22_inv %*% Q21)
    R_inv_list <- rep(list(R_q_inv), Q)
  } else { # compute Rinv block-wise
    R_inv_list = vector('list', Q)
    for(q in 1:Q){
      kappa_q = kappa[q]
      R_q_inv = C1 * (kappa_q^2 * Fmat + 2 * Gmat + kappa_q^(-2) * GFinvG)
      R_inv_list[[q]] = R_q_inv
    }
  }
  R_inv <- bdiag(R_inv_list)

  return(list(R_inv=R_inv, Fmat=Fmat, Gmat=Gmat, GFinvG=GFinvG))
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







