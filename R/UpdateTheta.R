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
#'
#' @return An updated list of parameter estimates, theta
#'
NULL

#' @rdname UpdateTheta
#' @export
#' @importFrom stats optimize
#' @importFrom INLA inla.qsample inla.qsolve inla.setOption
#' @import Matrix
UpdateTheta.spatial = function(template_mean, template_var, spde, BOLD, theta, C_diag, common_smoothness=TRUE, verbose=FALSE){

  Q = nrow(BOLD)
  V = ncol(BOLD)

  #initialize new objects
  theta_new = list(A = matrix(NA, Q, Q), nu0_sq = NA, kappa = rep(NA, Q))
  ICmean = matrix(0, Q, V) #subject IC mean
  ICvar = array(0, dim=c(Q,Q,V)) #subject IC var
  A_part1 = A_part2 = matrix(0, Q, Q) #two parts of product for A-hat (construct each looping over voxels)

  A = theta$A
  nu0_sq = theta$nu0_sq
  nu0C_inv = diag(1/(C_diag*nu0_sq))
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

  print('Updating Error Variance nu0_sq')

  Cinv = diag(1/C_diag)
  Cinv_A = Cinv %*% A_hat
  At_Cinv_A = t(A_hat) %*% Cinv %*% A_hat
  nu0sq_part1 = nu0sq_part2 = nu0sq_part3 = 0

  for(v in 1:V){

    y_v = BOLD[,v]

    nu0sq_part1 = nu0sq_part1 + t(y_v) %*% Cinv %*% y_v
    nu0sq_part2 = nu0sq_part2 + t(y_v) %*% Cinv_A %*% miu_s[,v]
    nu0sq_part3 = nu0sq_part3 + sum(diag(At_Cinv_A %*% miu_ssT[,,v]))

  }

  nu0sq_hat = 1/(Q*V)*(nu0sq_part1 - 2*nu0sq_part2 + nu0sq_part3)

  ##########################################
  ### E-STEP for kappa_q: SECOND POSTERIOR MOMENT OF delta_i
  ##########################################

  print('Updating SPDE Parameters kappa_q')

  #SPDE matrices, needed to construct R_q_inv
  F = spde$param.inla$M0
  G = 1/2*(spde$param.inla$M1 + t(spde$param.inla$M1))
  GFinvG = spde$param.inla$M2 #confirmed that this equals G %*% solve(F) %*% G

  #set up Sigma (QVxQV) as a sparse block diagonal matrix
  C1 = 1/(4*pi)
  Sigma_inv_list = vector('list', Q)
  for(q in 1:Q){
    D_q_inv = Diagonal(V, 1/sqrt(template_var[q,])) #sparse diagonal matrix
    kappa_q = theta$kappa[q]
    R_q_inv = C1 * (kappa_q^2 * F + 2 * G + kappa_q^(-2) * GFinvG)
    Sigma_q_inv = D_q_inv %*% R_q_inv %*% D_q_inv
    Sigma_inv_list[[q]] = Sigma_q_inv
  }
  Sigma_inv = bdiag(Sigma_inv_list)

  #set up P as a sparse matrix (see OneNote for illustration of this)
  cols = 1:(Q*V)
  rows_P1 = seq(1, (Q-1)*V+1, by=V)
  offset = rep(0:(V-1), each=Q)
  rows = rep(rows_P1, V) + offset
  P = sparseMatrix(i = rows, j = cols)

  #set up A, C and B = A' Cinv A as a sparse matrix
  ones = Diagonal(V)
  bigA = kronecker(ones, A)
  bigAC = kronecker(ones, t(A) %*% nu0C_inv)
  bigB = kronecker(ones, t(A) %*% nu0C_inv %*% A)

  # COMPUTE TERMS INVOLVING POSTERIOR MEAN AND VARIANCE OF DELTA_iq

  cov_delta_inv = t(P) %*% bigB %*% P + Sigma_inv

  # Two parts of Trace() term of MLE for kappa_q's
  # Part1 = Tr(mu_q' * K_q * mu_q), mu_q = T_q * cov_delta * P' * m --> use inla.qsolve to calculate cov_delta * P' * m
  # Part2 = Tr(K_q * T_q * cov_delta * T_q') --> use inla.qsample to draw from N(0, cov_delta), then estimate necessary elements of cov_delta

  ### Part 1 (solve system of linear equations)

  print('Computing part 1 of trace terms involving big matrix inverse')

  t0 <- Sys.time()
  yvec <- as.vector(BOLD) # [y(1),...,y(V)]
  s0vec <- as.vector(template_mean) # [s0(1),...,s0(V)]
  Pmvec <- t(P) %*% bigAC %*% (yvec - bigA %*% s0vec)
  inla.setOption(smtp="pardiso")
  #1 minute or less using pardiso!!
  cov_delta_Pm <- inla.qsolve(Q = cov_delta_inv, B=matrix(Pmvec, ncol=1), method='solve')
  Trace1_part1 <- rep(NA, Q) #Trace1 = Tr(Dq_inv F Dq_inv E[delta_delta_q])
  Trace2_part1 <- rep(NA, Q) #Trace2 = Tr(Dq_inv G F_inv G Dq_inv E[delta_delta_q])
  #3 seconds!
  for(q in 1:Q){

    # K_q = sparse matrix appearing in trace, D_q^(-1) * GFinvG * D_q^(-1)
    D_q_inv = Diagonal(V, 1/sqrt(template_var[q,]))
    K_q1 <- D_q_inv %*% F %*% D_q_inv
    K_q2 <- D_q_inv %*% GFinvG %*% D_q_inv

    # T_q = matrix that selects the qth block of size V from a matrix with V*Q rows
    e_q <- Matrix(0, nrow=1, ncol=Q, sparse=TRUE)
    e_q[1,q] <- 1
    T_q <- kronecker(e_q, ones)

    #compute mu_q = T_q * cov_delta * P' * m and mu_q' * K_q * mu_q (part 1 of trace)
    mu_q <- T_q %*% cov_delta_Pm
    Trace1_part1[q] <- t(mu_q) %*% K_q1 %*% mu_q
    Trace2_part1[q] <- t(mu_q) %*% K_q2 %*% mu_q
  }
  print(Sys.time() - t0)

  ### Part 2 (Monte Carlo)

  nsamp <- 1000 #number of samples
  print('Estimating part 2 of trace terms involving big matrix inverse')

  print(paste0('Drawing ',nsamp,' Monte Carlo samples'))
  #1-5 min
  print(system.time(musamp <- inla.qsample(nsamp, Q = cov_delta_inv, b=rep(0, V*Q), mu=rep(0, V*Q), num.threads=4))) #1.5 min for 100 samples, 5 min for 1000 samples

  #estimate diagonal blocks of cov_delta
  print('Computing trace of each block diagonal')
  t0 <- Sys.time()
  Trace1_part2 <- rep(NA, Q)
  Trace2_part2 <- rep(NA, Q)
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
    K_q1 <- D_q_inv %*% F %*% D_q_inv
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
    cov_delta_qq_sparse <- t(bigX_left_q) %*% bigX_right_q / (nsamp - 1)

    #compute Trace for this part of the matrix
    Trace1_part2[q] <- sum(K_q1 * cov_delta_qq_sparse) #equivalent to sum(diag(K_q1 %*% cov_delta_qq)) and FAST (note: K_q symmetric)
    Trace2_part2[q] <- sum(K_q2 * cov_delta_qq_sparse) #equivalent to sum(diag(K_q2 %*% cov_delta_qq)) and FAST (note: K_q symmetric)
  }
  print(Sys.time() - t0)

  # NUMERICALLY ESTIMATE MLEs for kappa_q's -- 10 MIN

  print("Performing numerical optimization for kappa_q's")

  if(common_smoothness==FALSE){
    kappa_opt <- rep(NA, Q)
    t0 <- Sys.time()
    for(q in 1:Q){
      if(verbose) print(paste('Optimization ',q,' of ',Q))
      kappa_opt_q <- optimize(Q2_kappa, lower=-10, upper=10, maximum=TRUE,
                              Fmat=F, Gmat=G, GFinvG=GFinvG,
                              bigTrace1=Trace1_part1[q] + Trace1_part2[q],
                              bigTrace2=Trace2_part1[q] + Trace2_part2[q])
      kappa_opt[q] <- exp(kappa_opt_q$maximum)
    }
  } else {
    if(verbose) print('Optimizing kappa')
    kappa_opt <- optimize(Q2_kappa, lower=-10, upper=10, maximum=TRUE,
                          Fmat=F, Gmat=G, GFinvG=GFinvG,
                          bigTrace1=sum(Trace1_part1 + Trace1_part2),
                          bigTrace2=sum(Trace2_part1 + Trace2_part2),
                          Q=Q) #to indicate common smoothness model
    kappa_opt <- rep(exp(kappa_opt$maximum), Q)

  }
  print(Sys.time() - t0)

  # RETURN NEW PARAMETER ESTIMATES

  theta_new$A <- A_hat
  theta_new$nu0_sq <- nu0sq_hat[1]
  theta_new$kappa <- kappa_opt
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

  print('Updating Error Variance nu0_sq')

  #use A-hat or A?

  Cinv = diag(1/C_diag)
  Cinv_A = Cinv %*% A_hat
  At_Cinv_A = t(A_hat) %*% Cinv %*% A_hat
  nu0sq_part1 = nu0sq_part2 = nu0sq_part3 = 0

  for(v in 1:V){

    y_v = BOLD[,v]
    nu0sq_part1 = nu0sq_part1 + t(y_v) %*% Cinv %*% y_v
    nu0sq_part2 = nu0sq_part2 + t(y_v) %*% Cinv_A %*% miu_s[,v]
    nu0sq_part3 = nu0sq_part3 + sum(diag(At_Cinv_A %*% miu_ssT[,,v]))
  }

  nu0sq_hat = 1/(Q*V)*(nu0sq_part1 - 2*nu0sq_part2 + nu0sq_part3)

  # RETURN NEW PARAMETER ESTIMATES

  theta_new$A <- A_hat
  theta_new$nu0_sq <- nu0sq_hat[1]
  return(theta_new)
}


