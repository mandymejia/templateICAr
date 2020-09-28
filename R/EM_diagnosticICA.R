#' @name EM_diagnosticICA
#' @rdname EM_diagnosticICA
#'
#' @title EM Algorithms for Diagnostic ICA Models
#'
#' @param template_mean (A list of G matrices, each VxL) template mean estimates for each group 1 to G, where L is the number of ICs, V is the number of data or mesh locations.
#' @param template_var (A list of G matrices, each VxL) template variance estimates for each group 1 to G
#' @param meshes Either NULL (assume spatial independence) or a list of objects of type \code{templateICA_mesh} created by \code{make_mesh()} (each list element corresponds to one brain structure)
#' @param BOLD  (VxL matrix) dimension-reduced fMRI data
#' @param theta0 (list) initial guess at parameter values: A (LxL mixing matrix) and (for spatial model only) kappa (SPDE smoothness parameter)
#' @param C_diag (Lx1) diagonal elements of residual covariance after dimension reduction
#' @param maxiter maximum number of EM iterations
#' @param epsilon smallest proportion change between iterations (e.g. .001)
#' @param verbose If TRUE, display progress of algorithm
#' @param ignore_determinant For spatial model only. If TRUE, ignore the normalizing constant in p(y|z) when computing posterior probabilities of z.
#'
#' @return  A list containing:
#' group_probs (posterior probabilities of group membership),
#' subICmean (estimates of subject-level ICs),
#' subICvar (variance of subject-level ICs) for non-spatial model only OR
#' subICcov (covariance of subject-level ICs, only computed for neighboring locations) and subICprec (precision of subject-level ICs) for spatial model only,
#' theta_MLE (list of final parameter estimates),
#' theta_path (path of parameter updates in SQUAREM, spatial model only),
#' success_flag (flag indicating convergence (\code{TRUE}) or not (\code{FALSE})),
#' error (change in parameter estimates at final iteration),
#' numiter (number of parameter update steps, approximately 3x the number of accelerated SQUAREM iterations),
#' squarem (output of SQUAREM function call),
#' and the template mean and variance used in model estimation
#'
#' @details \code{EM_dignosticICA.spatial} implements the expectation-maximization (EM) algorithm described in Mejia (2020+) for estimating
#' the probabilities of group membership, subject-level ICs and unknown parameters in the diagnostic ICA model with spatial priors on subject effects.
#'
#'
NULL


#' @rdname EM_diagnosticICA
#' @importFrom INLA inla.spde2.matern inla.qsolve
#' @importFrom Matrix Diagonal
#' @import SQUAREM
#'
EM_diagnosticICA.spatial = function(template_mean, template_var, meshes, BOLD, theta0, C_diag, maxiter=100, epsilon=0.001, verbose=TRUE, ignore_determinant=TRUE){

  nvox <- nrow(BOLD) #number of brain locations
  if(ncol(BOLD) > nvox) warning('More time points than data locations. Are you sure the data is oriented properly?')
  L <- ncol(template_mean[[1]]) #number of ICs
  G <- length(template_mean)
  if(L > nvox) stop('Cannot estimate more ICs than brain locations.')
  if(L != ncol(BOLD)) stop('Dimension-reduced BOLD data should have number of columns equal to number of ICs in templates')

  if(any(sapply(meshes, class) != 'templateICA_mesh')) stop('Each element of meshes argument should be of class templateICA_mesh. See help(make_mesh).')

  iter = 1
  theta = theta0
  success = 1
  for(g in 1:G) {
    num_smallvar <- sum(template_var[[g]] < 1e-6)
    if(num_smallvar>0){
      if(verbose) cat(paste0('Setting ',num_smallvar,' (',round(100*num_smallvar/length(template_var[[g]]),1),'%) very small variance values in group ',g,' template to ',1e-6,'.\n'))
      template_var[[g]][template_var[[g]] < 1e-6] = 1e-6 #to prevent problems when inverting covariance
    }
  }

  #pre-compute s0, D and D^{-1}*s0
  s0_vec_list <- vector('list', length=G)
  D_list <- vector('list', length=G)
  Dinv_s0_list <- vector('list', length=G)
  for(g in 1:G){
    s0_vec_list[[g]] = as.vector(template_mean[[g]]) #grouped by IC
    D_vec <- as.vector(sqrt(template_var[[g]])) #grouped by IC
    D_list[[g]] = Diagonal(nvox*L, D_vec)
    Dinv_s0_list[[g]] <- inla.qsolve(Q = D_list[[g]], B=matrix(s0_vec_list[[g]], ncol=1), method='solve')
  }

  #theta0 <- theta1 #last tested value of kappa0
  theta0$LL <- c(0,0) #log likelihood
  theta0_vec <- unlist(theta0[1:2]) #everything but LL
  names(theta0_vec)[1] <- 0 #store LL value in names of theta0_vec (required for squarem)


  t00000 <- Sys.time()
  result_squarem <- squarem(par=theta0_vec, fixptfn = UpdateThetaSQUAREM_diagnosticICA, objfn=LL_SQUAREM_diagnosticICA, control=list(trace=verbose, intermed=TRUE, tol=epsilon, maxiter=maxiter), template_mean, template_var, meshes, BOLD, C_diag, s0_vec_list, D_list, Dinv_s0_list, verbose=TRUE)
  if(verbose) print(Sys.time() - t00000)

  # err = 1000 #large initial value for difference between iterations
  # while(err > epsilon){
  #
  #   if(verbose) cat(paste0(' ~~~~~~~~~~~~~~~~~~~~~ ITERATION ', iter, ' ~~~~~~~~~~~~~~~~~~~~~ \n'))
  #
  #   t00 <- Sys.time()
  #   theta_new <- UpdateTheta_diagnosticICA.spatial(template_mean,
  #                                                  template_var,
  #                                                  meshes,
  #                                                  BOLD,
  #                                                  theta,
  #                                                  C_diag,
  #                                                  s0_vec_list,
  #                                                  D_list,
  #                                                  Dinv_s0_list,
  #                                                  verbose=verbose,
  #                                                  return_MAP=FALSE,
  #                                                  update=c('all'),
  #                                                  ignore_determinant=TRUE)
  #   if(verbose) print(Sys.time() - t00)
  #
  #   ### Compute change in parameters
  #
  #   A_old = theta$A
  #   A_new = theta_new$A
  #   kappa_old = theta$kappa
  #   kappa_new = theta_new$kappa
  #   #2-norm = largest eigenvalue = sqrt of largest eigenvalue of AA'
  #   err1 = norm(as.vector(A_new - A_old), type="2")/norm(as.vector(A_old), type="2")
  #   err2 = norm(kappa_new[1] - kappa_old[1], type="2")/norm(kappa_old[1], type="2")
  #   change1 = format(err1, digits=3, nsmall=3)
  #   change2 = format(err2, digits=3, nsmall=3)
  #   if(verbose) cat(paste0('Iteration ',iter, ': Difference is ',change1,' for A and ',change2,' for kappa\n'))
  #   err = max(err1, err2)
  #
  #   ### Move to next iteration
  #   theta <- theta_new
  #   iter = iter + 1
  #   if(iter > maxiter){
  #     success = 0;
  #     warning(paste0('Failed to converge within ', maxiter,' iterations'))
  #     break() #exit loop
  #   }
  # }

  path_A <- result_squarem$p.inter[,1:(L^2)]
  path_kappa <- result_squarem$p.inter[,(L^2+1)]
  path_LL <- result_squarem$p.inter[,ncol(result_squarem$p.inter)]
  theta_path <- list(A=path_A, kappa=path_kappa, LL=path_LL)

  theta_MLE <- theta0
  theta_MLE$A <- matrix(result_squarem$par[1:(L^2)], L, L)
  theta_MLE$kappa <- result_squarem$par[(L^2+1)]
  theta_MLE$LL <- as.numeric(names(result_squarem$par)[1])

  #success <- (result_squarem$convergence==0) #0 indicates convergence, 1 indicates failure to converge within maxiter
  numiter <- result_squarem$fpevals #number of parameter update steps (approximately 3x the number of SQUAREM iterations)

  ### Compute final posterior mean of subject ICs
  if(verbose) cat('Computing final posterior mean of subject ICs \n')
  MAP <- UpdateTheta_diagnosticICA.spatial(template_mean,
                                           template_var,
                                           meshes,
                                           BOLD,
                                           theta,
                                           C_diag,
                                           s0_vec_list,
                                           D_list,
                                           Dinv_s0_list,
                                           verbose=verbose,
                                           return_MAP=TRUE,
                                           ignore_determinant=TRUE)


  result <- list(group_probs = MAP$probs_z,
                 subjICmean=MAP$mu_s,
                 subjICcov=MAP$cov_s,
                 Omega = MAP$Omega_s_list,
                 theta_MLE=theta_MLE,
                 theta_path=theta_path,
                 numiter=numiter,
                 squarem = result_squarem,
                 template_mean = template_mean,
                 template_var=template_var)
  return(result)
}



#' @rdname EM_diagnosticICA
#'
EM_diagnosticICA.independent = function(template_mean, template_var, BOLD, theta0, C_diag, maxiter=100, epsilon=0.001, verbose=TRUE){

  nvox <- nrow(BOLD) #number of brain locations
  L <- ncol(template_mean[[1]]) #number of ICs
  G <- length(template_mean)
  if(L > nvox) stop('Cannot estimate more ICs than brain locations.')
  if(L != ncol(BOLD)) stop('Dimension-reduced BOLD data should have number of columns equal to number of ICs in templates')

  iter = 1
  theta = theta0
  success = 1
  for(g in 1:G) {
    num_smallvar <- sum(template_var[[g]] < 1e-6)
    if(num_smallvar>0){
      if(verbose) cat(paste0('Setting ',num_smallvar,' (',round(100*num_smallvar/length(template_var[[g]]),1),'%) very small variance values in group ',g,' template to ',1e-6,'.\n'))
      template_var[[g]][template_var[[g]] < 1e-6] = 1e-6 #to prevent problems when inverting covariance
    }
  }

  #make the group template variances equal (to avoid misclassification when one group has much larger variance)
  template_var_max <- template_var[[1]]
  for(g in 2:G){
    template_var_diff <- template_var[[g]] - template_var_max
    template_var_diff[template_var_diff < 0] <- 0 #places where max is already greater, do not change
    template_var_max <- template_var_max + template_var_diff #places where max is not greater, add to max
  }
  template_var_max[template_var_max < 1e-6] <- 1e-6
  # for(g in 1:G){
  #   template_var[[g]] <- template_var_max
  # }

  #get initial estimates of S to update pr_z
  MAP1 = UpdateTheta_diagnosticICA.independent(template_mean, template_var, template_var_max, BOLD, theta, C_diag, return_MAP=TRUE, verbose=verbose)
  dist1 <- colSums((MAP1$ICmean - template_mean[[1]])^2/(template_var_max))
  dist2 <- colSums((MAP1$ICmean - template_mean[[2]])^2/(template_var_max))
  #for 2 groups only:
  pr1 <- mean(dist1<dist2)
  pr2 <- 1-pr1
  pr_z = c(pr1, pr2)
  theta$pr_z <- pr_z
  print(paste0('Updated Initial Group Probabilities: ',paste(round(pr_z,3), collapse=', ')))


  err = 1000 #large initial value for difference between iterations
  while(err > epsilon){

    if(verbose) cat(paste0(' ~~~~~~~~~~~~~~~~~~~~~ ITERATION ', iter, ' ~~~~~~~~~~~~~~~~~~~~~ \n'))

    t00 <- Sys.time()
    theta_new = UpdateTheta_diagnosticICA.independent(template_mean, template_var, template_var_max, BOLD, theta, C_diag, return_MAP=FALSE, verbose=verbose)
    if(verbose) print(Sys.time() - t00)

    ### Compute change in parameters

    A_old = theta$A
    A_new = theta_new$A
    pr_old = theta$pr_z[1]
    pr_new = theta$pr_z[1]
    #2-norm = largest eigenvalue = sqrt of largest eigenvalue of AA'
    err = norm(as.vector(A_new - A_old), type="2")/norm(as.vector(A_old), type="2")
    err2 = abs(pr_new - pr_old)
    change = format(err, digits=3, nsmall=3)
    change2 = format(err2, digits=3, nsmall=3)
    if(verbose) cat(paste0('Iteration ',iter, ': Difference is ',change,' for A\n'))
    if(verbose) cat(paste0('Iteration ',iter, ': Difference is ',change2,' for pr1\n'))
    err = max(err, err2)

    ### Move to next iteration
    theta <- theta_new
    iter = iter + 1
    if(iter > maxiter){
      success = 0;
      warning(paste0('Failed to converge within ', maxiter,' iterations'))
      break() #exit loop
    }
  }

  MAP = UpdateTheta_diagnosticICA.independent(template_mean, template_var, template_var_max, BOLD, theta, C_diag, return_MAP=TRUE, verbose=verbose)

  result <- list(group_probs = MAP$group_probs,
                 subjICmean=MAP$ICmean,
                 subjICvar=MAP$ICvar,
                 theta_MLE=theta,
                 success_flag=success,
                 error=err,
                 numiter=iter-1,
                 template_mean = template_mean,
                 template_var = template_var)

  return(result)
}




#' @name UpdateTheta_diagnosticICA
#' @rdname UpdateTheta_diagnosticICA
#'
#' @title Parameter Estimates in EM Algorithm for Diagnostic ICA Model
#'
#' @param template_mean (A list of G matrices, each VxL) template mean estimates for each group 1 to G
#' @param template_var (A list of G matrices, each VxL) template variance estimates for each group 1 to G
#' @param template_var_max The maximum of the template variance across group
#' @param BOLD  (VxL matrix) dimension-reduced fMRI data
#' @param meshes Either NULL (assume spatial independence) or a list of objects of type \code{templateICA_mesh} created by \code{make_mesh()} (each list element corresponds to one brain structure)
#' @param theta (list) current parameter estimates
#' @param C_diag (Lx1) diagonal elements of residual covariance after dimension reduction
#' @param s0_vec_list List of vectorized template means (one for each group 1 to G)
#' @param D_list List of sparse diagonal matrices of template standard deviations (one for each group 1 to G)
#' @param Dinv_s0_list List of D^{-1} times s0_vec (one for each group 1 to G)
#' @param verbose If TRUE, display progress of algorithm
#' @param return_MAP If TRUE, return the posterior mean and precision of the latent fields and group membership probabilities instead of the parameter estimates
#' @param update Which parameters to update. Either "all", "A" or "kappa".
#' @param ignore_determinant For spatial model only. If TRUE, ignore the normalizing constant in p(y|z) when computing posterior probabilities of z.
#'
#' @return An updated list of parameter estimates, theta, OR if return_MAP=TRUE, the posterior mean and precision of the latent fields and posterior probabilities of group membership
#'
NULL

#' @rdname UpdateTheta_diagnosticICA
#' @importFrom stats optimize
#' @importFrom INLA inla.qsolve inla.qinv inla.setOption
#' @import Matrix
UpdateTheta_diagnosticICA.spatial = function(template_mean, template_var, meshes, BOLD, theta, C_diag, s0_vec_list, D_list, Dinv_s0_list, verbose=FALSE, return_MAP=FALSE, update=c('all','kappa','A'), ignore_determinant=TRUE){

  nvox = nrow(BOLD)
  L = ncol(BOLD)
  G = length(template_mean)

  #initialize new parameter values
  theta_new = theta

  #which parameters to update
  update_A <- (update[1] == 'all' | update[1] =='A')
  update_kappa <- (update[1] == 'all' | update[1] =='kappa')
  if(update_A) theta_new$A <- NA
  if(update_kappa) theta_new$kappa <- NA
  if(!update_A & !update_kappa & !return_MAP) stop('Please indicate at least one parameter to update (e.g. update="A") or set return_MAP=TRUE to return MAP estimates of S based on current parameter values.')

  #1. Compute mu_sy by solving for x in the system of equations Omega*x = m
  #   Compute m and Omega (using previous parameter estimates)
  #   Ingredients for Omega and m:
  #     s0_vec = "s_0" (long vector organized by ICs)
  #     y (long vector organized by locations)
  #     D, big diagonal matrix of template standard deviations
  #     R_inv, prior inverse correlation matrix for data locations (big block diagonal matrix, each block is sparse-ish)
  #     P (permutation matrix)

  ##########################################
  ### POSTERIOR MOMENTS OF s|z
  ##########################################

  if(verbose) cat('Computing Posterior Moments of S \n')

  y_vec = as.vector(t(BOLD)) #grouped by locations

  ### TO DO FOR WHOLE BRAIN MODEL: LOOP OVER MESHES, COMPUTE R_INV FOR EACH, AND COMBINE AS A BLOCK DIAGONAL.
  ### NEED TO VERIFY THAT LOCATIONS ARE GROUPED BY BRAIN STRUCTURE

  if(verbose) cat('...posterior precision \n') # less than 1 sec
  #Compute SPDE matrices (F, G, GFinvG) and Sigma_inv (LVxLV), a sparse block diagonal matrix
  stuff <- compute_R_inv(meshes, kappa=theta$kappa, C1=1/(4*pi))
  R_inv <- bdiag(stuff$Rinv)
  Fmat <- bdiag(stuff$Fmat)
  Gmat <- bdiag(stuff$Gmat)
  GFinvG <- bdiag(stuff$GFinvG)

  #set up P as a sparse matrix (see OneNote for illustration of this)
  P <- make_Pmat(L, nvox)

  ##########################################
  ### POSTERIOR PROBABILITIES OF z
  ##########################################

  #2. Compute posterior probabilities of group membership

  for_pr_zy <- rep(NA, G)
  pr_zy <- rep(NA, G)

  bigM <- bdiag_m(rep(list(theta$A), nvox)) #the matrix M_\otimes in the paper, changed notation from A to M
  bigM_inv <- bdiag_m(rep(list(solve(theta$A)), nvox)) #the matrix M_\otimes in the paper, changed notation from A to M
  bigC <- bdiag_m(rep(list(diag(C_diag)), nvox))
  bigC_inv <- bdiag_m(rep(list(diag(1/C_diag)), nvox))
  Pt_Minv_C_Minvt_P <- t(P) %*% bigM_inv %*% bigC %*% t(bigM_inv) %*% P
  for(g in 1:G){
    #print(g)
    mu_y_g <- bigM %*% P %*% s0_vec_list[[g]] #group s0 by locations
    d_g <- y_vec - mu_y_g

    #clever computation of exponential part
    #~12 sec for L=16, V=5246
    #~50 sec for L=16, V=8800
    B_g <- bigM %*% P %*% D_list[[g]]
    E_g <- R_inv + t(B_g) %*% bigC_inv %*% B_g
    firstpart <- t(d_g) %*% bigC_inv %*% B_g
    lastpart <- inla.qsolve(Q = E_g, B=t(firstpart), method='solve')
    exppart_g <- t(d_g) %*% bigC_inv %*% d_g - firstpart %*% lastpart #gives same answer as direct computation

    #direct(ish) computation of determinant part (~20 min for L=16, V=5246)
    if(!ignore_determinant){
      I_mat <- Diagonal(nvox*L, 1)
      D_inv_g <- solve(D_list[[g]])
      #bigmat <- I_mat + R_inv %*% D_inv_g %*% Pt_Minv_C_Minvt_P %*% D_inv_g #how to compute determinant? cholesky first?
      bigmat2 <- R_inv + R_inv %*% D_inv_g %*% Pt_Minv_C_Minvt_P %*% D_inv_g %*% R_inv # <- symmetric! (can use Cholesky) det(bigmat) = det(bigmat2)/det(R_inv)
      system.time(chol_big <- chol(bigmat2, pivot=TRUE))
      log_det_R_inv <- determinant(R_inv[1:nvox,1:nvox])$modulus[1] #<1 sec

      detpart1 <- 2*nvox*determinant(theta$A)$modulus[1]
      detpart2 <- sum(log(template_var[[g]]))
      detpart3 <- 2*sum(log(diag(chol_big)))
      detpart4 <- 2*L*log_det_R_inv

      detpart_g <- detpart1 + detpart2 + detpart3 - detpart4
    } else {
      detpart_g <- 0
    }

    for_pr_zy[g] <- -0.5 * (detpart_g + exppart_g[1]) #log p(y|z=g)
  }

  for(g in 1:G){
    pr_zy_inv_g <- sum(exp(for_pr_zy-for_pr_zy[g])) # 1/pr(z=1|y) = ( (e^M1 + e^M2 + e^M3) / e^M1 ) = e^(M1-M1) + e^(M2-M1) + e^(M3-M1) = 1 + e^(M2-M1) + e^(M3-M1)  <-- If any e^(Mk-Mg) Inf, the inverse will be zero so p(z=g|y)=0
    pr_zy[g] <- 1/pr_zy_inv_g
  }

  #fix numerical issues with very small values (values very close to 1 are numerically equal to 1, while values very close to zero are not)
  if(any(pr_zy==1)){
    pr_zy[pr_zy!=1] <- 0
  }

  if(verbose) cat(paste0('Posterior probabilities of z: ',paste(round(pr_zy,3), collapse=', '),'\n'))


  #1. Compute mu_s for each group
  if(verbose) cat('...posterior mean \n')
  #20 seconds per group with pardiso! (1 hour without) for L=16, V=5246
  #60 seconds for L=16, V=8800
  mu_s_list <- m_vec_list <- Omega_list <- Omega_inv_m_list <- vector('list', length=G)
  for(g in 1:G){
    if(pr_zy[g]==0) next
    stuff <- compute_mu_s(y_vec, D_list[[g]], Dinv_s0_list[[g]], R_inv, theta, P, C_diag)
    mu_s_list[[g]] <- stuff$mu
    m_vec_list[[g]] <- stuff$m
    Omega_list[[g]] <- stuff$Omega
    Omega_inv_m_list[[g]] <- stuff$Omega_inv_m
    if(verbose) print(summary(as.vector(stuff$mu)))
  }


  if(return_MAP){
    #compute posterior mean and covariance of s
    mu_s <- array(0, dim=dim(mu_s_list[[g]]))
    first_group <- TRUE
    for(g in 1:G){
      if(pr_zy[g]==0) next
      #sparsity structure of Omega_g is common across g, because D_g is only difference
      cov_sg <- (D_list[[g]] %*% inla.qinv(Omega_list[[g]]) %*% D_list[[g]])*pr_zy[g] #only calculates diagonal elements and non-zero elements of precision (correspond to neighboring locations in mesh)
      if(first_group) {
        mu_s <- mu_s_list[[g]]*pr_zy[g]
        cov_s <- cov_sg
      } else {
        mu_s <- mu_s + mu_s_list[[g]]*pr_zy[g]
        cov_s <- cov_s + cov_sg
      }
      first_group <- FALSE
    }
    return(list(probs_z = pr_zy, mu_s = mu_s, cov_s = cov_s, Omega_s_list = Omega_list))
    stop()
  }

  ##########################################
  ### MLE OF A
  ##########################################

  if(update_A){
    if(verbose) cat('Updating A \n')
    t00 <- Sys.time()


    #Compute first matrix appearing in A-hat
    A_part1 <- matrix(0, nrow=L, ncol=L)
    for(g in 1:G){
      if(pr_zy[g]==0) next
      Pmu_g <- P %*% mu_s_list[[g]]
      yPmu_g <- matrix(0, nrow=L, ncol=L)
      for(v in 1:nvox){
        inds_v <- (1:L) + (v-1)*L
        y_v <- y_vec[inds_v]
        Pmu_g_v <- Pmu_g[inds_v]
        yPmu_g_v <- outer(y_v, Pmu_g_v)
        yPmu_g <- yPmu_g + yPmu_g_v #sum up to later divide by nvox to average
      }
      A_part1 <- A_part1 + yPmu_g*pr_zy[g]
    }

    #Compute second matrix appearing in A-hat
    A_part2 <- matrix(0, nrow=L, ncol=L)
    #set up sparse indicator matrix for needed elements of Omega_PP_g^{-1}
    oneblock <- matrix(1, nrow=L, ncol=L)
    ind_mat <- bdiag_m(rep(list(oneblock), nvox))
    attributes(ind_mat)$x <- 0*attributes(ind_mat)$x
    #loop: ~7 mins per group for L=16, V=5246
    for(g in 1:G){
      if(pr_zy[g]==0) next
      Pmu_g <- P %*% mu_s_list[[g]]
      Omega_PP_g <- P %*% Omega_list[[g]] %*% t(P)
      D_PP_g <- P %*% D_list[[g]] %*% t(P)

      if(verbose) cat('..Calculating only non-sparse covariance terms \n')

      #augment Omega_PP with additional non-sparse locations needed for A_part2
      Omega_PP_g_aug <- Omega_PP_g + ind_mat
      system.time(Omega_PP_g_inv <- inla.qinv(Omega_PP_g_aug))

      T_mat_g <- matrix(0, L, L)
      for(v in 1:nvox){
        inds_v <- (1:L) + (v-1)*L

        #T1 matrix
        Omega_PP_g_inv_vv <- Omega_PP_g_inv[inds_v,inds_v]
        T1_g_vv <- D_PP_g[inds_v,inds_v] %*% Omega_PP_g_inv_vv %*% D_PP_g[inds_v,inds_v]

        #T2 matrix
        Pmu_g_v <- Pmu_g[inds_v]
        T2_g_vv <- outer(Pmu_g_v, Pmu_g_v)

        #T(v,v)
        T_g_vv <- T1_g_vv + T2_g_vv
        T_mat_g <- T_mat_g + T_g_vv
      }
      A_part2 <- A_part2 + T_mat_g*pr_zy[g]
    }

    A_hat <- A_part1 %*% solve(A_part2)
    theta_new$A <- as.matrix(A_hat)

    #compute log likelihood Q1(Theta) in paper

    #first part of Q1: sum_v 2 y(v)' C^{-1} A t(v) = sum_v Trace{ 2 C^{-1} A t(v) y(v)' } = Trace{ 2 C^{-1} A sum_v [ t(v) y(v)' ] }, where sum_v [ t(v) y(v)' ] = A_part1'
    LL1_part1 <- sum(diag( 2 * diag((1/C_diag)) %*% A_hat %*% t(A_part1) ))

    #second part of Q1: sum_v Trace{ A' C^{-1} A T(v,v) } = Trace{ A' C^{-1} A sum_v T(v,v) }, where sum_v T(v,v) is A_part2
    LL1_part2 <- sum(diag( t(A_hat) %*% diag((1/C_diag)) %*% A_hat %*% A_part2 ))

    theta_new$LL[1] <- LL1_part1 - LL1_part2

    if(verbose) print(Sys.time() - t00)

  }


  ##########################################
  ### E-STEP for kappa
  ##########################################

  if(update_kappa){

    #1 min for L=16, V=8800
    if(verbose) cat('Updating kappa  \n')
    t00 <- Sys.time()

    # Likelihood in terms of kappa's.
    # LL2 = LL2_part1 + sum_g pr(z=g) * [LL2_part2_g + LL2_part3_g]
    # LL2_part1 = L * log(det(R_inv))
    # LL2_part2_g = sum_ell [ Tr(R_inv * Omega_inv_ell) + Tr(R_inv * W_hat_ell) ]
    # LL2_part3_g = sum_ell [ s0_ell_g' * Dinv_g * R_inv * v_hat_ell_g ]
    #             v_g = 2 Omega_inv_g * m_g - Dinv_g * s0_g

    # Tr(R_ell_inv * Omega_inv_ell) --> Use inla.inv to compute necessary elements (non-zeros in R_ell_inv) of Omega_inv_ell for ell=1,...,L
    # Tr(R_ell_inv * W_hat_ell) --> W_hat_ell = outer(Omega_inv*m,Omega_inv*m)_ll, where Omega_inv*m is known. Just calculate necessary elements (non-zeros in R_ell_inv) of W_hat_ell for ell=1,...,L

    if(verbose) cat('..Computing trace terms in log-likelihood of kappa \n')

    ### SET UP FOR PART 2

    #1. Generate Monte Carlo samples for estimation of Omega_hat^{-1}_ll

    # if(verbose) cat(paste0('....drawing ',nsamp,' Monte Carlo samples for covariance estimation \n'))
    #
    # nsamp <- 5000 #number of samples
    # musamp <- inla.qsample(nsamp, Q = Omega, b=rep(0, V*L), mu=rep(0, V*L), num.threads=8) #sample from N(0, Omega^(-1)) to estimate Omega^(-1) (diagonal blocks only)
    #9 sec for 500 samples with V*L = 2530*3 = 12,642

    #2. Determine non-zero terms of R_ell_inv

    if(verbose) cat('....calculating only non-sparse covariance terms \n')

    #set up estimation of only terms necessary for trace computation
    nonzero_Rinv <- as.matrix(1*(R_inv[1:nvox,1:nvox] != 0)) #this is R_ell^{-1} for a single IC ell
    #diag(nonzero_Rinv) <- 1 #make sure all diagonal elements estimated (required for Trace 1)
    nonzero_cols <- which(R_inv[1:nvox,1:nvox] != 0, arr.ind = TRUE)[,2] #column indices of all non-zero locations in R_ell^{-1}

    #augment Omega with additional non-sparse locations
    ind_mat <- 1*(R_inv[1:nvox,1:nvox] != 0)
    attributes(ind_mat)$x <- 0*attributes(ind_mat)$x
    new_inds <- bdiag(rep(list(ind_mat), L))

    #Set up vectors needed for computation of only necessary elements of W_hat

    #x <- matrix(1, nrow=1, ncol=nvox) #placeholder for x in M = xx'
    #bigx_left <- diag(1, nvox) #KhatriRao(diag(1, nvox), x)
    #bigx_right <-  KhatriRao(nonzero_Rinv, x))
    #inds_left_x <- which(diag(1, nvox) != 0, arr.ind=TRUE) #1,1 2,2 3,3 etc. (bigx_left = I_V) <- just used to make a diagonal matrix below, don't actually need this
    inds_right_x <- which(nonzero_Rinv != 0, arr.ind=TRUE) #which(bigx_right != 0, arr.ind=TRUE)

    for(g in 1:G){

      if(pr_zy[g]==0) next

      #Omega_aug_g <- Omega_list[[g]] +  new_inds
      #compute elements of Omega_inv that are non-sparse in Omega and in R^{-1}
      Omega_inv_g <- inla.qinv(Omega_list[[g]]) #don't need any additional elements since Omega and R_inv have same sparsity structure

      if(verbose) cat('....computing necessary elements of RHS matrices in trace terms \n')
      #15 sec for L=16, nvox=5500

      #Compute Omega_inv_ll_g + W_hat_ll_g

      OplusW <- vector('list', length=G)
      for(ell in 1:L){
        inds_ell <- (1:nvox) + (ell-1)*nvox

        # COMPUTE OMEGA_INV_g[ell,ell] (NECESSARY ENTRIES ONLY)
        Omega_inv_ll_g <- Omega_inv_g[inds_ell,inds_ell]

        # COMPUTE W_hat_g[ell,ell] = Omega_inv_m_g[ell] * Omega_inv_m_g[ell]' (NECESSARY ENTRIES ONLY)
        # T_ell = matrix that selects the ellth block of size nvox from a matrix with nvox*L rows
        e_ell <- Matrix(0, nrow=1, ncol=L, sparse=TRUE)
        e_ell[1,ell] <- 1
        T_ell <- kronecker(e_ell, Diagonal(nvox))
        Omega_inv_m_ell_g <- T_ell %*% Omega_inv_m_list[[g]] #lth block of the vector Omega^{-1}m for group g

        #left-multiplication matrix
        bigx_left_ell_g <- Diagonal(nvox, x = as.vector(Omega_inv_m_ell_g)) #sparseMatrix(i = inds_left_x[,1], j = inds_left_x[,2], x = as.vector(Omega_inv_m_ell_g))
        #right-multiplication matrix
        x_repcols <- t(Omega_inv_m_ell_g)[,nonzero_cols] #repeat elements of Omega_inv_m_ell_g to fill in all non-zero cells of R^{-1}
        bigx_right_ell_g <- sparseMatrix(i = inds_right_x[,1], j = inds_right_x[,2], x = as.vector(x_repcols))
        #multiply together
        W_hat_ll_g <- t(bigx_left_ell_g) %*% bigx_right_ell_g

        # COMBINE OMEGA_INV[ell,ell] AND W_hat[ell,ell] --> two trace terms involving R_ell_inv can be combined

        OplusW_ll_g <- Omega_inv_ll_g + W_hat_ll_g
        if(ell==1) OplusW_g <- OplusW_ll_g else OplusW_g <- OplusW_g + OplusW_ll_g
      }
      OplusW[[g]] <- OplusW_g
    }


    ### SET UP FOR PART 3

    u_list <- v_list <- vector('list', length=G)
    for(g in 1:G){
      if(pr_zy[g]==0) next
      u_list[[g]] <- Dinv_s0_list[[g]]
      v_list[[g]] <- 2*Omega_inv_m_list[[g]] - u_list[[g]]
    }

    # NUMERICALLY SOLVE FOR MLE OF KAPPA

    if(verbose) cat("..Performing numerical optimization \n")
    #15 seconds for V=5246, L=16
    Amat <- bdiag(lapply(meshes, function(x) return(x$A)))
    system.time(kappa_opt <- optimize(LL2_kappa_diagnosticICA, lower=0, upper=5, maximum=TRUE,
                                      Amat=Amat, Fmat=Fmat, Gmat=Gmat, GFinvG=GFinvG, OplusW=OplusW, u=u_list, v=v_list, Q=L, pr_zy=pr_zy))
    LL2 <- kappa_opt$objective
    kappa_opt <- rep(kappa_opt$maximum, L)

    theta_new$kappa <- kappa_opt
    theta_new$LL[2] <- LL2

    if(verbose) print(Sys.time() - t00)

  }

  # RETURN NEW PARAMETER ESTIMATES

  return(theta_new)

}
#' @rdname UpdateTheta_diagnosticICA
UpdateTheta_diagnosticICA.independent = function(template_mean, template_var, template_var_max, BOLD, theta, C_diag, return_MAP=FALSE, verbose=TRUE){

  L <- ncol(BOLD)
  nvox <- nrow(BOLD)
  G <- length(template_mean)

  #initialize new objects
  theta_new = list(A = matrix(NA, L, L))
  A_part1 = A_part2 = matrix(0, L, L) #two parts of product for A-hat (construct each looping over voxels)

  A <- theta$A
  C <- diag(C_diag)
  C_inv <- diag(1/(C_diag))
  At_Cinv <- t(A) %*% C_inv
  At_Cinv_A <- At_Cinv %*% A

  ##########################################
  ### POSTERIOR PROBABILITIES OF z
  ##########################################

  # if(verbose) cat('Computing posterior probabilities of z (group membership) \n')
  #
  # exp_part <- rep(NA, G)
  # pr_zy <- rep(NA, G)
  # for(g in 1:G){
  #
  #   #exp_part1 <- L*nvox*log(2*pi)
  #   exp_part1 <- 0 #irrelevant as long as noninformative prior on z=g
  #   exp_part2 <- 0 #sum(log(template_var[[g]])) #sum over v,ell
  #
  #   exp_part3 <- 0
  #   for(v in 1:nvox){
  #
  #     y_v = BOLD[v,]
  #     s0_gv = template_mean[[g]][v,]
  #
  #     #mean and cov of y(v)|z
  #     mu_yz_v <- A %*% s0_gv
  #     #Sigma_yz_v <- A %*% diag(template_var[[g]][v,]) %*% t(A) + C
  #     Sigma_yz_v <- A %*% diag(template_var_max[v,]) %*% t(A) + C
  #
  #     exp_part3_v <- t(y_v - mu_yz_v) %*% solve(Sigma_yz_v) %*% (y_v - mu_yz_v)
  #     exp_part3 <- exp_part3 + exp_part3_v
  #
  #   }
  #
  #   exp_part[g] <- -0.5 * (exp_part1 + exp_part2 + exp_part3[1])
  # }
  #
  # for(g in 1:G){
  #   pr_zy_inv_g <- sum(exp(exp_part-exp_part[g])) # ( (e^M1 + e^M2 + e^M3) / e^M1 ) = e^(M1-M1) + e^(M2-M1) + e^(M3-M1) = 1 + e^(M2-M1) + e^(M3-M1)  <-- If any e^(Mk-Mg) Inf, the inverse will be zero so p(z=g|y)=0
  #   pr_zy[g] <- 1/pr_zy_inv_g
  # }


  # #randomly split locations into buckets to avoid infinite exponentials (this creates pseudo-independent observations)
  # random_vox <- sample(1:nvox, nvox, replace=FALSE) #randomly reorder voxels
  # size_buckets <- 1000
  # num_buckets <- ceiling(nvox/size_buckets) #number of buckets
  #
  # pr_zy_buckets <- matrix(NA, num_buckets, G)
  # for(b in 1:num_buckets){
  #
  #   first_vox_b <- size_buckets*(b-1) + 1
  #   last_vox_b <- min(first_vox_b + size_buckets - 1, nvox)
  #   vox_b <- random_vox[first_vox_b:last_vox_b]
  #
  #   exp_part_v <- rep(NA, G)
  #   pr_zy_v <- rep(NA, G)
  #   for(g in 1:G){
  #
  #     #exp_part1 <- L*nvox*log(2*pi)
  #     exp_part1 <- 0 #irrelevant as long as noninformative prior on z=g
  #     exp_part2 <- sum(log(template_var[[g]][vox_b,])) #sum over v,ell
  #
  #     exp_part3 <- 0
  #     for(v in vox_b){
  #
  #       y_v = BOLD[v,]
  #       s0_gv = template_mean[[g]][v,]
  #
  #       #mean and cov of y(v)|z
  #       mu_yz_v <- A %*% s0_gv
  #       Sigma_yz_v <- A %*% diag(template_var[[g]][v,]) %*% t(A) + C
  #
  #       exp_part3_v <- t(y_v - mu_yz_v) %*% solve(Sigma_yz_v) %*% (y_v - mu_yz_v)
  #       exp_part3 <- exp_part3 + exp_part3_v
  #
  #     }
  #
  #     exp_part_v[g] <- -0.5 * (exp_part1 + exp_part2 + exp_part3[1])
  #   }
  #
  #   for(g in 1:G){
  #     pr_zy_inv_g <- sum(exp(exp_part_v-exp_part_v[g])) # ( (e^M1 + e^M2 + e^M3) / e^M1 ) = e^(M1-M1) + e^(M2-M1) + e^(M3-M1) = 1 + e^(M2-M1) + e^(M3-M1)  <-- If any e^(Mk-Mg) Inf, the inverse will be zero so p(z=g|y)=0
  #     pr_zy_v[g] <- 1/pr_zy_inv_g
  #   }
  #   pr_zy_buckets[b,] <- pr_zy_v
  #
  # }
  # pr_zy <- colMeans(pr_zy_buckets)
  #
  #fix numerical issues with very small values (values very close to 1 are numerically equal to 1, while values very close to zero are not)
  # if(any(pr_zy==1)){
  #   pr_zy[pr_zy!=1] <- 0
  # }

  # if(verbose) print(paste(round(pr_zy,3), collapse=','))

  # if(is.infinite(exp(max(exp_part)-min(exp_part)))) {
  #   #this is for two groups, need to generalize for G>2
  #   pr_zy[which.max(exp_part)] <- 1
  #   pr_zy[which.min(exp_part)] <- 0
  # } else {
  #   exp_part <- exp_part - mean(exp_part) #this minimizes the most extreme magnitude in exp_part (to avoid )
  #   prod_pr_yz <- exp(exp_part)
  #   #compute p(z=g|y), g=1,...,G
  #   denom <- sum((1/G)*prod_pr_yz) #can change 1/G for unequal prior probs
  #   pr_zy <- (1/G)*prod_pr_yz / denom
  # }

  ##########################################
  ### POSTERIOR MOMENTS OF s
  ##########################################

  #instead of posterior probabilities of z, use MLE of group probabilities
  pr_zy <- theta$pr_z

  if(verbose) cat('Computing posterior moments of s (IC maps) \n')

  mu_sy <- array(NA, dim=c(nvox,L))
  mu_mut_sy <- array(NA, dim=c(nvox,L,L))
  var_sy <- array(NA, dim=c(nvox,L)) #only needed if return_MAP=TRUE
  for(v in 1:nvox){

    mu_sy_v <- rep(0, L)
    mu_mut_sy_v <- Sigma_sy_v <- array(0, dim=c(L,L))
    for(g in 1:G){

      if(pr_zy[g]==0) next
      D_gv_inv <- diag(1/template_var[[g]][v,])
      s0_gv <- template_mean[[g]][v,]
      Sigma_syz_gv <- solve(At_Cinv_A + D_gv_inv) #posterior covariance of s|z=g
      mu_syz_gv <- Sigma_syz_gv %*% (At_Cinv %*% BOLD[v,] + D_gv_inv %*% s0_gv) #posterior mean of s|z=g
      mu_mut_syz_gv <- mu_syz_gv %*% t(mu_syz_gv) + Sigma_syz_gv #posterior second moment of s|z=g

      mu_sy_v <- mu_sy_v + pr_zy[g]*mu_syz_gv #weighted sum over g=1..G --> posterior mean of s
      mu_mut_sy_v <- mu_mut_sy_v + pr_zy[g]*mu_mut_syz_gv
      Sigma_sy_v <- Sigma_sy_v + pr_zy[g]*Sigma_syz_gv #weighted sum over g=1..G --> posterior second moment of s

    }

    mu_sy[v,] <- mu_sy_v
    mu_mut_sy[v,,] <- mu_mut_sy_v
    var_sy[v,] <- diag(Sigma_sy_v) #only needed if return_MAP=TRUE

  }

  #check whether mu_mut_sy = mu_sy %*% t(mu_sy) + Sigma_sy

  if(return_MAP){
    result <- list(#group_probs = pr_zy,
                   ICmean = mu_sy,
                   ICvar = var_sy)
    return(result)
  }

  ##########################################
  ### MLE OF A
  ##########################################

  if(verbose) cat('Updating A \n')

  A_part1 <- A_part2 <- matrix(0, L, L)
  for(v in 1:nvox){

    y_v = BOLD[v,]

    ##########################################
    ### M-STEP FOR A: CONSTRUCT PARAMETER ESTIMATES
    ##########################################

    A_part1 = A_part1 + y_v %*% t(mu_sy[v,]) #LxL
    A_part2 = A_part2 + mu_mut_sy[v,,] #LxL

  }

  A_hat = (A_part1 %*% solve(A_part2))

  ##########################################
  ### ESTIMATE OF PR_Z
  ##########################################

  dist1 <- colSums((mu_sy - template_mean[[1]])^2/(template_var_max))
  dist2 <- colSums((mu_sy - template_mean[[2]])^2/(template_var_max))
  #for 2 groups only:
  pr1 <- mean(dist1<dist2)
  pr2 <- 1-pr1
  pr_z = c(pr1, pr2)
  print(paste0('Updated Group Probabilities: ',paste(round(pr_z,3), collapse=', ')))


  # RETURN NEW PARAMETER ESTIMATES

  theta_new$A <- A_hat
  theta_new$pr_z <- pr_z
  return(theta_new)
}

#' Computes part of log-likelihood involving kappa (or kappa_q) for numerical optimization
#'
#' @param kappa Value of kappa for which to compute log-likelihood
#' @param Amat Mesh projection matrix
#' @param Fmat Matrix used in computation of SPDE precision
#' @param Gmat Matrix used in computation of SPDE precision
#' @param GFinvG Matrix used in computation of SPDE precision
#' @param OplusW LIST OF sparse matrices containing estimated values of RHS of trace in part 2 of log-likelihood.  (one list element for each group g)
#' @param u LIST OF vectors needed for part 3 of log-likelihood (one list element for each group g)
#' @param v LIST OF vectors needed for part 3 of log-likelihood (one list element for each group g)
#' @param C1 For the unit variance case, \eqn{\tau^2 = C1/\kappa^2}, where \eqn{C1 = 1/(4\pi)} when \eqn{\alpha=2}, \eqn{\nu=1}, \eqn{d=2}
#' @param Q Equal to the number of ICs
#' @param pr_zy Current posterior probabilities of group membership z
#'
#' @import Matrix
#' @return Value of log-likelihood at logkappa
#'
#' @details This is the function to be maximized in order to determine the MLE for \eqn{\kappa} or the \eqn{\kappa_q}'s in the M-step of the EM algorithm in spatial
#' template ICA.  In the model where \eqn{\kappa_q} can be different for each IC \eqn{q}, the optimization function factorizes over the \eqn{\kappa_q}'s.  This function computes
#' the value of the part of the optimization function pertaining to one of the \eqn{\kappa_q}'s.
#'
LL2_kappa_diagnosticICA <- function(kappa, Amat, Fmat, Gmat, GFinvG, OplusW, u, v, C1 = 1/(4*pi), Q, pr_zy){

  G <- length(u)

  #get data and nodata indices of mesh vertices
  nmesh = ncol(Amat) #number of mesh locations
  inds_data = which(colSums(Amat) > 0)
  inds_nodata = setdiff(1:nmesh, inds_data)

  #COMPUTE R_q_inv FOR CURRENT VALUE OF kappa

  Qmat = C1*(kappa^2 * Fmat + 2 * Gmat + kappa^(-2) * GFinvG)
  Q11 = Qmat[inds_data,inds_data] # = Amat %*% Qmat %*% t(Amat)
  if(length(inds_nodata)>0){ #marginalize out non-data locations in mesh
    Q12 = Qmat[inds_data, inds_nodata]
    Q21 = Qmat[inds_nodata, inds_data]
    Q22 = Qmat[inds_nodata,inds_nodata]
    Q22_inv <- solve(Q22)
    R_inv = Q11 - (Q12 %*% Q22_inv %*% Q21)
  } else { #no non-data locations in mesh
    R_inv = Q11
  }


  #CALCULATE LOG LIKELIHOOD FOR KAPPA IN 3 PARTS


  ### PART 1: log(det(R_q_inv))

  chol_Rinv <- chol(R_inv) #R_inv = chol_Rinv' * chol_Rinv
  LL2_part1 <- 2*Q*sum(log(diag(chol_Rinv))) #log determinant of R_q_inv
  LL2_total <- LL2_part1

  for(g in 1:G){

    if(pr_zy[g]==0) next

    ### PART 2: Trace terms
    #OplusW = sum_q Omega_inv_qq + W_hat_qq
    LL2_part2_g <- sum(R_inv * OplusW[[g]]) #equivalent to sum_q(diag(R_inv_qq %*% OplusW[q])) and FAST (note: R_inv_qq symmetric)

    ### PART 3: u' * R_inv * v
    LL2_part3_g <- t(u[[g]]) %*% bdiag(rep(list(R_inv), Q)) %*% v[[g]]

    LL2_total <- LL2_total - LL2_part2_g*pr_zy[g] + LL2_part3_g[1]*pr_zy[g]

  }

  return(as.numeric(LL2_total))

}


#' Helper function for SQUAREM for estimating parameters
#'
#' @param theta_vec Vector of initial parameter values
#' @param template_mean Passed to UpdateTheta_diagnosticICA function
#' @param template_var Passed to UpdateTheta_diagnosticICA function
#' @param meshes Passed to UpdateTheta_diagnosticICA function
#' @param BOLD Passed to UpdateTheta_diagnosticICA function
#' @param C_diag Passed to UpdateTheta_diagnosticICA function
#' @param s0_vec_list Passed to UpdateTheta_diagnosticICA function
#' @param D_list Passed to UpdateTheta_diagnosticICA function
#' @param Dinv_s0_list Passed to UpdateTheta_diagnosticICA function
#' @param verbose Passed to UpdateTheta_diagnosticICA function
#'
#' @return Vector of updated parameter values
#'
UpdateThetaSQUAREM_diagnosticICA <- function(theta_vec,
                                             template_mean,
                                             template_var,
                                             meshes,
                                             BOLD,
                                             C_diag,
                                             s0_vec_list,
                                             D_list,
                                             Dinv_s0_list,
                                             verbose){

  L = ncol(template_mean[[1]])

  #convert theta vector to list format
  A <- theta_vec[1:(L^2)]
  kappa <- rep(theta_vec[(L^2+1)],L)
  theta <- list(A = matrix(A, nrow=L, ncol=L),
                kappa = kappa)

  #update theta parameters
  if(verbose) cat('~~~~~~~~~~~ UPDATING PARAMETER ESTIMATES ~~~~~~~~~~~ \n')
  theta_new = UpdateTheta_diagnosticICA.spatial(template_mean,
                                                template_var,
                                                meshes,
                                                BOLD,
                                                theta,
                                                C_diag,
                                                s0_vec_list,
                                                D_list,
                                                Dinv_s0_list,
                                                verbose=verbose,
                                                return_MAP=FALSE,
                                                update='all',
                                                ignore_determinant=TRUE)
  #convert theta_new list to vector format
  #theta_new$A <- as.matrix(theta_new$A)
  theta_new_vec <- unlist(theta_new[1:2]) #everything but LL
  names(theta_new_vec)[1] <- sum(theta_new$LL)
  return(theta_new_vec)
}

#' Helper function for SQUAREM for extracting the negative of the log likelihood
#'
#' @param theta_vec Vector of current parameter values
#' @param template_mean Not used, but squarem will return error without
#' @param template_var  Not used, but squarem will return error without
#' @param meshes  Not used, but squarem will return error without
#' @param BOLD  Not used, but squarem will return error without
#' @param C_diag  Not used, but squarem will return error without
#' @param s0_vec_list  Not used, but squarem will return error without
#' @param D_list  Not used, but squarem will return error without
#' @param Dinv_s0_list  Not used, but squarem will return error without
#' @param verbose  Not used, but squarem will return error without
#'
#' @return Negative log-likelihood given current values of parameters
#'
LL_SQUAREM_diagnosticICA <- function(theta_vec,
                                     template_mean,
                                     template_var,
                                     meshes,
                                     BOLD,
                                     C_diag,
                                     s0_vec_list,
                                     D_list,
                                     Dinv_s0_list,
                                     verbose){

  LL <- names(theta_vec)[1]
  #print(LL)
  LL <- as.numeric(LL)
  print(LL)
  return(-1*LL)

}




