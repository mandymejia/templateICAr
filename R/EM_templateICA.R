#' @name EM_templateICA
#' @rdname EM_templateICA
#'
#' @title EM Algorithms for Template ICA Models
#'
#' @param template_mean (\eqn{V \times Q} matrix) mean maps for each IC in template,
#'  where \eqn{Q} is the number of ICs, \eqn{V=nvox} is the number of data locations.
#' @param template_var  (\eqn{V \times Q} matrix) between-subject variance maps for each IC in template
#' @param meshes \code{NULL} for spatial independence model, otherwise a list of
#'  objects of class "templateICA_mesh" containing the triangular mesh (see
#'  \code{\link{make_mesh}}) for each brain structure.
#' @param BOLD  (\eqn{V \times Q} matrix) dimension-reduced fMRI data
#' @param theta0 (list) initial guess at parameter values: A (\eqn{QxQ} mixing matrix),
#'  nu0_sq (residual variance from first level) and (for spatial model only)
#'  kappa (SPDE smoothness parameter for each IC map)
#' @param C_diag (\eqn{Qx1}) diagonal elements of matrix proportional to
#'  residual variance.
# @param common_smoothness If \code{TRUE}, use the common smoothness version
#'  of the spatial template ICA model, which assumes that all IC's have the
#'  same smoothness parameter, \eqn{\kappa}
#' @param maxiter Maximum number of EM iterations. Default: 100.
#' @param usePar Parallelize the computation over voxels? Default: \code{FALSE}. Can be the number of cores
#'  to use or \code{TRUE}, which will use the number on the PC minus two. Not implemented yet for spatial
#'  template ICA.
#' @param epsilon Smallest proportion change between iterations. Default: 0.01.
#' @param verbose If \code{TRUE}, display progress of algorithm. Default: \code{FALSE}.
#'
#' @return  A list: theta (list of final parameter estimates), subICmean
#'  (estimates of subject-level ICs), subICvar (variance of subject-level ICs,
#'  for non-spatial model) or subjICcov (covariance matrix of subject-level ICs,
#'  for spatial model -- note that only diagonal and values for neighbors are
#'  computed), and success (flag indicating convergence (\code{TRUE}) or not
#'  (\code{FALSE}))
#'
#' @details \code{EM_templateICA.spatial} implements the expectation-maximization
#'  (EM) algorithm described in Mejia et al. (2019+) for estimating the
#'  subject-level ICs and unknown parameters in the template ICA model with
#'  spatial priors on subject effects.
#'
#'  In both models, if original fMRI timeseries has covariance
#'  \eqn{\sigma^2 I_T}, the prewhitened timeseries achieved by premultiplying
#'  by (\eqn{QxT}) matrix \eqn{H} from PCA has diagonal covariance
#'  \eqn{\sigma^2HH'}, so C_diag is \eqn{diag(HH')}.
#'
#'
NULL

#' @rdname EM_templateICA
# @importFrom INLA inla.spde2.matern inla.qsolve
#' @importFrom Matrix Diagonal
#' @importFrom SQUAREM squarem
#'
EM_templateICA.spatial <- function(
  template_mean, template_var, meshes, BOLD,
  theta0, C_diag, maxiter=100,  usePar=FALSE, epsilon=0.01, verbose=FALSE){

  INLA_check()

  if(!all.equal(dim(template_var), dim(template_mean))) stop('The dimensions of template_mean and template_var must match.')

  ntime <- ncol(BOLD) #length of timeseries
  nvox <- nrow(BOLD) #number of data locations
  if(ntime > nvox) warning('More time points than data locations. Are you sure the data is oriented properly?')
  if(nrow(template_mean) != nvox) stop('Templates and BOLD must have the same number of data locations (columns).')

  Q <- ncol(template_mean) #number of ICs
  if(Q > nvox) stop('Cannot estimate more ICs than data locations.')
  if(Q > ntime) stop('Cannot estimate more ICs than time points.')

  template_var[template_var < 1e-6] <- 1e-6 #to prevent problems when inverting covariance

  #pre-compute s0, D and D^{-1}*s0
  s0_vec <- as.vector(template_mean) #grouped by IC
  D_vec <- as.vector(sqrt(template_var)) #grouped by IC
  D <- Diagonal(nvox*Q, D_vec)
  Dinv_s0 <- INLA::inla.qsolve(Q = D, B=matrix(s0_vec, ncol=1), method='solve')

  # ### REFINE STARTING VALUE FOR KAPPA
  #
  # if(verbose) cat('Refining starting value for kappa \n')
  #
  # # Determine direction of change:
  # # Positive change --> search for kappa_max, set kappa_min to kappa1.
  # # Negative change --> search for kappa_min, set kappa_max to kappa1.
  # kappa_min <- kappa_max <- theta0$kappa[1]
  # theta1 <- UpdateTheta_templateICA.spatial(template_mean, template_var, meshes, BOLD, theta0, C_diag, s0_vec, D, Dinv_s0, common_smoothness=TRUE, verbose=verbose, update='kappa')
  # kappa_diff0 <- theta1$kappa[1] - theta0$kappa[1]
  # theta <- theta0
  #
  # kappa_diff <- kappa_diff0
  # if(kappa_diff0 < 0){
  #
  #   if(verbose) cat('...Kappa decreasing, finding lower bound for kappa search \n ')
  #
  #   kappa_min <- kappa_min/2
  #   while(kappa_diff < 0){
  #     if(verbose) cat(paste0('... testing kappa = ',round(kappa_min,3),'\n '))
  #     theta$kappa <- rep(kappa_min, Q)
  #     theta1 <- UpdateTheta_templateICA.spatial(template_mean, template_var, meshes, BOLD, theta, C_diag, s0_vec, D, Dinv_s0, common_smoothness=TRUE, verbose=verbose, update='kappa')
  #     kappa_diff <- theta1$kappa[1] - theta$kappa[1]
  #     if(kappa_diff > 0) {
  #       #set minimum and stop here
  #       kappa_min <- theta1$kappa[1]
  #       break
  #     } else {
  #       #set new value for kappa
  #       kappa_min <- kappa_min/2
  #     }
  #   }
  # } else if(kappa_diff0 > 0){
  #
  #   if(verbose) cat('...Kappa increasing, finding upper bound for kappa search \n ')
  #
  #   kappa_max <- kappa_max*2
  #   while(kappa_diff > 0){
  #     if(verbose) cat(paste0('... testing kappa = ',round(kappa_max, 3),'\n '))
  #     theta$kappa <- rep(kappa_max, Q)
  #     theta1 <-  UpdateTheta_templateICA.spatial(template_mean, template_var, meshes, BOLD, theta, C_diag, s0_vec, D, Dinv_s0, common_smoothness=TRUE, verbose=FALSE, update='kappa')
  #     kappa_diff <- theta1$kappa[1] - theta$kappa[1]
  #     if(kappa_diff < 0) {
  #       #set maximum and stop here
  #       kappa_max <- theta1$kappa[1]
  #       break
  #     } else {
  #       #set new value for kappa
  #       kappa_max <- kappa_max*2
  #     }
  #   }
  # }
  #
  # #use binary search until convergence
  # if(verbose) cat('...Starting binary search for starting value of kappa \n ')
  # kappa_test <- (kappa_min + kappa_max)/2
  # kappa_change <- 1
  # while(kappa_change > epsilon){
  #   if(verbose) cat(paste0('... testing kappa = ',round(kappa_test, 3),'\n '))
  #   theta$kappa <- rep(kappa_test, Q)
  #   theta1 <-  UpdateTheta_templateICA.spatial(template_mean, template_var, meshes, BOLD, theta, C_diag, s0_vec, D, Dinv_s0, common_smoothness=TRUE, verbose=FALSE, update='kappa')
  #   kappa_diff <- theta1$kappa[1] - theta$kappa[1] #which direction is the estimate of kappa moving in?
  #   if(kappa_diff > 0) {
  #     kappa_min <- theta1$kappa[1]  #reset minimum to current value
  #     kappa_test <- (theta1$kappa[1] + kappa_max)/2 #go halfway to max
  #   } else {
  #     kappa_max <- theta1$kappa[1]  #reset maximum to current value
  #     kappa_test <- (theta1$kappa[1] + kappa_min)/2 #go halfway to min
  #   }
  #
  #   #how much different is the next value of kappa to be tested versus the current one?
  #   kappa_change <- abs((kappa_test - theta1$kappa[1])/theta1$kappa[1])
  # }


  ### RUN SQUAREM ALGORITHM UNTIL CONVERGENCE

  #theta0 <- theta1 #last tested value of kappa0
  theta0$LL <- c(0,0) #log likelihood
  theta0_vec <- unlist(theta0[1:2]) #everything but LL
  names(theta0_vec)[1] <- 0 #store LL value in names of theta0_vec (required for squarem)

  t00000 <- Sys.time()
  saveRDS(list(
    par=theta0_vec, fixptfn = UpdateThetaSQUAREM_templateICA, objfn=LL_SQUAREM,
    control=list(trace=verbose, intermed=TRUE, tol=epsilon, maxiter=maxiter),
    tmean=template_mean, tvar=template_var, meshes=meshes,
    BOLD=BOLD, C_diag=C_diag, s0_vec=s0_vec, D=D, Dinv_s0=Dinv_s0, verbose=TRUE
  ), "tICA_spatial_pre_squarem1")
  result_squarem <- squarem(
    par=theta0_vec, fixptfn = UpdateThetaSQUAREM_templateICA, objfn=LL_SQUAREM,
    control=list(trace=verbose, intermed=TRUE, tol=epsilon, maxiter=maxiter),
    template_mean, template_var, meshes,
    BOLD, C_diag, s0_vec, D, Dinv_s0, verbose=TRUE
  )
  if(verbose) print(Sys.time() - t00000)

  path_A <- result_squarem$p.inter[,1:(Q^2)]
  path_kappa <- result_squarem$p.inter[,(Q^2)+(1:Q)]
  path_LL <- result_squarem$p.inter[,ncol(result_squarem$p.inter)]
  theta_path <- list(A=path_A, kappa=path_kappa, LL=path_LL)

  theta_MLE <- theta0
  theta_MLE$A <- matrix(result_squarem$par[1:(Q^2)], Q, Q)
  theta_MLE$kappa <- result_squarem$par[(Q^2)+(1:Q)]
  theta_MLE$LL <- as.numeric(names(result_squarem$par)[1])

  #success <- (result_squarem$convergence==0) #0 indicates convergence, 1 indicates failure to converge within maxiter
  numiter <- result_squarem$fpevals #number of parameter update steps (approximately 3x the number of SQUAREM iterations)

  ### Compute final posterior mean of subject ICs
  if(verbose) cat('Computing final posterior mean of subject ICs \n')
  mu_cov_s <- UpdateTheta_templateICA.spatial(template_mean,
                                             template_var,
                                             meshes,
                                             BOLD,
                                             theta_MLE[1:2],
                                             C_diag,
                                             s0_vec,
                                             D,
                                             Dinv_s0,
                                             #common_smoothness=common_smoothness,
                                             verbose=verbose,
                                             return_MAP=TRUE)

  result <- list(subjICmean=mu_cov_s$mu_s,
                 subjICcov=mu_cov_s$cov_s,
                 Omega = mu_cov_s$Omega_s,
                 theta_MLE=theta_MLE,
                 theta_path=theta_path,
                 numiter=numiter,
                 squarem = result_squarem,
                 template_mean = template_mean,
                 template_var=template_var)
  return(result)
}

#' @rdname EM_templateICA
EM_templateICA.independent <- function(
  template_mean, template_var, BOLD, theta0, C_diag, maxiter=100, epsilon=0.01, usePar=FALSE, verbose){

  if(!all.equal(dim(template_var), dim(template_mean))) stop('The dimensions of template_mean and template_var must match.')

  ntime <- ncol(BOLD) #length of timeseries
  nvox <- nrow(BOLD) #number of brain locations
  if(ntime > nvox) warning('More time points than brain locations. Are you sure?')
  if(nrow(template_mean) != nvox) stop('Templates and BOLD must have the same number of brain locations (columns).')

  Q <- ncol(template_mean) #number of ICs
  if(Q > nvox) stop('Cannot estimate more ICs than brain locations.')
  if(Q > ntime) stop('Cannot estimate more ICs than time points.')

  iter <- 1
  theta <- theta0
  success <- 1
  template_var[template_var < 1e-6] <- 1e-6 #to prevent problems when inverting covariance

  err <- 1000 #large initial value for difference between iterations
  while(err > epsilon){

    if(verbose) cat(paste0(' ~~~~~~~~~~~~~~~~~~~~~ ITERATION ', iter, ' ~~~~~~~~~~~~~~~~~~~~~ \n'))

    t00 <- Sys.time()
    theta_new = UpdateTheta_templateICA.independent(template_mean, template_var, BOLD, theta, C_diag, verbose=verbose)
    if(verbose) print(Sys.time() - t00)

    ### Compute change in parameters

    # A_old <- theta$A
    # A_new <- theta_new$A
    #2-norm <- largest eigenvalue <- sqrt of largest eigenvalue of AA'
    A_change <- norm(as.vector(theta_new$A - theta$A), type="2") / norm(as.vector(theta$A), type="2")

    # nu0_sq_old <- theta$nu0_sq
    # nu0_sq_new <- theta_new$nu0_sq
    nu0_sq_change <- abs(theta_new$nu0_sq - theta$nu0_sq)/theta$nu0_sq

    change <- c(A_change, nu0_sq_change)
    err <- max(change)
    change <- format(change, digits=3, nsmall=3)
    if(verbose) cat(paste0('Iteration ',iter, ': Difference is ',change[1],' for A, ',change[2],' for nu0_sq \n'))

    ### Move to next iteration
    theta <- theta_new
    iter <- iter + 1
    if(iter > maxiter){
      success <- 0
      warning(paste0('Failed to converge within ', maxiter,' iterations'))
      break() #exit loop
    }
  }

  ### Compute final posterior mean of subject ICs

  #A = theta$A
  At_nu0Cinv <- crossprod(theta$A, diag(1/(C_diag*theta$nu0_sq)))
  At_nu0Cinv_A <- At_nu0Cinv %*% theta$A

  if (usePar) {
    if (!requireNamespace("foreach", quietly = TRUE)) {
      stop(
        "Package \"foreach\" needed to parallel loop over voxels. Please install it.",
        call. = FALSE
      )
    }

    `%dopar%` <- foreach::`%dopar%`
    q <- foreach::foreach(v = seq(nvox), .combine=rbind) %dopar% {
      y_v <- BOLD[v,]
      s0_v <- template_mean[v,]
      E_v_inv <- diag(1/template_var[v,])
      Sigma_s_v <- solve(E_v_inv + At_nu0Cinv_A)
      c(Sigma_s_v	%*% (At_nu0Cinv %*% y_v + E_v_inv %*% s0_v), diag(Sigma_s_v))
    }
    miu_s <- q[,seq(Q)]
    var_s <- q[,seq(Q+1, Q*2)]
  } else {
    miu_s <- matrix(NA, nrow=nvox, ncol=Q)
    var_s <- matrix(NA, nrow=nvox, ncol=Q)
    for (v in seq(nvox)) {
      y_v <- BOLD[v,]
      s0_v <- template_mean[v,]
      E_v_inv <- diag(1/template_var[v,])
      Sigma_s_v <- solve(E_v_inv + At_nu0Cinv_A)
      miu_s[v,] <- Sigma_s_v	%*% (At_nu0Cinv %*% y_v + E_v_inv %*% s0_v) #Qx1
      var_s[v,] <- diag(Sigma_s_v)
    }
  }

  result <- list(subjICmean=miu_s,
                 subjICse=sqrt(var_s),
                 theta_MLE=theta,
                 success_flag=success,
                 error=err,
                 numiter=iter-1,
                 template_mean = template_mean,
                 template_var = template_var)
  return(result)
}


#' @name UpdateTheta_templateICA
#' @rdname UpdateTheta_templateICA
#'
#' @title Parameter Estimates in EM Algorithm for Template ICA Model
#'
#' @param template_mean (\eqn{V \times Q} matrix) mean maps for each IC in template
#' @param template_var (\eqn{V \times Q} matrix) between-subject variance maps for each IC in template
#' @param meshes \code{NULL} for spatial independence model, otherwise a list of
#'  objects of class "templateICA_mesh" containing the triangular mesh (see
#'  \code{\link{make_mesh}}) for each brain structure.
#' @param BOLD  (\eqn{V \times Q} matrix) dimension-reduced fMRI data
#' @param theta (list) current parameter estimates
#' @param C_diag \eqn{(Qx1)} diagonal elements of residual covariance after dimension reduction
#' @param s0_vec Vectorized template means
#' @param D Sparse diagonal matrix of template standard deviations
#' @param Dinv_s0 The inverse of D times s0_vec
# @param common_smoothness If \code{TRUE}. use the common smoothness version of the spatial template ICA model, which assumes that all IC's have the same smoothness parameter, \eqn{\kappa}
#' @param verbose If \code{TRUE}, display progress of algorithm. Default: \code{FALSE}.
#' @param return_MAP If \code{TRUE}. return the posterior mean and precision of
#'  the latent fields instead of the parameter estimates. Default: \code{FALSE}.
#' @param update Which parameters to update. Either \code{"all"}, \code{"A"} or \code{"kappa"}.
#'
#' @return An updated list of parameter estimates, theta, OR if return_MAP=TRUE, the posterior mean and precision of the latent fields
#'
NULL

#' @rdname UpdateTheta_templateICA
#'
#' @importFrom stats optimize
# @importFrom INLA inla.qsolve inla.qinv inla.setOption
#' @importFrom Matrix Matrix sparseMatrix
UpdateTheta_templateICA.spatial <- function(template_mean, template_var, meshes, BOLD, theta, C_diag, s0_vec, D, Dinv_s0, verbose=FALSE, return_MAP=FALSE, update=c('all','kappa','A')){

  INLA_check()

  Q <- ncol(template_mean)
  nvox <- nrow(BOLD)
  ntime <- ncol(BOLD)
  #spde <- mesh$spde

  #initialize new parameter values
  theta_new <- theta

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

  if(verbose) cat('Computing Posterior Moments of S \n')

  y_vec = as.vector(t(BOLD)) #grouped by locations

  #if(verbose) cat('...posterior precision \n') # less than 1 sec
  #Compute SPDE matrices (F, G, GFinvG) and Sigma_inv (QVxQV), a sparse block diagonal matrix
  stuff <- compute_R_inv(meshes, kappa=theta$kappa, C1=1/(4*pi))
  R_inv <- bdiag(rep(list(bdiag(stuff$Rinv)), Q)) #block diagonal over components
  Fmat <- bdiag(rep(list(bdiag(stuff$Fmat)), Q)) #block diagonal over components
  Gmat <- bdiag(rep(list(bdiag(stuff$Gmat)), Q)) #block diagonal over components
  GFinvG <- bdiag(rep(list(stuff$GFinvG), Q))  #block diagonal over components
  Amat <- stuff$Amat

  #set up P as a sparse matrix (see OneNote for illustration of this)
  P <- make_Pmat(Q, nvox)

  #1. Compute mu_s
  #if(verbose) cat('...posterior mean \n') #20 seconds with pardiso! (1 hour without)
  stuff <- compute_mu_s(y_vec, D, Dinv_s0, R_inv, theta, P, C_diag)
  mu_s <- stuff$mu
  #m_vec <- stuff$m
  Omega <- stuff$Omega
  Omega_inv_m <- stuff$Omega_inv_m

  #if(verbose) print(summary(as.vector(mu_s)))

  if(return_MAP){
    cov_s = (D %*% INLA::inla.qinv(Omega) %*% D) #only computes diagonal elements and elements that are non-zero in precision of s (corresponding to neighboring locations in mesh)
    return(list(mu_s = mu_s, cov_s = cov_s, Omega_s = Omega))
    stop()
  }

  if(update_A){
    if(verbose) cat('Updating A \n')
    #t00 <- Sys.time()

    #Compute first matrix appearing in A-hat
    Pmu_vec <- P %*% mu_s
    yPmu <- matrix(0, nrow=ntime, ncol=Q)
    for(v in 1:nvox){
      inds_y <- (1:ntime) + (v-1)*ntime
      inds_Pmu <- (1:Q) + (v-1)*Q
      y_v <- y_vec[inds_y]
      Pmu_v <- Pmu_vec[inds_Pmu]
      yPmu_v <- outer(y_v, Pmu_v)
      yPmu <- yPmu + yPmu_v #sum up to later divide by nvox to average
    }

    #Compute second matrix appearing in A-hat
    Omega_PP <- P %*% Omega %*% t(P)
    D_PP <- P %*% D %*% t(P)

    #if(verbose) cat('..Calculating only non-sparse covariance terms \n')

    #set up sparse indicator matrix for needed elements of Omega_PP^{-1}
    oneblock <- matrix(1, nrow=Q, ncol=Q)
    ind_mat <- bdiag_m(rep(list(oneblock), nvox))

    #augment Omega_PP with additional non-sparse locations
    attributes(ind_mat)$x <- 0*attributes(ind_mat)$x
    Omega_PP_aug <- Omega_PP + ind_mat
    Omega_PP_inv <- INLA::inla.qinv(Omega_PP_aug)

    T_mat <- matrix(0, Q, Q)
    for(v in 1:nvox){
      inds <- (1:Q) + (v-1)*Q

      #T1 matrix
      Omega_PP_inv_vv <- Omega_PP_inv[inds,inds]
      T1_vv <- D_PP[inds,inds] %*% Omega_PP_inv_vv %*% D_PP[inds,inds]

      #T2 matrix
      Pmu_v <- Pmu_vec[inds]
      T2_vv <- outer(Pmu_v, Pmu_v)

      #T(v,v)
      T_vv <- T1_vv + T2_vv
      T_mat <- T_mat + T_vv
    }

    A_hat <- t(solve(t(T_mat), t(yPmu))) # = yPmu %*% T_mat_inv
    theta_new$A <- A_hat

    #first part of Q1: sum_v 2 y(v)' C^{-1} A t(v) = sum_v Trace{ 2 C^{-1} A t(v) y(v)' } = Trace{ 2 C^{-1} A sum_v [ t(v) y(v)' ] }, where sum_v [ t(v) y(v)' ] = A_part1'
    LL1_part1 <- sum(diag( 2 * diag((1/C_diag)) %*% A_hat %*% t(yPmu) ))
    # LL1_part1 <- 0
    # for(v in 1:nvox){
    #   inds_y <- (1:ntime) + (v-1)*ntime
    #   inds_Pmu <- (1:Q) + (v-1)*Q
    #   y_v <- y_vec[inds_y]
    #   y_v_Cinv <- y_v * (1/C_diag)
    #   Pmu_v <- Pmu_vec[inds_Pmu] #t(v) vector in paper
    #
    #   LL1_v <- 2 * t(y_v_Cinv) %*% A_hat %*% Pmu_v
    #   LL1_part1 <- LL1_part1 + as.numeric(LL1_v)
    # }

    #
    #second part of Q1: sum_v Trace{ A' C^{-1} A T(v,v) } = Trace{ A' C^{-1} A sum_v T(v,v) }, where sum_v T(v,v) is A_part2
    LL1_part2 <- sum(diag( t(A_hat) %*% diag((1/C_diag)) %*% A_hat %*% T_mat ))
    #
    theta_new$LL[1] <- LL1_part1 - LL1_part2
    #print(LL1_part1)
    #print(LL1_part2*2)

  }


  #keep value of nu0sq_hat from PCA-based dimension reduction
  nu0sq_hat <- theta$nu0_sq
  theta_new$nu0_sq <- nu0sq_hat[1]

  ##########################################
  ### E-STEP for kappa
  ##########################################

  if(update_kappa){

    if(verbose) cat('Updating kappa  \n')

    # Likelihood in terms of kappa_q's.
    # LL2_part1 = log(det(R_q_inv))
    # LL2_part2 = Tr(R_q_inv * Omega_inv_q) + Tr(R_q_inv * W_hat_q)
    # LL2_part3 = u_q' * R_q_inv * v_hat_q
    #             u_q = Dinv * s0
    #             v_q = 2 Omega_inv * m - Dinv * s0

    # Tr(R_q_inv * Omega_inv_q) --> Use INLA::inla.inv to compute necessary elements (non-zeros in R_q_inv) of Omega_inv_q for q=1,...,Q
    # Tr(R_q_inv * W_hat_q) --> W_hat_q = outer(Omega_inv*m,Omega_inv*m)_qq, where Omega_inv*m is known. Just calculate necessary elements (non-zeros in R_q_inv) of W_hat_q for q=1,...,Q

    if(verbose) cat('..Computing trace terms in log-likelihood of kappa \n')

    ### SET UP FOR PART 2

    #1. Generate Monte Carlo samples for estimation of Omega_hat^{-1}_qq

    # if(verbose) cat(paste0('....drawing ',nsamp,' Monte Carlo samples for covariance estimation \n'))
    #
    # nsamp <- 5000 #number of samples
    # musamp <- INLA::inla.qsample(nsamp, Q = Omega, b=rep(0, V*Q), mu=rep(0, V*Q), num.threads=8) #sample from N(0, Omega^(-1)) to estimate Omega^(-1) (diagonal blocks only)
    #9 sec for 500 samples with V*Q = 2530*3 = 12,642

    #2. Determine non-zero terms of R_q_inv

    #if(verbose) cat('....calculating only non-sparse covariance terms \n')

    #set up estimation of only terms necessary for trace computation
    nonzero_Rinv <- as.matrix(1*(R_inv[1:nvox,1:nvox] != 0))
    #diag(nonzero_Rinv) <- 1 #make sure all diagonal elements estimated (required for Trace 1)
    nonzero_cols <- which(R_inv[1:nvox,1:nvox] != 0, arr.ind = TRUE)[,2] #column indices of non-zero locations

    #augment Omega with additional non-sparse locations
    ind_mat <- 1*(R_inv[1:nvox,1:nvox] != 0)
    attributes(ind_mat)$x <- 0*attributes(ind_mat)$x
    #Omega_aug <- Omega + bdiag(rep(list(ind_mat), Q)) #augmentation not needed here

    #compute elements of Omega_inv that are non-sparse in Omega and in R^{-1}
    Omega_inv <- INLA::inla.qinv(Omega)  #don't need any additional elements since Omega includes all non-zero elements of R_inv
    #~1 min for V=10000, L=16

    #if(verbose) cat('....setting up helper objects for trace computation \n') #15 sec (Q=16, nvox=5500)

    # #3. Set up matrices needed for computation of only necessary elements of Omega_hat^{-1}
    #
    # X <- matrix(1, nrow=nsamp, ncol=nvox) #placeholder for X in Omega^(-1) = XX'
    # bigX_left <- KhatriRao(diag(1, nvox), X) # -- 13 SEC (do one time instead of KhatriRao(diag(1, nvox), Xctr_q) for each q)
    # bigX_right <- KhatriRao(nonzero_Rinv, X) # -- 160 SEC (do one time instead of KhatriRao(as.matrix(1*(K_q != 0)), Xctr_q) for each q)
    # inds_left_X <- which(bigX_left != 0, arr.ind=TRUE)
    # inds_right_X <- which(bigX_right != 0, arr.ind=TRUE)
    #
    #4. Set up vectors needed for computation of only necessary elements of W_hat

    #x <- matrix(1, nrow=1, ncol=nvox) #placeholder for x in M = xx'
    #bigx_left <- diag(1, nvox) #KhatriRao(diag(1, nvox), x)
    #bigx_right <-  KhatriRao(nonzero_Rinv, x))
    #inds_left_x <- which(diag(1, nvox) != 0, arr.ind=TRUE) #1,1 2,2 3,3 etc. (bigx_left = I_V) <- just used to make a diagonal matrix below, don't actually need this
    inds_right_x <- which(nonzero_Rinv != 0, arr.ind=TRUE) #which(bigx_right != 0, arr.ind=TRUE)

    #if(verbose) cat('....computing necessary elements of RHS matrices in trace terms \n')
    #15 sec for Q=16, nvox=5500

    #if(!common_smoothness) OplusW <- vector('list', length=Q) #Omega_inv_qq + W_hat_qq

    for(q in 1:Q){
      #if(verbose) cat(paste('......block ',q,' of ',Q,' \n'))
      inds_q <- (1:nvox) + (q-1)*nvox

      # COMPUTE OMEGA_INV[q,q] (NECESSARY ENTRIES ONLY)

      # Xctr_q <- scale(t(musamp[inds_q,]), scale=FALSE)
      #
      # #left-multiplication matrix
      # vals_left_X <- as.vector(Xctr_q)
      # bigX_left_q <- sparseMatrix(i = inds_left_X[,1], j = inds_left_X[,2], x = vals_left_X)
      # #right-multiplication matrix
      # X_repcols <- Xctr_q[,nonzero_cols]
      # vals_right_X <- as.vector(X_repcols)
      # bigX_right_q <- sparseMatrix(i = inds_right_X[,1], j = inds_right_X[,2], x = vals_right_X)
      # #multiply together
      # Omega_inv_qq <- t(bigX_left_q) %*% bigX_right_q / (nsamp - 1)

      Omega_inv_qq <- Omega_inv[inds_q,inds_q]

      # COMPUTE W_hat[q,q] = Omega_inv_m[q] * Omega_inv_m[q' (NECESSARY ENTRIES ONLY)

      # T_q = matrix that selects the qth block of size nvox from a matrix with nvox*Q rows
      e_q <- Matrix(0, nrow=1, ncol=Q, sparse=TRUE)
      e_q[1,q] <- 1
      T_q <- kronecker(e_q, Diagonal(nvox))
      Omega_inv_m_q <- T_q %*% Omega_inv_m

      #left-multiplication matrix
      bigx_left_q <- Diagonal(nvox, x = as.vector(Omega_inv_m_q)) #sparseMatrix(i = inds_left_x[,1], j = inds_left_x[,2], x = as.vector(Omega_inv_m_ell_g))
      #right-multiplication matrix
      x_repcols <- t(Omega_inv_m_q)[,nonzero_cols] #repeat elements of Omega_inv_m_q to fill in all non-zero cells of R^{-1}
      bigx_right_q <- sparseMatrix(i = inds_right_x[,1], j = inds_right_x[,2], x = as.vector(x_repcols))
      #multiply together
      W_hat_qq <- t(bigx_left_q) %*% bigx_right_q

      # COMBINE OMEGA_INV[q,q] AND W_hat[q,q] --> two trace terms involving R_q_inv can be combined

      OplusW_qq <- Omega_inv_qq + W_hat_qq
      #if(common_smoothness){
        if(q==1) OplusW <- OplusW_qq else OplusW <- OplusW + OplusW_qq
      #} else {
      #  OplusW[[q]] <- OplusW_qq
      #}
    }


    ### SET UP FOR PART 3

    u <- Dinv_s0
    v <- 2*Omega_inv_m - u

    # NUMERICALLY SOLVE FOR MLE OF KAPPA

    if(verbose) cat("..Performing numerical optimization \n")
    #~8 seconds for nvox=5500, Q=3 (comon smoothness)
    #~30 seconds for nvox=10000, Q=16

    # nmeshes <- length(meshes)
    # Amat_list <- vector('list', length=nmeshes)
    # for(m in 1:nmeshes){
    #   Amat_list[[m]] <- meshes[[m]]$A
    # }
    # Amat <- bdiag(Amat_list)

    # kappa_vals <- seq(0.2,0.7,0.05)
    # LL_vals <- matrix(NA, 3, length(kappa_vals))
    # for(ii in 1:length(kappa_vals)){
    #   LL_vals[,ii] <- LL2_kappa(kappa_vals[ii], Amat=Amat, Fmat=Fmat, Gmat=Gmat, GFinvG=GFinvG, OplusW=OplusW, u=u, v=v, Q=Q)
    # }

    #if(common_smoothness){
      kappa_opt <- optimize(LL2_kappa, lower=0, upper=5, maximum=TRUE,
                            Amat=Amat, Fmat=Fmat, Gmat=Gmat, GFinvG=GFinvG, OplusW=OplusW, u=u, v=v, Q=Q) #Q=Q to indicate common smoothness model to the LL2_kappa function
      LL2 <- kappa_opt$objective
      kappa_opt <- rep(kappa_opt$maximum, Q)
    # } else {
    #   kappa_opt <- LL2 <- rep(NA, Q)
    #   for(q in 1:Q){
    #     if(verbose) cat(paste('Optimization ',q,' of ',Q,' \n'))
    #     inds_q <- (1:nvox) + (q-1)*nvox
    #     kappa_opt_q <- optimize(templateICAr::LL2_kappa, lower=0, upper=5, maximum=TRUE,
    #                             Amat=Amat, Fmat=Fmat, Gmat=Gmat, GFinvG=GFinvG, OplusW=OplusW[[q]], u=u[inds_q], v=v[inds_q], Q=NULL)
    #     LL2[q] <- kappa_opt_q$objective
    #     kappa_opt[q] <- (kappa_opt_q$maximum)
    #   }
    #   LL2 <- sum(LL2)
    # }

    theta_new$kappa <- kappa_opt
    theta_new$LL[2] <- LL2

  }

  # RETURN NEW PARAMETER ESTIMATES

  return(theta_new)

}

#' @rdname UpdateTheta_templateICA
UpdateTheta_templateICA.independent <- function(template_mean, template_var, BOLD, theta, C_diag, verbose){

  nQ <- ncol(template_mean)
  nV <- nrow(BOLD)
  nT <- ncol(BOLD)

  #initialize new objects
  theta_new <- list(A = matrix(NA, nT, nQ), nu0_sq = NA)
  #two parts of product for A-hat (construct each looping over voxels)
  A_part1 <- matrix(0, nT, nQ)
  A_part2 <- matrix(0, nQ, nQ)

  A <- theta$A # TxQ
  nu0_sq <- theta$nu0_sq # const
  nu0C_inv <- diag(1/(C_diag*nu0_sq)) #Sigma0_inv in matlab code # TxT
  At_nu0Cinv <- crossprod(A, nu0C_inv) # QxT
  At_nu0Cinv_A <- At_nu0Cinv %*% A # QxQ

  if(verbose) cat('Updating A \n')

  #store posterior moments for M-step of nu0_sq
  # miu_s <- matrix(NA, nrow=nV, ncol=nQ) # not used anymore
  # miu_ssT <- array(NA, dim=c(nV, nQ, nQ)) # not used anymore

  for (vv in 1:nV) {
    y_v <- BOLD[vv,] # T
    s0_v <- template_mean[vv,] # Q

    ##########################################
    ### E-STEP FOR A AND nu0^2: POSTERIOR MOMENTS OF s_i(v)
    ##########################################

    E_v_inv <- diag(1/template_var[vv,]) # QxQ
    Sigma_s_v <- solve(E_v_inv + At_nu0Cinv_A) # QxQ
    miu_s_v <- Sigma_s_v	%*% (At_nu0Cinv %*% y_v + E_v_inv %*% s0_v) #Qx1
    miu_ssT_v <- tcrossprod(miu_s_v) + Sigma_s_v #QxQ
    # miu_s[vv,] <- miu_s_v #save for M-step of nu0_sq
    # miu_ssT[vv,,] <- miu_ssT_v #save for M-step of nu0_sq

    ##########################################
    ### M-STEP FOR A: CONSTRUCT PARAMETER ESTIMATES
    ##########################################

    A_part1 <- A_part1 + tcrossprod(as.matrix(y_v), miu_s_v) #TxQ + (Tx1)*(1xQ)
    A_part2 <- A_part2 + miu_ssT_v #QxQ
  }

  #A_hat <- orthonorm(A_part1 %*% solve(A_part2))
  A_hat <- (A_part1 %*% solve(A_part2))

  ##########################################
  ### M-STEP FOR nu0^2: CONSTRUCT PARAMETER ESTIMATES
  ##########################################

  # cat('Updating Error Variance nu0_sq \n')
  #
  # #use A-hat or A?
  #
  # Cinv <- diag(1/C_diag)
  # Cinv_A <- Cinv %*% A_hat
  # At_Cinv_A <- t(A_hat) %*% Cinv %*% A_hat
  # nu0sq_part1 <- nu0sq_part2 <- nu0sq_part3 <- 0
  #
  # for(vv in 1:nV){
  #
  #   y_v <- BOLD[,vv]
  #   nu0sq_part1 <- nu0sq_part1 + t(y_v) %*% Cinv %*% y_v
  #   nu0sq_part2 <- nu0sq_part2 + t(y_v) %*% Cinv_A %*% miu_s[,vv]
  #   nu0sq_part3 <- nu0sq_part3 + sum(diag(At_Cinv_A %*% miu_ssT[,,vv]))
  # }
  #
  # nu0sq_hat <- 1/(nQ*nV)*(nu0sq_part1 - 2*nu0sq_part2 + nu0sq_part3)
  nu0sq_hat <- theta$nu0_sq


  # RETURN NEW PARAMETER ESTIMATES
  theta_new$A <- A_hat
  theta_new$nu0_sq <- nu0sq_hat[1]
  return(theta_new)
}



#' Compute posterior mean and precision of s
#'
#' Compute posterior mean and precision matrix of s
#'
#' @param y_vec Vectorized, dimension-reduced fMRI data, grouped by locations. A vector of length \eqn{QV}.
#' @param D Sparse diagonal matrix of template standard deviations
#' @param Dinv_s0 The inverse of D times s0_vec
#' @param R_inv Estimate of inverse spatial correlation matrix (sparse)
#' @param theta List of current parameter estimates
#' @param P Permutation matrix for regrouping by locations (instead of by ICs.)
#' @param C_diag Diagonals of residual covariance of the first level model. A vector of length \eqn{Q}.
#'
#' @importFrom Matrix Diagonal
# @importFrom INLA inla.qsolve
#'
#' @return A list containing the posterior mean \eqn{\mu} (mu) and precision
#'  \eqn{\Omega} (Omega) of s=(s1,...,sQ), along with the supporting vector m,
#'  where \eqn{\mu = \Omega^{-1}m}.
#'
#' @keywords internal
#'
compute_mu_s <- function(y_vec, D, Dinv_s0, R_inv, theta, P, C_diag){

  INLA_check()

  # ntime <- length(C_diag) # not used
  Q <- ncol(theta$A)
  nvox <- nrow(P)/Q

  A <- theta$A
  if('nu0_sq' %in% names(theta)) nu0_sq <- theta$nu0_sq else nu0_sq <- 1

  #set up B, C, and products thereof
  ones <- Diagonal(nvox)
  B <- kronecker(ones, A)
  nu0C_inv <- kronecker(ones, diag(1/(C_diag*nu0_sq)))
  Pt_Bt_nu0C_inv <- t(P) %*% t(B) %*% nu0C_inv

  #compute m (using current parameter estimates in theta)
  m1_vec <- D %*% Pt_Bt_nu0C_inv %*% y_vec
  m2_vec <- R_inv %*% Dinv_s0
  m_vec <- m1_vec + m2_vec

  # Compute Omega (using current parameter estimates in theta)
  Omega <- R_inv + D %*% Pt_Bt_nu0C_inv %*% B %*% P %*% D

  # Compute mu_s|y by solving for x in the system of equations Omega*x = m
  Omega_inv_m <- INLA::inla.qsolve(Q = Omega, B=m_vec, method='solve') # <-- slowest part (up to 1 hour for Q=16, nvox=5500), but with inla.setOption(smtp='pardiso') goes down to 20 seconds!
  mu_s <- D %*% Omega_inv_m

  return(list(mu=mu_s, m=m_vec, Omega=Omega, Omega_inv_m=Omega_inv_m))
}


#' Compute SPDE and prior precision matrices for S
#'
#' Compute SPDE matrices (F, G and GFinvG) and prior precision matrix for S
#'
#' @param meshes \code{NULL} for spatial independence model, otherwise a list of
#'  objects of class "templateICA_mesh" containing the triangular mesh (see
#'  \code{\link{make_mesh}}) for each brain structure.
#' @param kappa Current estimates of SPDE parameter kappa for each latent field
#' @param C1 Constant, equal to \eqn{1/(4*pi)} for a 2-dimensional field with alpha=2
#' @param rm_extra If \code{TRUE}. remove extra (non-data) vertices from the mesh for greater computational efficiency
#'
#' @importFrom Matrix bdiag
#' @importFrom stats var
#'
#' @return A list containing R inverse and SPDE matrices
#'
#' @keywords internal
#'
compute_R_inv <- function(meshes, kappa, C1=1/(4*pi), rm_extra=FALSE){

  L <- length(kappa)
  if(length(kappa) > 1) { if(var(kappa)==0) onekappa <- TRUE }
  if(length(kappa)==1) onekappa <- TRUE
  if(onekappa) kappa <- kappa[1]

  if(!inherits(meshes, "list")) stop('meshes argument must be a list of objects of class "templateICA_mesh"')

  #SPDE matrices, needed to construct R_l_inv
  nmeshes <- length(meshes)
  Fmat_all <- Gmat_all <- GFinvG_all <- Rinv_all <- Amat_all <- vector('list', length=nmeshes)
  for(m in 1:nmeshes){
    mesh <- meshes[[m]]
    spde <- mesh$spde
    Fmat <- spde$param.inla$M0
    Gmat <- 1/2*(spde$param.inla$M1 + Matrix::t(spde$param.inla$M1))
    GFinvG <- spde$param.inla$M2 #this equals G %*% solve(F) %*% G

    #get data and nodata indices of mesh vertices
    Amat <- mesh$A #n_loc x n_mesh
    N <- ncol(mesh$A) #number of mesh locations
    inds_data <- which(colSums(Amat) > 0)
    inds_nodata <- setdiff(1:N, inds_data)

    #check if there are no non-data locations in mesh
    if(length(inds_nodata)==0) rm_extra <- TRUE #no need to marginalize out non-data locations since there are none

    #set up R^{-1} (VxV) as a sparse matrix

    #if(onekappa) #just compute Rinv once
      Qmat <- C1*(kappa^2 * Fmat + 2 * Gmat + kappa^(-2) * GFinvG)
      Q11 <- Qmat[inds_data,inds_data] # <- Amat %*% Qmat %*% t(Amat)
      if(rm_extra==FALSE){
        #if(verbose) cat('marginalizing out non-data locations in R\n')
        Q12 <- Qmat[inds_data, inds_nodata]
        Q21 <- Qmat[inds_nodata, inds_data]
        Q22 <- Qmat[inds_nodata,inds_nodata]
        Q22_inv <- solve(Q22)
        R_l_inv <- Q11 - (Q12 %*% Q22_inv %*% Q21)
      } else {
        R_l_inv <- Q11
      }
      Rinv <- R_l_inv
      #Rinv_list <- list(R_l_inv) #rep(list(R_l_inv), L)
    #} else { # compute Rinv block-wise
      ## TO DO: FIX THIS OR REMOVE NON-COMMON SMOOTHNESS OPTION
      # Rinv_list <- vector('list', L)
      # for(l in 1:L){
      #   kappa_l <- kappa[l]
      #   R_l_inv <- C1 * (kappa_l^2 * Fmat + 2 * Gmat + kappa_l^(-2) * GFinvG)
      #   Rinv_list[[l]] <- R_l_inv
      # }
    #}
    #Rinv <- bdiag(Rinv_list) #block diagonal over components

    Fmat_all[[m]] <- Fmat
    Gmat_all[[m]] <- Gmat
    GFinvG_all[[m]] <- GFinvG
    Rinv_all[[m]] <- Rinv
    Amat_all[[m]] <- Amat
  } #end loop over hemispheres

  # block diagonal over hemispheres
  return(list(Rinv=bdiag(Rinv_all),
              Fmat=bdiag(Fmat_all),
              Gmat=bdiag(Gmat_all),
              GFinvG=bdiag(GFinvG_all),
              Amat=bdiag(Amat_all)))
}

#' Make permutation matrix
#'
#' Create permutation matrix P to regroup elements of s by locations instead of by ICs
#'
#' @param Q The number of template ICs
#' @param nvox The number of spatial locations (\eqn{V})
#'
#' @details If s=(s1,...,sQ) is grouped by ICs 1,...Q, then Ps=(s(1),...,s(nvox)) is grouped by locations 1,...,nvox
#'
#' @importFrom Matrix sparseMatrix
#'
#' @return P Permutation matrix size \eqn{QVxQV}
#'
#' @keywords internal
#'
make_Pmat <- function(Q, nvox){
  cols <- 1:(Q*nvox)
  rows_P1 <- seq(1, (Q-1)*nvox+1, by=nvox)
  offset <- rep(0:(nvox-1), each=Q)
  rows <- rep(rows_P1, nvox) + offset
  t(sparseMatrix(i = rows, j = cols))
}


#' Block diagonal matrix
#'
#' Construct block diagonal matrix of many small blocks
#'
#' @description Code for function provided in examples of bdiag function from Matrix package
#'
#' @param lmat List of matrices
#'
#' @import Matrix
#' @importFrom methods new
#'
#' @return A sparse matrix obtained by combining the arguments into a block diagonal matrix
#'
#' @keywords internal
#'
bdiag_m <- function(lmat) {
  ## Copyright (C) 2016 Martin Maechler, ETH Zurich
  if(!length(lmat)) return(new("dgCMatrix"))
  stopifnot(is.list(lmat), is.matrix(lmat[[1]]),
            (k <- (d <- dim(lmat[[1]]))[1]) == d[2], # k x k
            all(vapply(lmat, dim, integer(2)) == k)) # all of them
  N <- length(lmat)
  if(N * k > .Machine$integer.max)
    stop("resulting matrix too large; would be  M x M, with M=", N*k)
  M <- as.integer(N * k)
  ## result: an   M x M  matrix
  new("dgCMatrix", Dim = c(M,M),
      ## 'i :' maybe there's a faster way (w/o matrix indexing), but elegant?
      i = as.vector(matrix(0L:(M-1L), nrow=k)[, rep(seq_len(N), each=k)]),
      p = k * 0L:M,
      x = as.double(unlist(lmat, recursive=FALSE, use.names=FALSE)))
}


#' Update theta SQUAREM
#'
#' Helper function for SQUAREM for estimating parameters
#'
#' @param theta_vec Vector of initial parameter values
#' @param template_mean Passed to UpdateTheta_templateICA function
#' @param template_var Passed to UpdateTheta_templateICA function
#' @param meshes Passed to UpdateTheta_templateICA function
#' @param BOLD Passed to UpdateTheta_templateICA function
#' @param C_diag Passed to UpdateTheta_templateICA function
#' @param s0_vec Passed to UpdateTheta_templateICA function
#' @param D Passed to UpdateTheta_templateICA function
#' @param Dinv_s0 Passed to UpdateTheta_templateICA function
# @param common_smoothness Passed to UpdateTheta_templateICA function
#' @param verbose Passed to UpdateTheta_templateICA function
#'
#' @return Vector of updated parameter values
#'
#' @keywords internal
#'
UpdateThetaSQUAREM_templateICA <- function(theta_vec, template_mean, template_var, meshes, BOLD, C_diag, s0_vec, D, Dinv_s0, verbose){

  ntime <- ncol(BOLD) #length of timeseries
  nvox <- nrow(BOLD) #number of data locations
  Q <- ncol(template_mean)

  #convert theta vector to list format
  A <- theta_vec[1:(Q*ntime)]
  #nu0_sq <- theta_vec[(Q^2)+1]
  kappa <- theta_vec[(Q*ntime)+(1:Q)]
  theta <- list(A = matrix(A, nrow=ntime, ncol=Q),
                #nu0_sq = nu0_sq,
                kappa = kappa)

  #update theta parameters
  if(verbose) cat('~~~~~~~~~~~ UPDATING PARAMETER ESTIMATES ~~~~~~~~~~~ \n')
  theta_new <- UpdateTheta_templateICA.spatial(
    template_mean,
    template_var,
    meshes,
    BOLD,
    theta,
    C_diag,
    s0_vec,
    D,
    Dinv_s0,
    #common_smoothness=common_smoothness,
    verbose=verbose
  )

  #convert theta_new list to vector format
  theta_new$A <- as.matrix(theta_new$A)
  theta_new_vec <- unlist(theta_new[1:2]) #everything but LL
  names(theta_new_vec)[1] <- sum(theta_new$LL)
  return(theta_new_vec)
}


#' Log-likelihood SQUAREM
#'
#' Helper function for SQUAREM for extracting log likelihood
#'
#' @param theta_vec Vector of current parameter values
#' @param template_mean Not used, but squarem will return error without
#' @param template_var  Not used, but squarem will return error without
#' @param meshes  Not used, but squarem will return error without
#' @param BOLD  Not used, but squarem will return error without
#' @param C_diag  Not used, but squarem will return error without
#' @param s0_vec  Not used, but squarem will return error without
#' @param D  Not used, but squarem will return error without
#' @param Dinv_s0  Not used, but squarem will return error without
# @param common_smoothness  Not used, but squarem will return error without
#' @param verbose  Not used, but squarem will return error without
#'
#' @return Negative log-likelihood given current values of parameters
#'
#' @keywords internal
#'
LL_SQUAREM <- function(theta_vec, template_mean, template_var, meshes, BOLD, C_diag, s0_vec, D, Dinv_s0, verbose){

  LL <- as.numeric(names(theta_vec)[1])
  return(-1*LL)

}


#' Compute part of appa log-likelihood
#'
#' Compute part of log-likelihood involving kappa (or kappa_q) for numerical optimization
#'
#' @param kappa Value of kappa for which to compute log-likelihood
#' @param Amat Mesh projection matrix
#' @param Fmat Matrix used in computation of SPDE precision
#' @param Gmat Matrix used in computation of SPDE precision
#' @param GFinvG Matrix used in computation of SPDE precision
#' @param OplusW Sparse matrix containing estimated values of RHS of trace in part 2 of log-likelihood. In common smoothness model, represents the sum over q=1,...,Q.
#' @param u Vector needed for part 3 of log-likelihood
#' @param v Vector needed for part 3 of log-likelihood
#' @param C1 For the unit variance case, \eqn{\tau^2 = C1/\kappa^2}, where \eqn{C1 = 1/(4\pi)} when \eqn{\alpha=2}, \eqn{\nu=1}, \eqn{d=2}
#' @param Q Equal to the number of ICs for the common smoothness model, or NULL for the IC-specific smoothness model
#'
#' @importFrom Matrix bdiag
#'
#' @return Value of log-likelihood at logkappa
#'
#' @keywords internal
#'
#' @details This is the function to be maximized in order to determine the MLE for \eqn{\kappa} or the \eqn{\kappa_q}'s in the M-step of the EM algorithm in spatial
#' template ICA.  In the model where \eqn{\kappa_q} can be different for each IC \eqn{q}, the optimization function factorizes over the \eqn{\kappa_q}'s.  This function computes
#' the value of the part of the optimization function pertaining to one of the \eqn{\kappa_q}'s.
#'
LL2_kappa <- function(kappa, Amat, Fmat, Gmat, GFinvG, OplusW, u, v, C1 = 1/(4*pi), Q=NULL){

  #get data and nodata indices of mesh vertices
  N <- ncol(Amat) #number of mesh locations
  inds_data <- which(colSums(Amat) > 0)
  inds_nodata <- setdiff(1:N, inds_data)

  #COMPUTE R_q_inv FOR CURRENT VALUE OF kappa

  Qmat <- C1*(kappa^2 * Fmat + 2 * Gmat + kappa^(-2) * GFinvG)
  Q11 <- Qmat[inds_data,inds_data] # <- Amat %*% Qmat %*% t(Amat)
  if(length(inds_nodata)>0){ #marginalize out non-data locations in mesh
    Q12 <- Qmat[inds_data, inds_nodata]
    Q21 <- Qmat[inds_nodata, inds_data]
    Q22 <- Qmat[inds_nodata,inds_nodata]
    Q22_inv <- solve(Q22)
    R_q_inv <- Q11 - (Q12 %*% Q22_inv %*% Q21)
  } else { #no non-data locations in mesh
    R_q_inv <- Q11
  }


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
    LL2_part3 <- t(u) %*% bdiag(rep(list(R_q_inv), Q)) %*% v #for common smoothness model
  } else {
    LL2_part3 <- t(u) %*% R_q_inv %*% v
  }

  result <- LL2_part1 - LL2_part2 + LL2_part3
  return(as.numeric(result))
}
