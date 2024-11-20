#' VB_FCtemplateICA
#'
#' VB Algorithm for FC Template ICA Model
#'
#' @param template_mean (\eqn{V \times Q} matrix) mean maps for each IC in the
#'  template, where \eqn{Q} is the number of ICs, and \eqn{V=nvox} is the number
#'  of data locations.
#' @param template_var  (\eqn{V \times Q} matrix) between-subject variance maps
#'  for each IC in the template.
#' @param template_FC (list) Parameters of functional connectivity template.
#' @param method_FC Variational Bayes (VB) method for FC template ICA model:
#'  \code{"VB1"} (default) uses a conjugate Inverse-Wishart prior for the cor(A);
#'  \code{"VB2"} draws samples from p(cor(A)) to emulate the population distribution
#'  using a combination of Cholesky, SVD, and random pivoting.
#' @param nsamp_u For VB1, the number of samples to generate from u ~ Gamma, where
#' A is Gaussian conditional on u. Default: \code{10000}.
#' @param CI_FC Level of posterior credible interval to construct for each FC element.
#' Default: \code{0.95}.
#' @param return_FC_samp Should the FC samples (\eqn{QxQxK}) be returned, where K
#' is the number of posterior samples generated? May be a large object.  For VB1,
#' K is equal to nsamp_u. For VB2, K is equal to the number of samples from p(G)
#' contained in template_FC.
#' @param prior_params Alpha and beta parameters of IG prior on \eqn{\tau^2}
#'  (error variance). Default: \code{0.001} for both.
#' @param BOLD (\eqn{V \times T} matrix) preprocessed fMRI data.
#' @param A0,S0,S0_var Initial guesses at latent variables: \code{A} (\eqn{TxQ}
#'  mixing matrix), \code{S} (\eqn{QxV} matrix of spatial ICs), and
#'  variance matrix \code{S0_var}.
#' @param maxiter Maximum number of VB iterations. Default: \code{100}.
#' @param miniter Minimum number of VB iterations. Default: \code{3}.
#' @param epsilon Smallest proportion change in parameter estimates between iterations.
#'  Default: \code{0.001}.
# @param eps_inter Intermediate values of epsilon at which to save results (used
#  to assess benefit of more stringent convergence rules). Default:
#  \code{NULL} (do not save). These values should be in decreasing order
#  (larger to smaller error) and all values should be between zero and \code{epsilon}.
#' @param usePar Parallelize the computation? Default: \code{FALSE}. Can be the
#' number of cores to use or \code{TRUE}, which will use the number available minus two.
#' @param PW Prewhiten to account for residual autocorrelation?
# @param tESS_correction Take into account effective sample size due to temporal autocorrelation?
# @param sESS_correction Take into account effective sample size due to spatial correlation?
#' @param verbose If \code{TRUE}, display progress of algorithm. Default: \code{FALSE}.
#'
#' @return A list of computed values, including the final parameter estimates.
#'
#' @importFrom fMRItools is_posNum
#' @importFrom abind abind
#' @importFrom stats acf
# @importFrom INLA inla.mesh.create
#' @keywords internal
#'
VB_FCtemplateICA <- function(
  template_mean, #VxQ
  template_var, #VxQ
  template_FC,
  method_FC = c('VB1','VB2'),
  nsamp_u = 10000,
  CI_FC = 0.95,
  return_FC_samp=FALSE,
  prior_params=c(0.001, 0.001),
  BOLD, #VxT
  TR=NULL,
  A0, S0, S0_var,
  maxiter=100,
  miniter=3,
  epsilon=0.001,
  #eps_inter=NULL,
  usePar=TRUE,
  PW=FALSE,
  #tESS_correction=FALSE,
  #sESS_correction=FALSE,
  verbose=FALSE){

  stopifnot(length(prior_params)==2)
  stopifnot(is.numeric(prior_params))
  stopifnot(is_posNum(maxiter, "numeric") && maxiter==round(maxiter))
  stopifnot(is_posNum(miniter, "numeric") && miniter==round(miniter))
  stopifnot(miniter <= maxiter)
  stopifnot(is_posNum(epsilon))
  # if (!is.null(eps_inter)) {
  #   stopifnot(is.numeric(eps_inter) && all(diff(eps_inter) < 0))
  #   stopifnot(eps_inter[length(eps_inter)]>0 && eps_inter[1]>epsilon)
  # }
  method_FC <- match.arg(method_FC, c("VB1", "VB2"))

  if (!all.equal(dim(template_var), dim(template_mean))) stop('The dimensions of template_mean and template_var must match.')
  ntime <- ncol(BOLD) #length of timeseries
  nvox <- nrow(BOLD) #number of brain locations
  if (ntime > nvox) warning('More time points than brain locations. Are you sure?')
  if (nrow(template_mean) != nvox) stop('Templates and BOLD must have the same number of brain locations (columns).')

  if(method_FC == 'VB1') nICs <- nrow(template_FC$psi)   #number of ICs
  if(method_FC == 'VB2') nICs <- nrow(template_FC$FC_samp_mean) #number of ICs
  if (ncol(template_mean) != nICs) stop('template_FC is incompatible with template_mean & template_var. The number of ICs in each must match.')
  if (nICs > nvox) stop('Cannot estimate more ICs than brain locations.')
  if (nICs > ntime) stop('Cannot estimate more ICs than time points.')

  iter <- 1
  success <- 1
  template_var[template_var < 1e-6] <- 1e-6 #to prevent problems when inverting covariance (this should rarely/never happen with NN variance)


  #0. (Optional) Compute temporal effective sample size (ESS) to correct sums over t

  #if(!is.numeric(TR)) TR <- 1 #assume 1-second temporal resolution if TR not available
  #nDCT <- fMRItools:::Hz2dct(ntime, TR, 0.05)
  #dct <- dct_bases(ntime, nDCT)
  #A00 <- nuisance_regression(A0, dct)

  #consider residual temporal dependence
  E0 <- BOLD - t(A0 %*% t(S0)) #VxT

  do_PW <- PW; rm(PW)
  if(do_PW){
    set.seed(1234)
    rvox <- sample(1:nvox, 100)
    #tESS <- rep(NA, 100)
    Cor_E <- matrix(0, ntime, ntime)
    for(k in 1:100){
      vox_k <- rvox[k]
      ar_order <- 6
      ar_k <- ar.yw(E0[k,], order.max = ar_order, aic = FALSE)$ar
      A_k <- toeplitz(c(1,-ar_k,rep(0, ntime - ar_order - 1)))
      A_k[upper.tri(A_k)] <- 0
      Cov_k <- solve(tcrossprod(A_k))
      Cor_k <- cov2cor(Cov_k)
      #acf_k <- acf(E0[k,], lag.max=ntime, plot=FALSE)$acf
      #Cor_k <- toeplitz(c(acf_k))
      #tESS[k] <- sum(diag(Cor_k))^2/sum(Cor_k^2) #Trace definition from Bretherton et al 1999 https://doi.org/10.1175/1520-0442(1999)012%3C1990:TENOSD%3E2.0.CO;2

      Cor_E <- Cor_E + Cor_k
    }
    Cor_E <- Cor_E/100
    svd_Cor_E <- svd(Cor_E, nv = FALSE)
    PW <- svd_Cor_E$u %*% diag(1/sqrt(svd_Cor_E$d)) %*% t(svd_Cor_E$u)

    # tESS <- mean(tESS)
    # if(verbose) cat(paste0('\t Temporal ESS accounting for autocorrelation: ', round(tESS), '\n'))
    # if(tESS >= ntime){
    #   tESS <- NULL
    #   tESS_correction <- FALSE
    #   if(verbose) cat(paste0('\t Skipping temporal ESS correction since estimated ESS >= ntime \n'))
    # }
  } else {
    PW <- diag(nrow=ntime, ncol=ntime)
  }

  # if(is.logical(sESS_correction)){
  #   if(!sESS_correction) do_sESS <- FALSE else stop('sESS_correction must be FALSE or a list of length 3.')
  # } else {
  #   do_sESS <- TRUE
  #   #in this case, sESS_correction must be a list of length 3
  #   if(!is.list(sESS_correction)){
  #     if(length(sESS_correction) != 3){
  #       stop('If sESS_correction is not FALSE, it must be a list of length 3.')
  #     }
  #   }
  # }
  #
  # if(do_sESS){
  #   P <- sESS_correction[[1]]
  #   FV <- sESS_correction[[2]]
  #   ind <- sESS_correction[[3]]
  #
  #   if (!requireNamespace("INLA", quietly = TRUE)) {
  #     stop(
  #       "Package \"INLA\" needed to for spatial ESS correction ",
  #       "Please install it at https://www.r-inla.org/download-install.",
  #       call. = FALSE
  #     )
  #   }
  #  mesh <- INLA::inla.mesh.create(loc = P, tv = FV)
  #
  #  rtime <- round(seq(1, ntime, length.out=12)); rtime <- rtime[2:11] #evenly spaced time intervals
  #  sESS <- rep(NA, 10)
  #   for(k in 1:10){
  #     time_k <- rtime[k]
  #     sESS[k] <- estimate.ESS(mesh, Y = E0[,time_k], ind = ind, trace=TRUE)
  #   }
  #   sESS <- mean(sESS)
  #   if(verbose) cat(paste0('\t Spatial ESS accounting for spatial correlation: ', round(sESS), '\n'))
  #   if(sESS >= nvox){
  #     sESS <- NULL
  #     sESS_correction <- FALSE
  #     if(verbose) cat(paste0('\t Skipping spatial ESS correction since estimated ESS >= nvox \n'))
  #   }
  # } else {
  #   sESS <- NULL
  # }

  #1. Compute initial estimates of posterior moments for A, S, V(S), tau^2

  mu_A <- PW %*% A0 #TxQ
  mu_S <- t(S0) #QxV
  cov_S <- array(0, dim = c(nICs, nICs, nvox)) #QxQxV
  for (v in 1:nvox) { cov_S[,,v] <- diag(S0_var[v,]) }
  E_SSt <- apply(cov_S, 1:2, sum) + mu_S %*% t(mu_S)
  E_SSt0 <- E_SSt #without ESS correction
  #if(do_sESS){ E_SSt <- E_SSt * sESS / nvox }

  #mu_tau2 <- apply(BOLD - t(mu_A %*% mu_S),1,var) #Vx1
  mu_tau2 <- mean((BOLD %*% PW - t(mu_A %*% mu_S))^2) #scalar
  BOLD <- PW %*% t(BOLD) #make the BOLD TxV to match the paper

  #pre-compute some stuff
  D_inv <- 1/template_var
  D_inv_S <- D_inv * template_mean
  BOLD2 <- sum(BOLD^2) #sum over all v,t

  #sample a bunch of Gamma variates for p(a|u)
  if(method_FC == 'VB1'){
    nu_a <- template_FC$nu + 1 - nICs
    u_samps <- stats::rgamma(nsamp_u, shape = nu_a/2, rate = nu_a/2)
  }

  # #save intermediate results
  # save_inter <- !is.null(eps_inter)
  # if(save_inter){
  #   results_inter <- vector('list', length(eps_inter))
  #   names(results_inter) <- paste0('epsilon_',eps_inter)
  #   next_eps <- eps_inter[1]
  # } else {
  #   results_inter <- NULL
  # }

  #2. Iteratively update approximate posteriors

  err <- 1000 #large initial value for difference between iterations
  #ELBO_vals <- rep(NA, maxiter) #keep track of ELBO at each iteration (convergence criterion)
  while (err > epsilon | iter <= miniter) {

    #if(verbose) cat(paste0('\t ~~~~~~~~~~~ VB ITERATION ', iter, ' ~~~~~~~~~~~ \n'))
    #t00 <- Sys.time()

    #a. UPDATE A

    #this function will return mu_A and E[A'A], which is needed for tau2 and for S

    if(method_FC == 'VB1'){

      #15 sec Q=25, V=90k, T=1200!
      A_new <- update_A(mu_tau2, mu_S, E_SSt,
                      BOLD, ntime, nICs, nvox,
                      u_samps, template_FC, final=FALSE,
                      usePar=FALSE) #for this function parallelization doesn't do much, and it runs very quickly anyway

    } else {

      #15 min for 50,000 samples from p(G) (Q=25, V=90k, T=1200)
      #13 min without parallelization
      A_new <- update_A_Chol(mu_tau2, mu_S, E_SSt,
                             BOLD, ntime, nICs, nvox,
                             template_FC, final=FALSE, exact=FALSE,
                             usePar=usePar,
                             verbose=verbose)
    }

    #1. Look at eigen and err, add to Notes, then comment out those lines
    #2. Run update_A() with final=TRUE and compare intervals

    mu_A_old <- mu_A
    mu_A <- A_new$mu_A
    change_A <- mean(abs(mu_A - mu_A_old)/(abs(mu_A_old)+0.1)) #add 0.1 to denominator to avoid dividing by zero
    E_AtA <- A_new$E_AtA
    E_AtA0 <- E_AtA #without ESS correction
    #if(tESS_correction){ E_AtA <- E_AtA * tESS / ntime }

    #b. UPDATE S
    S_new <- update_S(mu_tau2, mu_A, E_AtA, D_inv, D_inv_S, BOLD, nICs, nvox, ntime)
    change_S <- mean(abs(S_new$mu_S - mu_S)/(abs(mu_S)+0.1)) #add 0.1 to denominator to avoid dividing by zero
    mu_S <- S_new$mu_S
    E_SSt <- S_new$E_SSt
    E_SSt0 <- E_SSt #without ESS correction
    #if(do_sESS){ E_SSt <- E_SSt * sESS / nvox }

    #d. UPDATE tau^2

    tau2_new <- update_tau2(BOLD, BOLD2,
                            mu_A, E_AtA0, mu_S, E_SSt0,
                            prior_params, ntime, nvox)
    beta1_nvox <- tau2_new$beta1 #beta1 / nvox
    alpha1_nvox <- tau2_new$alpha1 #alpha1 / nvox
    mu_tau2_new <- beta1_nvox/(alpha1_nvox - 1/nvox) #numerator and denominator divided by nvox for computational stability
    change_tau2 <- abs(mu_tau2_new - mu_tau2)/(mu_tau2+0.1) #add 0.1 to denominator to avoid dividing by zero
    #change_tau2 <- sqrt(crossprod(c(mu_tau2_new - mu_tau2))) #same as in SQUAREM
    mu_tau2 <- mu_tau2_new

    #if (verbose) print(Sys.time() - t00)

    #change in estimates
    change <- c(change_A, change_S, change_tau2)
    err <- max(change)
    change <- format(change, digits=3, nsmall=3)
    if(verbose) cat(paste0('\t Iteration ',iter, ': Difference is ',change[1],' for A, ',change[2],' for S, ',change[3],' for tau2 \n'))

    # #ELBO
    # ELBO_vals[iter] <- compute_ELBO(mu_S, cov_S, cov_A, cov_alpha, template_mean, template_var, ntime)
    # if(iter == 1) err2 <- (ELBO_vals[iter] - ELBO_init)/ELBO_init
    # if(iter > 1) err2 <- (ELBO_vals[iter] - ELBO_vals[iter - 1])/ELBO_vals[iter - 1]
    # if(verbose) cat(paste0('Iteration ',iter, ': Change in ELBO = ',round(err2, 7),', Change in Params = ', round(err, 7),'\n'))

    # #Save intermediate result?
    # if(save_inter){
    #   ##only consider convergence for positive change (no longer relevant because err based on parameters, not ELBO)
    #   #if(err > 0){
    #   if(err < max(eps_inter)){ #if we have reached one of the intermediate convergence thresholds, save results
    #     which_eps <- max(which(err < eps_inter)) #most stringent convergence level met
    #     if(is.null(results_inter[[which_eps]])){ #save intermediate result at this convergence level if we haven't already
    #       results_inter[[which_eps]] <- list(S = t(mu_S), A = mu_A, tau2_mean = mu_tau2_new, error=err, numiter=iter)
    #     }
    #     #}
    #   }
    # }

    ### Move to next iteration
    iter <- iter + 1
    if (iter > maxiter) {
      success <- 0
      warning(paste0('Failed to converge within ', maxiter,' iterations'))
      break() #exit loop
    }
  } #end iterations

  ### Inference on cov_A

  if(method_FC == 'VB1'){
    A_new <- update_A(mu_tau2, mu_S, E_SSt,
                      BOLD, ntime, nICs, nvox,
                      u_samps, template_FC, final=TRUE, usePar=FALSE)

  } else {
    A_new <- update_A_Chol(mu_tau2, mu_S, E_SSt,
                           BOLD, ntime, nICs, nvox,
                           template_FC, final=TRUE, exact=TRUE,
                           usePar=usePar,
                           verbose=verbose)
  }

  mu_A <- A_new$mu_A
  #inference on FC(A) using samples from (a_t|Y,u) (VB1) or (a_t|Y,G) (VB2)
  FC_samp <- A_new$FC_samp
  FC_samp <- abind::abind(FC_samp, along=3)
  FC_mean <- apply(FC_samp, 1:2, mean, na.rm=TRUE)
  alpha <- 1 - CI_FC
  FC_LB <- apply(FC_samp, 1:2, quantile, alpha/2, na.rm=TRUE) #LB of 95% credible interval
  FC_UB <- apply(FC_samp, 1:2, quantile, 1-alpha/2, na.rm=TRUE) #UB of 95% credible interval

  #obtain SE(S)
  E_AtA <- A_new$E_AtA
  E_AtA0 <- E_AtA #without ESS correction
  #if(tESS_correction){ E_AtA <- E_AtA * tESS / ntime }
  S_new <- update_S(mu_tau2, mu_A, E_AtA, D_inv, D_inv_S, BOLD, nICs, nvox, ntime, final=TRUE)
  mu_S <- S_new$mu_S
  cov_S <- S_new$cov_S
  cov_S_list <- lapply(seq(dim(cov_S)[3]), function(x) cov_S[ , , x])
  subjICse <- sqrt(sapply(cov_S_list, diag))

  list(
    subjICmean = t(mu_S),
    subjICse = t(subjICse),
    S_cov = cov_S,
    A = mu_A,
    FC = list(mean = FC_mean, LB = FC_LB, UB = FC_UB, level = CI_FC),
    tau2_mean = mu_tau2,
    success_flag=success,
    error=err,
    numiter=iter-1,
    #ELBO=ELBO_vals[1:(iter-1)],
    #results_inter = results_inter,
    template_mean = template_mean,
    template_var = template_var,
    template_FC = template_FC,
    method_FC = method_FC,
    PW = PW
  )
}

#' Update A for VB FC Template ICA using IW prior on FC
#'
#' @param mu_tau2,mu_S,E_SSt Most recent estimates of posterior
#' moments for these variables.
#' @param BOLD (\eqn{T \times V} matrix) preprocessed fMRI data.
#' @param ntime,nICs,nvox Number of timepoints, number of ICs, number of locations
#' @param u_samps Samples from Gamma for multivariate t prior on A
#' @param template_FC IW parameters for FC template
# @param sESS Spatial effective sample size
#' @param final If TRUE, return cov_A. Default: \code{FALSE}
#' @param usePar Parallelize the computation? Default: \code{FALSE}. Can be the
#' number of cores to use or \code{TRUE}, which will use the number available minus two.
#'
#' @return List of length two: \code{mu_A} (TxQ) and \code{E_AtA} (QxQ).
#'
#' @keywords internal
update_A <- function(
  mu_tau2, mu_S, E_SSt,
  BOLD, ntime, nICs, nvox,
  u_samps, template_FC,
  #sESS=NULL,
  usePar=TRUE,
  final=FALSE){

  nu_a <- template_FC$nu + 1 - nICs
  psi0_inv <- solve(template_FC$psi)
  nu_a_psi0_inv <- nu_a * psi0_inv

  nU <- length(u_samps)
  E_Cov_A_u <- matrix(0, nICs, nICs) # estimate via MC using u_samps
  E_A_u <- array(0, dim=c(length(u_samps), nICs, ntime)) #for estimating Cov(E(a|u))
  tmp1 <- (1/mu_tau2) * E_SSt #QxQ
  tmp2 <- (1/mu_tau2) * (mu_S %*% t(BOLD)) #QxT (a sum over v)
  #if(!is.null(sESS)) tmp2 <- tmp2 * sESS/nvox #correct the sum over v for ESS
  if(final) { FC_samp <- array(0, dim=c(nICs, nICs, nU)) }

  #define function to compute E[a_t|u], Cov(a_t|u), and draw a sample if final=TRUE
  A_u <- function(u_val, final){

    tmp1i <- tmp1 + u_val * nu_a_psi0_inv
    # tmp1i_chol <- chol(tmp1i) #gives an upper triangular matrix
    # tmp1i_chol_inv <- backsolve(tmp1i_chol, x = diag(nICs))
    # Cov_A_ui2 <- (tmp1i_chol_inv %*% t(tmp1i_chol_inv))
    Cov_A_ui <- solve(tmp1i)
    #all.equal(Cov_A_ui, Cov_A_ui2) #TRUE
    #collect terms to compute Cov(E[A|u])
    mu_A_ui <- Cov_A_ui %*% tmp2 #QxT

    if(final==TRUE) {
      Cov_A_ui_chol <- chol(Cov_A_ui) #UT cholesky factor of covariance matrix
      Z_samp <- matrix(rnorm(nICs*ntime, mean = 0, sd = 1), nrow=nICs, ncol=ntime)
      A_u_samp <- (t(Cov_A_ui_chol) %*% Z_samp) + mu_A_ui #more efficient than mvrnorm because cov(a_t|G_k)=V_k is same for all a_t, t=1,...,T
      FC_samp <- cor(t(A_u_samp))
    } else {
      FC_samp <- NULL
    }

    return(list(mu = mu_A_ui, cov = Cov_A_ui, FC_samp = FC_samp))
  }

  if(!usePar){

    for(ii in 1:nU){
      ui <- u_samps[ii]
      #iteratively compute E_u[Cov(A|u)]
      A_ui <- A_u(ui, final) #list(mu = mu_A_ui, cov = Cov_A_ui, FC_samp = FC_samp)

      E_Cov_A_u <- E_Cov_A_u + A_ui$cov
      E_A_u[ii,,] <- A_ui$mu
      if(final==TRUE) FC_samp[,,ii] <- A_ui$FC_samp
    }

  } else {

    #parallel version with foreach
    `%dopar%` <- foreach::`%dopar%`
    A_ui_list <- foreach::foreach(ii = 1:nU) %dopar% {
      A_u(u_samps[ii], final)
    }
    E_Cov_A_u <- matrix(rowSums(sapply(A_ui_list, function(x) x$cov)), nrow = nICs, ncol = nICs)
    E_A_u <- abind::abind(lapply(A_ui_list, function(x) x$mu), along=0)
    if(final==TRUE) FC_samp <- abind::abind(lapply(A_ui_list, function(x) x$FC_samp), along=3)
  }

  E_Cov_A_u <- E_Cov_A_u/nU
  mu_A <- t(E_Cov_A_u %*% tmp2) #TxQ
  Cov_A_part1 <- E_Cov_A_u #common over all t
  Cov_A_part2 <- apply(E_A_u, 3, cov, simplify=TRUE) #QxQ cov matrices get vectorized
  Cov_A_part2_sum <- matrix(apply(Cov_A_part2, 1, sum), nICs, nICs) #sum over t for E[A'A], reform into matrix
  Cov_A_sum <- Cov_A_part1 * ntime + Cov_A_part2_sum

  ### Constrain each column of A to have var=1 and mean=0
  sd_A <- apply(mu_A, 2, sd)
  D_A <- diag(1/sd_A) #use this below to correct A and cov(A) for unit-var A
  mu_A <- scale(mu_A)
  Cov_A_sum <- D_A %*% Cov_A_sum %*% D_A
  E_AtA <- Cov_A_sum + t(mu_A) %*% mu_A

  result <- list(mu_A = mu_A, E_AtA = E_AtA)
  if(final) { result$FC_samp <- FC_samp } #this contains Cor(A), where a_t is a sample from a_t|u for each u~Gamma

  return(result)
}

#' Update A for VB FC Template ICA using Cholesky prior for FC
#'
#' @param mu_tau2,mu_S,E_SSt Most recent estimates of posterior
#' moments for these variables.
#' @param BOLD (\eqn{T \times V} matrix) preprocessed fMRI data.
#' @param ntime,nICs,nvox Number of timepoints, number of ICs, number of locations
#' @param template_FC IW parameters for FC template
# @param sESS Spatial effective sample size
#' @param final If TRUE, return cov_A. Default: \code{FALSE}
#' @param exact If FALSE, at intermediate steps (final = \code{FALSE}) will attempt to
#' approximate V_k as described in the appendix of Mejia et al. 2024. Setting
#' \code{exact = TRUE} may be helpful when using spatial ESS correction,
#' which can lead to a lower rate of usable prior samples for the approximation.
#' @param usePar Parallelize the computation? Default: \code{FALSE}. Can be the
#' number of cores to use or \code{TRUE}, which will use the number available minus two.
#' @param verbose If \code{TRUE}, display progress of algorithm. Default: \code{FALSE}.
#'
#' @return List of length two: \code{mu_A} (TxQ) and \code{E_AtA} (QxQ).
#' @importFrom stats sd cor
#'
#' @keywords internal
update_A_Chol <- function(mu_tau2, mu_S, E_SSt,
                          BOLD, ntime, nICs, nvox,
                          template_FC, #sESS=NULL,
                          final=FALSE, exact=FALSE,
                          usePar=TRUE, verbose=FALSE){

  #step 1: approximate V_k for each FC sample G_k
  #step 2: Compute E[a_t|G] for each sample G (this requires V_k)
  #step 3: Estimate mean and cov of E[a_t|G] over G to get E[a_t] and V(a_t)
  #step 4: Compute E[A'A]
  #step 5: If final=TRUE, generate samples of a_t -> A -> Cov(A) for inference

  #preliminaries
  Emat <- 1/mu_tau2 * E_SSt
  E_chol <- chol(Emat) #UT Cholesky factor R: E = R'R
  E_chol_inv <- backsolve(E_chol, x = diag(nICs)) #R^(-1)
  E_inv <- tcrossprod(E_chol_inv) # E^(-1) = R^(-1) R^(-T)
  maxeig_E_inv <- eigen(E_inv, only.values = TRUE)$values[1]

  #more preliminaries
  tmp <- 1/mu_tau2 * (mu_S %*% t(BOLD)) #QxT -- ingredient for E[a_t|G] (a sum over v)
  #if(!is.null(sESS)) tmp <- tmp * sESS / nvox

  #organize max eigenvalues by pivot
  nP <- length(template_FC$pivots)
  maxeig_mat <- matrix(template_FC$FC_samp_maxeig, ncol=nP)

  # max_eigen <- err <- c()
  #collect samples of FC matrices
  # if(final==TRUE) { FC_samp <- array(NA, dim = c(nICs, nICs, nG))}

  fun_pivot <- function(pp, final, exact=FALSE){

    #setup
    max_eigen <- c()
    o_p <- order(template_FC$pivots[[pp]])
    nK <- nrow(template_FC$FC_samp_cholinv[[pp]])
    if(final) { FC_samp <- array(NA, dim = c(nICs, nICs, nK)) } else {FC_samp <- NULL}

    #loop over samples
    mu_A_G_pp <- array(NA, dim = c(nICs, ntime, nK)) # collect E[a_t|G] for samples G
    V_p_mean <- array(0, dim = c(nICs, nICs)) #to sequentially compute the mean of V_k

    #first, check how many samples have max eigenvalue < 1 (required for approximation of inverse)?
    which_valid <- rep(TRUE, nK) #keep track of how many valid samples
    if(!final | !exact){
      for(kk in 1:nK){
        eig_pk <- maxeig_E_inv * maxeig_mat[kk,pp] #avoid the need to run eigen() during model estimation
        max_eigen <- c(max_eigen, eig_pk)
        if(eig_pk >= 1) {
          which_valid[kk] <- FALSE
          #next() #will not use this sample of G for mu_A and cov_A
        }
      } #end loop over k
      if(mean(which_valid) < 0.9) {
        warning(paste0('< 90% of Cholesky samples from pivot ', pp, ' cannot be used for approximation of V_k. Setting exact = TRUE.'))
        exact <- TRUE
        which_valid <- rep(TRUE, nK)
      }
    }
    nK_valid <- sum(which_valid) #this will be nK when exact = TRUE

    for(kk in 1:nK){

      #a) obtain R_k^(-1) where G_k^(-1) = R_k^(-1)R_k^(-T). Note that R_k^(-1) is not actually UT because it has been un-pivoted.

      R_pk_inv_UT <- template_FC$FC_samp_cholinv[[pp]][kk,] #vectorized inverse of pivoted Cholesky UT factor
      R_pk_inv <- (UT2mat(R_pk_inv_UT))[o_p,] #un-pivot by permuting rows (not columns because inverse)
      G_pk_inv <- tcrossprod(R_pk_inv) #this is G_k^(-1)
      #R_k <- UT2mat(template_FC$Chol_samp[[pp]][kk,])
      #G_k <- crossprod(R_k[,o_p]) #permute columns of UT Cholesky factor

      #c) approximate or calculate V_k

      if(final | exact){

        ### For final estimation or when exact=TRUE, compute V_k exactly
        V_pk <- solve(Emat + G_pk_inv)

      } else {

        ### For intermediate estimation, approximate V_k unless exact = TRUE
        #b) compute the max eigenvalue of E^(-1) %*% R_k^(-1) %*% R_k^(-T) and save
        #B_pk <- t(E_chol_inv) %*% G_pk_inv %*% E_chol_inv
        #eig_pk <- eigen(B_pk, only.values = TRUE)$values[1]

        if(!which_valid[kk]) next()
        B_pk <- t(E_chol_inv) %*% G_pk_inv %*% E_chol_inv
        V_pk <- E_chol_inv %*% ( diag(nICs) - B_pk + (B_pk %*% B_pk) ) %*% t(E_chol_inv) #use approximation (I+B)^(-1) = I - A + A^2
      }

      V_p_mean <- V_p_mean + V_pk #sum over samples k (later will divide by nK_valid)

      #d) save E[a_t|G] to estimate its mean and covariance over G (do this within pivots for computational efficiency)

      mu_pk <- V_pk %*% tmp # QxT
      mu_A_G_pp[,,kk] <- mu_pk

      #x) if final=TRUE, take a sample of (a_t|G) ~ N(mu_pk, V_pk) for each G_k --> sample of A|G --> sample of cov(A)|G
      if(final==TRUE) {
        V_pk_chol <- chol(V_pk)
        Z_samp <- matrix(rnorm(nICs*ntime, mean = 0, sd = 1), nrow=nICs, ncol=ntime)
        A_G_samp <- (t(V_pk_chol) %*% Z_samp) + mu_pk #more efficient than mvrnorm because cov(a_t|G_k)=V_k is same for all a_t, t=1,...,T
        FC_samp[,,kk] <- cor(t(A_G_samp))
      }

    } #end loop over samples for current pivot

    #if(nK_valid <= 1) stop('Only one or zero Cholesky samples from current pivot can be used, contact developer')
    #if(nK_valid < nK*0.5) { warning(paste0('Over half of Cholesky samples from pivot ', pp, ' cannot be used for approximation of V_k. Setting exact = TRUE.'))

    V_p_mean <- V_p_mean/nK_valid

    # e) estimate E_G[E[a_t|G]] and Cov_G(E[a_t|G]) for the current pivot

    #E_G[E[a_t|G]]
    mu_A_pp <- apply(mu_A_G_pp, 1:2, mean, na.rm=TRUE) #set na.rm=TRUE to ignore the NAs due to samples with max eig > 1
    #mu_A_bypivot[[pp]] <- mu_A_pp

    #Cov_G(E[a_t|G])
    mu_A_G_pp_ctr <- mu_A_G_pp - array(rep(mu_A_pp, nK), dim=dim(mu_A_G_pp)) #center across samples
    tmp1 <- apply(mu_A_G_pp_ctr, 2:3, tcrossprod, simplify=TRUE) #compute xx' for each time point and sample, where x is Qx1 (each QxQ matrix will be vectorized)
    tmp2 <- apply(tmp1, 1:2, sum, na.rm=TRUE)/(nK_valid - 1) #sum over samples, divide by K-1 --> Q*Q x T
    cov_A_G_pp <- array(tmp2, dim=c(nICs, nICs, ntime)) #reshape to unvectorize cov for each time point --> Q x Q x T
    cov_A_sum_pp <- V_p_mean * ntime + apply(cov_A_G_pp, 1:2, sum) #sum over t

    return(list(mu = mu_A_pp,
                cov = cov_A_sum_pp,
                FC_samp = FC_samp, #NULL if !final
                nK_valid = nK_valid,
                max_eigen = max_eigen,
                exact = exact))
  }

  #loop over pivots
  #if(verbose) cat(paste0('\t Looping over ', nP, ' Cholesky pivots: '))

  if(!usePar){

    mu_A_bypivot <- cov_A_sum_bypivot <- vector('list', length = nP)
    if(final) FC_samp <- vector('list', length = nP)
    nK_valid_all <- 0
    max_eigen <- c()
    for(pp in 1:nP){

      if(verbose & (pp/10 == round(pp/10))) cat(paste0(pp,' '))

      #do function for current pivot
      result_pp <- fun_pivot(pp, final = final, exact = exact)

      #aggregate results
      nK_valid_all <- nK_valid_all + result_pp$nK_valid
      mu_A_bypivot[[pp]] <- result_pp$mu
      cov_A_sum_bypivot[[pp]] <- result_pp$cov
      if(final) FC_samp[[pp]] <- result_pp$FC_samp
      max_eigen <- c(max_eigen, result_pp$max_eigen)

    } #end loop over pivots

    if(final) FC_samp <- abind::abind(FC_samp, 3)

  } else {

    #parallel version with foreach
    `%dopar%` <- foreach::`%dopar%`
    result_list <- foreach::foreach(pp = 1:nP) %dopar% {
      fun_pivot(pp, final)
    }
    nK_valid_all <- sum(sapply(result_list, function(x) x$nK_valid))
    mu_A_bypivot <- lapply(result_list, function(x) x$mu)
    cov_A_sum_bypivot <- lapply(result_list, function(x) x$cov)
    if(final) FC_samp <- lapply(result_list, function(x) x$FC_samp)
    max_eigen <- unlist(lapply(result_list, function(x) x$max_eigen))

  }

  nG <- sum(sapply(template_FC$FC_samp_cholinv, nrow)) #total number of samples
  if(verbose & (nK_valid_all < nG)) cat(paste0('(',nG - nK_valid_all, ' samples dropped among ', nG, ') \n'))

  # f) average over pivots to get mu_A, cov_A
  mu_A <- t(Reduce('+', mu_A_bypivot)/nP)
  Cov_A_sum <- Reduce('+', cov_A_sum_bypivot)/nP

  # g) rescale so mu_A has mean 0 and var 1 before computing E[A'A]
  sd_A <- apply(mu_A, 2, sd)
  D_A <- diag(1/sd_A) #use this to correct A and cov(A) for unit-var A
  mu_A <- scale(mu_A)
  Cov_A_sum <- D_A %*% Cov_A_sum %*% D_A
  E_AtA <- Cov_A_sum + t(mu_A) %*% mu_A

  result <- list(mu_A = mu_A,
                 E_AtA = E_AtA,
                 max_eigen = max_eigen,
                 #err = err,
                 cov_A_sum_bypivot = cov_A_sum_bypivot,
                 nK_valid_all = nK_valid_all)
  if(final) { result$FC_samp <- FC_samp } #this contains Cor(A), where a_t is a sample from a_t|G for each G
  return(result)

}



#' Update S for VB FC Template ICA
#'
#' @param mu_tau2,mu_A,E_AtA, Most recent posterior estimates
#' @param D_inv,D_inv_S Some pre-computed quantities.
#' @param BOLD (\eqn{T \times V} matrix) preprocessed fMRI data.
#' @param nICs,nvox,ntime Number of ICs, number of data locations, and time series length
# @param tESS Effective sample size
#' @param final If TRUE, return cov_S. Default: \code{FALSE}
#'
#' @return List of length three: \code{mu_S} (QxV), \code{E_SSt} (QxQ), and \code{cov_S} (QxQxV).
#'
#' @keywords internal
update_S <- function(
  mu_tau2, mu_A, E_AtA,
  D_inv, D_inv_S,
  BOLD,
  nICs, nvox, ntime,
  #tESS = NULL,
  final = FALSE){

  tmp1 <- (1/mu_tau2) * E_AtA

  cov_S <- array(0, dim = c(nICs, nICs, nvox))
  for(v in 1:nvox){ cov_S[,,v] <- solve(tmp1 + diag(D_inv[v,])) }

  mu_S <- array(0, dim = c(nICs, nvox)) #QxV
  tmp2 <- (1/mu_tau2) * (t(mu_A) %*% BOLD) #this is a sum over t=1,...,T
  #if(!is.null(tESS)) tmp2 <- tmp2 * tESS / ntime
  for(v in 1:nvox){
    #sum over t=1...T part
    #D_v_inv <- diag(1/template_var[v,])
    mu_S[,v] <- cov_S[,,v] %*% (tmp2[,v] + D_inv_S[v,] )
  }

  E_SSt <- apply(cov_S, 1:2, sum) + mu_S %*% t(mu_S)

  result <- list(mu_S=mu_S, E_SSt=E_SSt)
  if(final) result$cov_S <- cov_S
  return(result)
}


#' Update tau for VB FC Template ICA
#'
#' @param BOLD (\eqn{T \times V} matrix) preprocessed fMRI data.
#' @param BOLD2 A precomputed quantity, \code{sum(BOLD^2)}
#' @param mu_A,E_AtA,mu_S,E_SSt Current posterior estimates
#' @param prior_params Alpha and beta parameters of IG prior on \eqn{\tau^2}
#'  (error variance).
#' @param ntime,nvox Number of timepoints in data and the number of data
#'  locations.
# @param tESS Effective sample size
#'
#' @return List of length two: \code{alpha1} and \code{beta1}.
#'
#' @keywords internal
update_tau2 <- function(
  BOLD, BOLD2, mu_A, E_AtA, mu_S, E_SSt, prior_params, ntime, nvox){ #}, tESS = NULL){

  alpha0 <- prior_params[1]
  beta0 <- prior_params[2]

  #these are alpha and beta divided by nvox for computational stability
  # if(!is.null(tESS)){
  #   alpha1_nvox <- alpha0/nvox + tESS/2
  #   beta1_nvox <- beta0/nvox +
  #     (1/2)*BOLD2/nvox * (tESS / ntime) -
  #     sum(BOLD * mu_A %*% mu_S/nvox) * (tESS / ntime) +
  #     (1/2)*sum(diag(E_AtA %*% E_SSt)/nvox)
  # } else {
    alpha1_nvox <- alpha0/nvox + ntime/2
    beta1_nvox <- beta0/nvox +
      (1/2)*BOLD2/nvox -
      sum(BOLD * mu_A %*% mu_S/nvox) +
      (1/2)*sum(diag(E_AtA %*% E_SSt)/nvox)
  #}

  # #pre-compute sum over t inside of trace (it doesn't involve v)
  # tmp2a <- t(mu_A) %*% mu_A + ntime*cov_A
  #
  # tmp1 <- t(BOLD) %*% mu_A
  # for(v in 1:nvox){
  #   #sums over t=1,...,T in beta_v
  #   tmp1v <- crossprod(tmp1[v,], mu_S[,v])
  #   tmp2b <- tcrossprod(mu_S[,v]) + cov_S[,,v]
  #   Tr_v <- sum(diag(tmp2a %*% tmp2b)) #TO DO: make trace more efficient with element-wise multiplication and colSums ?
  #   beta1[v] <- beta0 + (1/2)*BOLD2[v] - tmp1v + (1/2)*Tr_v
  # }

  list(alpha1=alpha1_nvox, beta1=beta1_nvox)
}

#' Bdiag m2
#'
#' @param mat a k x k 'matrix'
#' @param N how many times to repeat \code{mat}
#' @return a sparse (N*k x N*k) matrix of class \code{"dgCMatrix"}.
bdiag_m2 <- function(mat, N) {
  # Fast version of Matrix :: .bdiag() -- for the case of *many identical*  (k x k) matrices:
  k <- nrow(mat)
  if(N * k > .Machine$integer.max)
    stop("resulting matrix too large; would be  M x M, with M=", N*k)
  M <- as.integer(N * k)
  x <- rep(c(mat), N)
  ## result: an   M x M  matrix
  new("dgCMatrix", Dim = c(M,M),
      ## 'i :' maybe there's a faster way (w/o matrix indexing), but elegant?
      i = as.vector(matrix(0L:(M-1L), nrow=k)[, rep(seq_len(N), each=k)]),
      p = k * 0L:M,
      x = as.double(x))
}

#' Estimation of effective sample size
#'
#' @param mesh INLA mesh
#' @param Y data
#' @param ind index of the data locations in the mesh
# @param use_inla should INLA be used to compute standard deviations? If FALSE, the
# standard deviations is approximated as being constant over the locations
#' @param trace If \code{TRUE}, a formula based on the trace is used, otherwise the inverse correlation is used
#' @details
#' The functions computes the effective sample size as \eqn{trace(Q^{-1})/trace(Q*Q)}
#' as in Bretherton et al. (1999), Journal of Climate.
#'
#' @return Estimate of the effective sample size
estimate.ESS <- function(mesh, Y, ind = NULL, trace=FALSE) {

  fem <- INLA::inla.mesh.fem(mesh)
  #fem <- fmesher:::fm_fem(mesh)

  res <- optim(c(log(1),log(1)), lik, Y = Y, C = fem$c0, G = fem$g1, ind = ind)
  Q <- exp(res$par[2])^2*(exp(res$par[1])^2*fem$c0 + fem$g1)
  if(!is.null(ind)) {
    ind2 <- setdiff(1:dim(Q)[1], ind)
    Q <- Q[ind,ind] - Q[ind,ind2]%*%solve(Q[ind2,ind2],Q[ind2,ind])
  }

  Sigma <- solve(Q)
  if(trace){
    return(sum(diag(Sigma))^2/sum(Sigma^2))
  } else {
    Cor <- cov2cor(Sigma)
    return(sum(Cor))
  }
}

#' Compute likelihood in SPDE model for ESS estimation
#'
#' @param theta Value of hyperparameters
#' @param Y Data vector
#' @param G SPDE G matrix
#' @param C SPDE C matrix
#' @param ind Indices of data locations in the mesh
#'
#' @return Log likelihood value
#'
lik <- function(theta, Y, G,C,ind = NULL) {
  kappa <- exp(theta[1])
  tau <- exp(theta[2])
  Q <- tau^2*(kappa^2*C + G)

  if(!is.null(ind)) {
    ind2 <- setdiff(1:dim(Q)[1], ind)
    Q <- Q[ind,ind] - Q[ind,ind2]%*%solve(Q[ind2,ind2],Q[ind2,ind])
  }
  R <- chol(Q)
  ld <- c(determinant(R, logarithm = TRUE, sqrt = TRUE)$modulus)
  lik <- as.double(ld - 0.5*t(Y)%*%Q%*%Y)

  return(-lik)
}


