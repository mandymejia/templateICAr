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
#' @param eps_inter Intermediate values of epsilon at which to save results (used
#'  to assess benefit of more stringent convergence rules). Default:
#'  \code{NULL} (do not save). These values should be in decreasing order
#'  (larger to smaller error) and all values should be between zero and
#'  \code{epsilon}.
#' @param verbose If \code{TRUE}, display progress of algorithm.
#' Default: \code{FALSE}.
#'
#' @return A list of computed values, including the final parameter estimates.
#'
#' @importFrom fMRItools is_posNum
#' @keywords internal
#'
VB_FCtemplateICA <- function(
  template_mean, #VxQ
  template_var, #VxQ
  template_FC,
  method_FC = c('VB1','VB2'),
  prior_params=c(0.001, 0.001),
  BOLD, #VxT
  A0, S0, S0_var,
  maxiter=100,
  miniter=3,
  epsilon=0.001,
  eps_inter=NULL,
  verbose=FALSE){

  stopifnot(length(prior_params)==2)
  stopifnot(is.numeric(prior_params))
  stopifnot(is_posNum(maxiter, "numeric") && maxiter==round(maxiter))
  stopifnot(is_posNum(miniter, "numeric") && miniter==round(miniter))
  stopifnot(miniter <= maxiter)
  stopifnot(is_posNum(epsilon))
  if (!is.null(eps_inter)) {
    stopifnot(is.numeric(eps_inter) && all(diff(eps_inter) < 0))
    stopifnot(eps_inter[length(eps_inter)]>0 && eps_inter[1]>epsilon)
  }
  method_FC <- match.arg(method_FC, c("VB1", "VB2"))

  if (!all.equal(dim(template_var), dim(template_mean))) stop('The dimensions of template_mean and template_var must match.')
  ntime <- ncol(BOLD) #length of timeseries
  nvox <- nrow(BOLD) #number of brain locations
  if (ntime > nvox) warning('More time points than brain locations. Are you sure?')
  if (nrow(template_mean) != nvox) stop('Templates and BOLD must have the same number of brain locations (columns).')

  nICs <- nrow(template_FC$psi)   #number of ICs
  if (ncol(template_mean) != nICs) stop('template_FC is incompatible with template_mean & template_var. Check number of ICs in each.')
  if (nICs > nvox) stop('Cannot estimate more ICs than brain locations.')
  if (nICs > ntime) stop('Cannot estimate more ICs than time points.')

  iter <- 1
  success <- 1
  template_var[template_var < 1e-6] <- 1e-6 #to prevent problems when inverting covariance (this should rarely/never happen with NN variance)

  #1. Compute initial estimates of posterior moments for G, A, S, V(S), tau^2

  mu_A <- A0 #TxQ
  mu_S <- t(S0) #QxV
  cov_S <- array(0, dim = c(nICs, nICs, nvox)) #QxQxV
  for (v in 1:nvox) { cov_S[,,v] <- diag(S0_var[v,]) }
  E_SSt <- apply(cov_S, 1:2, sum) + mu_S %*% t(mu_S)

  #mu_tau2 <- apply(BOLD - t(mu_A %*% mu_S),1,var) #Vx1
  mu_tau2 <- mean((BOLD - t(mu_A %*% mu_S))^2) #scalar
  BOLD <- t(BOLD) #make the BOLD TxV to match the paper

  #pre-compute some stuff
  D_inv <- 1/template_var
  D_inv_S <- D_inv * template_mean
  BOLD2 <- sum(BOLD^2) #sum over all v,t

  #sample a bunch of Gamma variates for p(a|u)
  nu_a <- template_FC$nu + 1 - nICs
  u_samps <- stats::rgamma(10000, shape = nu_a/2, rate = nu_a/2)

  #save intermediate results
  save_inter <- !is.null(eps_inter)
  if(save_inter){
    results_inter <- vector('list', length(eps_inter))
    names(results_inter) <- paste0('epsilon_',eps_inter)
    next_eps <- eps_inter[1]
  } else {
    results_inter <- NULL
  }

  #2. Iteratively update approximate posteriors

  err <- 1000 #large initial value for difference between iterations
  #ELBO_vals <- rep(NA, maxiter) #keep track of ELBO at each iteration (convergence criterion)
  while (err > epsilon & iter <= miniter) {

    if(verbose) cat(paste0(' ~~~~~~~~~~~~~~~~~~~~~ VB ITERATION ', iter, ' ~~~~~~~~~~~~~~~~~~~~~ \n'))
    t00 <- Sys.time()

    #a. UPDATE A

    #this function will return mu_A and E[A'A], which is needed for tau2 and for S

    #only ~5s for Q=25, V=90k, T=1200!
    if(method_FC == 'VB1'){
      A_new <- update_A(mu_tau2, mu_S, E_SSt,
                      BOLD, ntime, nICs,
                      u_samps, template_FC, final=FALSE)

      mu_A_old <- mu_A
      mu_A <- A_new$mu_A
      change_A <- mean(abs(mu_A - mu_A_old)/(abs(mu_A_old)+0.1)) #add 0.1 to denominator to avoid dividing by zero
      E_AtA <- A_new$E_AtA
    } else {

      # [TO DO] implement Cholesky prior-based estimation of G

    }

    #b. UPDATE S
    S_new <- update_S(mu_tau2, mu_A, E_AtA, D_inv, D_inv_S, BOLD, nICs, nvox)
    change_S <- mean(abs(S_new$mu_S - mu_S)/(abs(mu_S)+0.1)) #add 0.1 to denominator to avoid dividing by zero
    mu_S <- S_new$mu_S
    E_SSt <- S_new$E_SSt

    #d. UPDATE tau^2

    tau2_new <- update_tau2(BOLD, BOLD2,
                            mu_A, E_AtA, mu_S, E_SSt,
                            prior_params, ntime, nvox)
    beta1_nvox <- tau2_new$beta1 #beta1 / nvox
    alpha1_nvox <- tau2_new$alpha1 #alpha1 / nvox
    mu_tau2_new <- beta1_nvox/(alpha1_nvox - 1/nvox)
    change_tau2 <- abs(mu_tau2_new - mu_tau2)/(mu_tau2+0.1) #add 0.1 to denominator to avoid dividing by zero
    #change_tau2 <- sqrt(crossprod(c(mu_tau2_new - mu_tau2))) #same as in SQUAREM
    mu_tau2 <- mu_tau2_new

    if (verbose) print(Sys.time() - t00)

    #change in estimates
    change <- c(change_A, change_S, change_tau2)
    err <- max(change)
    change <- format(change, digits=3, nsmall=3)
    if(verbose) cat(paste0('Iteration ',iter, ':l0 Difference is ',change[1],' for A, ',change[2],' for S, ',change[3],' for tau2 \n'))

    # #ELBO
    # ELBO_vals[iter] <- compute_ELBO(mu_S, cov_S, cov_A, cov_alpha, template_mean, template_var, ntime)
    # if(iter == 1) err2 <- (ELBO_vals[iter] - ELBO_init)/ELBO_init
    # if(iter > 1) err2 <- (ELBO_vals[iter] - ELBO_vals[iter - 1])/ELBO_vals[iter - 1]
    # if(verbose) cat(paste0('Iteration ',iter, ': Change in ELBO = ',round(err2, 7),', Change in Params = ', round(err, 7),'\n'))

    #Save intermediate result?
    if(save_inter){
      ##only consider convergence for positive change (no longer relevant because err based on parameters, not ELBO)
      #if(err > 0){
      if(err < max(eps_inter)){ #if we have reached one of the intermediate convergence thresholds, save results
        which_eps <- max(which(err < eps_inter)) #most stringent convergence level met
        if(is.null(results_inter[[which_eps]])){ #save intermediate result at this convergence level if we haven't already
          results_inter[[which_eps]] <- list(S = t(mu_S), A = mu_A, tau2_mean = mu_tau2_new, error=err, numiter=iter)
        }
        #}
      }
    }

    ### Move to next iteration
    iter <- iter + 1
    if (iter > maxiter) {
      success <- 0
      warning(paste0('Failed to converge within ', maxiter,' iterations'))
      break() #exit loop
    }
  } #end iterations

  ### Inference on cov_A

  A_new <- update_A(mu_tau2, mu_S, E_SSt,
                    BOLD, ntime, nICs,
                    u_samps, template_FC, final=TRUE)
  mu_A <- A_new$mu_A
  #inference on Cov(A) using samples Cov(A|Y,u), u ~ Gamma
  Cov_A_samp <- A_new$Cov_A_samp
  Cov_A_mean <- apply(Cov_A_samp, 1:2, mean)
  Cov_A_LB <- apply(Cov_A_samp, 1:2, quantile, 0.025) #LB of 95% credible interval
  Cov_A_UB <- apply(Cov_A_samp, 1:2, quantile, 0.975) #UB of 95% credible interval

  #obtain SE(S)
  S_new <- update_S(mu_tau2, mu_A, A_new$E_AtA, D_inv, D_inv_S, BOLD, nICs, nvox, final=TRUE)
  mu_S <- S_new$mu_S
  cov_S <- S_new$cov_S
  cov_S_list <- lapply(seq(dim(cov_S)[3]), function(x) cov_S[ , , x])
  subjICse <- sqrt(sapply(cov_S_list, diag))

  list(
    subjICmean = t(mu_S),
    subjICse = t(subjICse),
    S_cov = cov_S,
    A = mu_A,
    A_cov = list(mean = Cov_A_mean, LB95 = Cov_A_LB, UB95 = Cov_A_UB),
    tau2_mean = mu_tau2,
    success_flag=success,
    error=err,
    numiter=iter-1,
    #ELBO=ELBO_vals[1:(iter-1)],
    results_inter = results_inter,
    template_mean = template_mean,
    template_var = template_var,
    template_FC = template_FC,
    method_FC = method_FC
  )
}

#' Update A for VB FC Template ICA
#'
#' @param mu_tau2,mu_S,E_SSt Most recent estimates of posterior
#' moments for these variables.
#' @param BOLD (\eqn{T \times V} matrix) preprocessed fMRI data.
#' @param ntime,nICs Number of timepoints, number of ICs
#' @param u_samps Samples from Gamma for multivariate t prior on A
#' @param template_FC IW parameters for FC template
#' @param final If TRUE, return cov_A. Default: \code{FALSE}
#'
#' @return List of length two: \code{mu_A} (TxQ) and \code{E_AtA} (QxQ).
#'
#' @keywords internal
update_A <- function(
  mu_tau2, mu_S, E_SSt,
  BOLD, ntime, nICs,
  u_samps, template_FC,
  final=FALSE){

  nu_a <- template_FC$nu + 1 - nICs
  psi0_inv <- solve(template_FC$psi)
  nu_a_psi0_inv <- nu_a * psi0_inv

  nU <- length(u_samps)
  E_Cov_A_u <- matrix(0, nICs, nICs) # estimate via MC using u_samps
  E_A_u <- array(0, dim=c(length(u_samps), nICs, ntime)) #for estimating Cov(E(a|u))
  tmp1 <- (1/mu_tau2) * E_SSt #QxQ
  tmp2 <- (1/mu_tau2) * mu_S %*% t(BOLD) #QxT
  if(final) {
    Cov_A_avg <- cov(t(tmp2))
    Cov_A_samp <- array(0, dim=c(nICs, nICs, nU))
  }

  for(ii in 1:nU){
    ui <- u_samps[ii]
    #iteratively compute E_u[Cov(A|u)]
    tmp1i <- tmp1 + ui * nu_a_psi0_inv
    # tmp1i_chol <- chol(tmp1i) #gives an upper triangular matrix
    # tmp1i_chol_inv <- backsolve(tmp1i_chol, x = diag(nICs))
    # Cov_A_ui2 <- (tmp1i_chol_inv %*% t(tmp1i_chol_inv))
    Cov_A_ui <- solve(tmp1i)
    #all.equal(Cov_A_ui, Cov_A_ui2) #TRUE
    E_Cov_A_u <- E_Cov_A_u + Cov_A_ui
    #collect terms to compute Cov(E[A|u])
    E_A_u[ii,,] <- Cov_A_ui %*% tmp2
    if(final) Cov_A_samp[,,ii] <- Cov_A_ui %*% Cov_A_avg %*% Cov_A_ui

    # #draw a sample from A --> Cov(A) (QxQ)
    # sampAi <- matrix(rnorm(ntime*nICs), nICs, ntime) #QxT matrix
    # sampAi_covA <- tmp1i_chol_inv %*% sampAi #multiply each Qx1 N(0,1) vec by A, where AAt = Cov(A|u)
    # #Cov(sampAi_covA[,1]) = A * I * t(A)
  }
  E_Cov_A_u <- E_Cov_A_u/length(u_samps)
  mu_A <- t(E_Cov_A_u %*% tmp2) #TxQ
  Cov_A_part1 <- E_Cov_A_u #common over all t
  Cov_A_part2 <- apply(E_A_u, 3, cov, simplify=TRUE)
  Cov_A_part2_sum <- matrix(apply(Cov_A_part2, 1, sum), nICs, nICs) #sum over t for E[A'A]
  Cov_A_sum <- Cov_A_part1 * ntime + Cov_A_part2_sum

  ### Constrain each column of A to have var=1 and mean=0
  sd_A <- apply(mu_A, 2, sd)
  D_A <- diag(1/sd_A) #use this below to correct A and cov(A) for unit-var A
  mu_A <- scale(mu_A)
  Cov_A_sum <- D_A %*% Cov_A_sum %*% D_A
  E_AtA <- Cov_A_sum + t(mu_A) %*% mu_A

  result <- list(mu_A = mu_A,
                 E_AtA = E_AtA)
  if(final) {
    #result$cov_A <- Cov_A #this contains Cov(a_t) for t=1,...,T
    Cov_A_samp <- apply(Cov_A_samp, 3, function(x) D_A %*% x %*% D_A)
    Cov_A_samp <- array(Cov_A_samp, dim=c(nICs, nICs, nU))
    result$Cov_A_samp <- Cov_A_samp #this contains Cov(a|u) for samples u ~ Gamma
  }

  return(result)
}

#' Update S for VB FC Template ICA
#'
#' @param mu_tau2,mu_A,E_AtA, Most recent posterior estimates
#' @param D_inv,D_inv_S Some pre-computed quantities.
#' @param BOLD (\eqn{T \times V} matrix) preprocessed fMRI data.
#' @param nICs,nvox Number of ICs, and the number of data locations.
#' @param final If TRUE, return cov_S. Default: \code{FALSE}
#'
#' @return List of length three: \code{mu_S} (QxV), \code{E_SSt} (QxQ), and \code{cov_S} (QxQxV).
#'
#' @keywords internal
update_S <- function(
  mu_tau2, mu_A, E_AtA,
  D_inv, D_inv_S,
  BOLD,
  nICs, nvox,
  final = FALSE){

  tmp1 <- (1/mu_tau2) * E_AtA

  cov_S <- array(0, dim = c(nICs, nICs, nvox))
  for(v in 1:nvox){
    cov_S[,,v] <- solve(tmp1 + diag(D_inv[v,]))
  }

  mu_S <- array(0, dim = c(nICs, nvox)) #QxV
  tmp2 <- (1/mu_tau2) * (t(mu_A) %*% BOLD)
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
#'
#' @return List of length two: \code{alpha1} and \code{beta1}.
#'
#' @keywords internal
update_tau2 <- function(
  BOLD, BOLD2, mu_A, E_AtA, mu_S, E_SSt, prior_params, ntime, nvox){

  alpha0 <- prior_params[1]
  beta0 <- prior_params[2]

  alpha1_nvox <- alpha0/nvox + ntime/2 #multiplied by nvox for computational stability
  beta1_nvox <- beta0/nvox +
    (1/2)*BOLD2/nvox -
    sum(BOLD * mu_A %*% mu_S/nvox) +
    (1/2)*sum(diag(E_AtA %*% E_SSt)/nvox)

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

##' Fast version of Matrix :: .bdiag() -- for the case of *many*  (k x k) matrices:
##' @param lmat list(<mat1>, <mat2>, ....., <mat_N>)  where each mat_j is a  k x k 'matrix'
##' @return a sparse (N*k x N*k) matrix of class  \code{"\linkS4class{dgCMatrix}"}.
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

##' Fast version of Matrix :: .bdiag() -- for the case of *many identical*  (k x k) matrices:
##' @param mat a  k x k 'matrix'
##' @param N how many times to repeat \code{mat}
##' @return a sparse (N*k x N*k) matrix of class  \code{"\linkS4class{dgCMatrix}"}.
bdiag_m2 <- function(mat, N) {
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

