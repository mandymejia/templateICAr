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

  mu_A <- scale(A0, center=FALSE) #TxQ
  mu_S <- t(S0) #QxV
  cov_S <- array(0, dim = c(nICs, nICs, nvox)) #QxQxV
  for (v in 1:nvox) { cov_S[,,v] <- diag(S0_var[v,]) }
  mu_tau2 <- apply(BOLD - t(mu_A %*% mu_S),1,var) #Vx1
  mu_G <- cov(A0)
  BOLD <- t(BOLD) #make the BOLD TxV to match the paper

  #pre-compute some stuff
  D_inv <- 1/template_var
  D_inv_S <- D_inv * template_mean
  BOLD2 <- sum(BOLD^2) #sum over all v,t

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
  while (err > epsilon | iter <= miniter) {

    if(verbose) cat(paste0(' ~~~~~~~~~~~~~~~~~~~~~ VB ITERATION ', iter, ' ~~~~~~~~~~~~~~~~~~~~~ \n'))
    t00 <- Sys.time()

    #a. UPDATE A

    A_new <- update_A(mu_tau2, mu_S, cov_S, mu_G, #mu_alpha,
                      BOLD, ntime, nICs, nvox)
    mu_A_old <- mu_A
    mu_A <- A_new[[1]]
    cov_A <- A_new[[2]]
    #constrain each column of A to have var=1
    sd_A <- apply(mu_A, 2, sd)
    mu_A <- scale(mu_A, center=FALSE)
    cov_A <- diag(1/sd_A) %*% cov_A %*% diag(1/sd_A)
    A_new <- list(mu_A=mu_A, cov_A=cov_A)
    change_A <- mean(abs(mu_A - mu_A_old)/(abs(mu_A_old)+0.1)) #add 0.1 to denominator to avoid dividing by zero
    #change_A <- sqrt(crossprod(c(mu_A - mu_A_old))) #same as in SQUAREM

    #b. UPDATE S
    S_new <- update_S(mu_tau2, mu_A, cov_A, D_inv, D_inv_S, BOLD, ntime, nICs, nvox)
    change_S <- mean(abs(S_new[[1]] - mu_S)/(abs(mu_S)+0.1)) #add 0.1 to denominator to avoid dividing by zero
    #change_S <- sqrt(crossprod(c(S_new[[1]] - mu_S))) #same as in SQUAREM
    mu_S <- S_new[[1]]
    cov_S <- S_new[[2]]

    #c. UPDATE G
    if(method_FC == 'VB1'){
      G_new <- update_G_IW(mu_A, cov_A, #mu_alpha, cov_alpha,
                        template_FC, ntime, nICs)
      psi_G <- G_new[[2]]
      nu_G <- G_new[[1]]
      mu_G_new <- psi_G / (nu_G - nICs - 1)
      change_G <- mean(abs(mu_G - mu_G_new)/(abs(mu_G)+0.1)) #add 0.1 to denominator to avoid dividing by zero
      #change_G <- sqrt(crossprod(c(mu_G_new - mu_G))) #same as in SQUAREM
      mu_G <- mu_G_new
    } else {

    # [TO DO] implement Cholesky prior-based estimation of G

    }

    #d. UPDATE tau^2
    tau2_new <- update_tau2(BOLD, BOLD2, mu_A, mu_S,
                            cov_A = ntime*cov_A, #sum over t
                            cov_S = apply(cov_S, 1:2, sum), #sum over v
                            prior_params, ntime, nvox)
    beta1 <- tau2_new[[2]] #scalar
    alpha1 <- tau2_new[[1]] #scalar
    mu_tau2_new <- beta1/(alpha1 - 1)
    change_tau2 <- abs(mu_tau2_new - mu_tau2)/(mu_tau2+0.1) #add 0.1 to denominator to avoid dividing by zero
    #change_tau2 <- sqrt(crossprod(c(mu_tau2_new - mu_tau2))) #same as in SQUAREM
    mu_tau2 <- mu_tau2_new

    if (verbose) print(Sys.time() - t00)

    #change in estimates
    change <- c(change_A, change_S, change_G, change_tau2)
    err <- max(change)
    change <- format(change, digits=3, nsmall=3)
    if(verbose) cat(paste0('Iteration ',iter, ':l0 Difference is ',change[1],' for A, ',change[2],' for S, ',change[3],' for G, ',change[4],' for tau2 \n'))

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
          results_inter[[which_eps]] <- list(S = t(mu_S), A = mu_A, G_mean = mu_G, tau2_mean = mu_tau2_new, error=err, numiter=iter)
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

  cov_S_list <- lapply(seq(dim(cov_S)[3]), function(x) cov_S[ , , x])
  subjICse <- sqrt(sapply(cov_S_list, diag))

  list(
    subjICmean = t(mu_S),
    subjICse = t(subjICse),
    S_cov = cov_S,
    A = mu_A,
    A_cov = cov_A,
    G_mean = mu_G,
    #alpha_mean = mu_alpha,
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
#' @param mu_tau2,mu_S,cov_S,mu_G Most recent estimates of posterior
#' moments for these variables.
#' @param BOLD (\eqn{T \times V} matrix) preprocessed fMRI data.
#' @param ntime,nICs,nvox Number of timepoints in data, number of ICs, and
#'  the number of data locations.
#'
#' @return List of length two: \code{mu_A} (TxQ) and \code{cov_A} (QxQ).
#'
#' @keywords internal
update_A <- function(
  mu_tau2, mu_S, cov_S, mu_G,
  BOLD,
  ntime, nICs, nvox){

  #cov_A (QxQ) -- common across t=1,...,T
  G_inv <- solve(mu_G)
  tmp1 <- mu_S/sqrt(mu_tau2)
  tmp1 <- tcrossprod(tmp1)
  tmp2 <- cov_S/mu_tau2
  tmp2 <- apply(tmp2, c(1,2), sum)
  tmp <- tmp1 + tmp2
  cov_A <- solve(tmp + G_inv)

  #mu_A
  tmp1 <- BOLD/mu_tau2
  tmp2 <- tmp1 %*% t(mu_S)
  mu_A <- array(0, dim = c(ntime, nICs)) #TxQ
  for(t in 1:ntime) mu_A[t,] <- cov_A %*% tmp2[t,]

  list(mu_A=mu_A, cov_A=cov_A)
}

#' Update S for VB FC Template ICA
#'
#' @param mu_tau2,mu_A,cov_A, Most recent estimates of posterior moments for
#'  these variables.
#' @param D_inv,D_inv_S Some pre-computed quantities.
#' #' @param BOLD (\eqn{T \times V} matrix) preprocessed fMRI data.
#' @param ntime,nICs,nvox Number of timepoints in data, number of ICs, and
#'  the number of data locations.
#'
#' @return List of length two: \code{mu_S} (QxV) and \code{cov_S} (QxQxV).
#'
#' @keywords internal
update_S <- function(
  mu_tau2, mu_A, cov_A,
  D_inv, D_inv_S,
  BOLD,
  ntime, nICs, nvox){

  cov_S <- array(0, dim = c(nICs, nICs, nvox))
  #sum over t=1...T part
  tmp <- t(mu_A) %*% mu_A + ntime * cov_A
  for(v in 1:nvox){
    cov_S[,,v] <- solve((1/mu_tau2) * tmp + diag(D_inv[v,]))
  }

  mu_S <- array(0, dim = c(nICs, nvox)) #QxV
  tmp <- t(BOLD) %*% mu_A
  for(v in 1:nvox){
    #sum over t=1...T part
    #D_v_inv <- diag(1/template_var[v,])
    mu_S[,v] <- cov_S[,,v] %*% ( (1/mu_tau2) * tmp[v,] + D_inv_S[v,] )
  }

  list(mu_S=mu_S, cov_S=cov_S)
}

#' Update G for VB FC Template ICA
#'
#' @param mu_A,cov_A Most recent estimates of posterior
#' moments for these variables.
#' @param template_FC (list) Parameters of functional connectivity template.
#' @param ntime,nICs Number of timepoints in data and number of ICs.
#'
#' @return List of length two: \code{nu1} and \code{psi1}.
#'
#' @keywords internal
update_G_IW <- function(
  mu_A, cov_A, template_FC, ntime, nICs){

  nu0 <- template_FC$nu
  psi0 <- template_FC$psi

  nu1 <- nu0 + ntime
  #sums over t=1,...,T in Psi_G
  tmp1 <- t(mu_A) %*% mu_A + ntime*cov_A
  #tmp2 <- colSums(mu_A)
  psi1 <- psi0 + tmp1 #-
    #2*tcrossprod(mu_alpha, tmp2) +
    #ntime * (tcrossprod(mu_alpha) + cov_alpha)

  list(nu1=nu1, psi1=psi1)
}

#' Update tau for VB FC Template ICA
#'
#' @param BOLD (\eqn{T \times V} matrix) preprocessed fMRI data.
#' @param BOLD2 A precomputed quantity, \code{sum(BOLD^2)}
#' @param mu_A (\eqn{T \times Q} matrix) Current estimate of A
#' @param mu_S (\eqn{Q \times V} matrix) Current estimate of S
#' @param cov_A (\eqn{Q \times Q} matrix) Current estimate of sum_t Cov(a_t)
#' @param cov_S (\eqn{Q \times Q} matrix) Current estimate of sum_v Cov(s_v)
#' @param prior_params Alpha and beta parameters of IG prior on \eqn{\tau^2}
#'  (error variance).
#' @param ntime,nvox Number of timepoints in data and the number of data
#'  locations.
#'
#' @return List of length two: \code{alpha1} and \code{beta1}.
#'
#' @keywords internal
update_tau2 <- function(
  BOLD, BOLD2, mu_A, mu_S, cov_A, cov_S, prior_params, ntime, nvox){

  alpha0 <- prior_params[1]
  beta0 <- prior_params[2]

  E_aat <- cov_A + t(mu_A) %*% mu_A # sum_t E[a_t a_t'] = (sum_t Cov(a_t)) + mu(A)'mu(A)
  E_sst <- cov_S + mu_S %*% t(mu_S)

  alpha1 <- alpha0 + ntime*nvox/2
  beta1 <- beta0 +
    (1/2)*BOLD2 +
    sum(BOLD * mu_A %*% mu_S) +
    sum(diag(E_aat %*% E_sst))


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

  list(alpha1=alpha1, beta1=beta1)
}

