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
#' @param prior_params Alpha and beta parameters of IG prior on \eqn{\tau^2}
#'  (error variance). Default: \code{0.001} for both.
#' @param BOLD (\eqn{V \times T} matrix) preprocessed fMRI data.
#' @param A0,S0,S0_var Initial guesses at latent variables: \code{A} (\eqn{TxQ}
#'  mixing matrix), \code{S} (\eqn{QxV} matrix of spatial ICs), and
#'  covariance matrix \code{S0_var}.
#' @param maxiter Maximum number of VB iterations. Default: \code{100}.
#' @param miniter Minimum number of VB iterations. Default: \code{5}.
#' @param epsilon Smallest proportion change in ELBO between iterations.
#' Default: \code{0.01} percent.
#' @param verbose If \code{TRUE}, display progress of algorithm.
#' Default: \code{FALSE}.
#'
#' @return Large list of results.
#'
#' @keywords internal
#'
VB_FCtemplateICA <- function(
  template_mean, #VxQ
  template_var, #VxQ
  template_FC,
  prior_params=c(0.001, 0.001),
  BOLD, #VxT
  A0, S0, S0_var,
  maxiter=100,
  miniter=5,
  epsilon=0.0001,
  verbose=FALSE){

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
  for (v in 1:nvox) { cov_S[,,nvox] <- diag(S0_var[v,]) }
  mu_tau2 <- apply(BOLD - t(mu_A %*% mu_S),1,var) #Vx1
  mu_alpha <- colMeans(A0)
  mu_G <- cov(A0)
  cov_alpha <- 1/ntime * mu_G
  BOLD <- t(BOLD) #make the BOLD TxV to match the paper

  #pre-compute some stuff
  D_inv <- 1/template_var
  D_inv_S <- D_inv * template_mean
  BOLD2_v <- colSums(BOLD^2) #sum over t=1,...,T

  #2. Iteratively update approximate posteriors

  err <- 1000 #large initial value for difference between iterations
  ELBO_vals <- rep(NA, maxiter)
  while (err > epsilon | iter <= miniter) {

    if(verbose) cat(paste0(' ~~~~~~~~~~~~~~~~~~~~~ VB ITERATION ', iter, ' ~~~~~~~~~~~~~~~~~~~~~ \n'))
    t00 <- Sys.time()

    #a. UPDATE A

    A_new <- update_A(mu_tau2, mu_S, cov_S, mu_G, mu_alpha, BOLD, ntime, nICs, nvox)
    change_A <- mean(abs(A_new[[1]] - mu_A)/(abs(mu_A)+0.1)) #add 0.1 to denominator to avoid dividing by zero
    mu_A <- A_new[[1]]
    cov_A <- A_new[[2]]

    #b. UPDATE S
    S_new <- update_S(mu_tau2, mu_A, cov_A, D_inv, D_inv_S, BOLD, ntime, nICs, nvox)
    change_S <- mean(abs(S_new[[1]] - mu_S)/(abs(mu_S)+0.1)) #add 0.1 to denominator to avoid dividing by zero
    mu_S <- S_new[[1]]
    cov_S <- S_new[[2]]

    #c1. UPDATE alpha
    mu_alpha <- colMeans(mu_A)
    cov_alpha <- mu_G/ntime

    #c2. UPDATE G
    G_new <- update_G(mu_A, cov_A, mu_alpha, cov_alpha, template_FC, ntime, nICs)
    psi_G <- G_new[[2]]
    nu_G <- G_new[[1]]
    mu_G_new <- psi_G / nu_G
    change_G <- mean(abs(mu_G - mu_G_new)/(abs(mu_G)+0.1)) #add 0.1 to denominator to avoid dividing by zero
    mu_G <- mu_G_new

    #d. UPDATE tau^2
    tau2_new <- update_tau2(BOLD, BOLD2_v, mu_A, mu_S, cov_A, cov_S, prior_params, ntime, nvox)
    beta1 <- tau2_new[[2]] #length V
    alpha1 <- tau2_new[[1]] #scalar
    mu_tau2_new <- beta1/(alpha1 - 1)
    change_tau2 <- mean(abs(mu_tau2_new - mu_tau2)/(mu_tau2+0.1)) #add 0.1 to denominator to avoid dividing by zero
    mu_tau2 <- mu_tau2_new

    if (verbose) print(Sys.time() - t00)

    #change in estimates
    change <- c(change_A, change_S, change_G, change_tau2)
    err <- max(change)
    change <- format(change, digits=3, nsmall=3)
    if(verbose) cat(paste0('Iteration ',iter, ': Proportional Difference is ',change[1],' for A, ',change[2],' for S, ',change[3],' for G, ',change[4],' for tau2 \n'))

    #ELBO
    ELBO_vals[iter] <- ELBO(mu_S, cov_S, cov_A, cov_alpha, template_mean, template_var, ntime)
    if(iter > 1) {
      err <- (ELBO_vals[iter] - ELBO_vals[iter - 1])/ELBO_vals[iter - 1]
      if(verbose) cat(paste0('Iteration ',iter, ': Proportional Change in ELBO is ',round(100*err, 3),'% \n'))
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
    success_flag=success,
    error=err,
    numiter=iter-1,
    template_mean,
    template_var,
    template_FC
  )
}

#' Update A for VB FC Template ICA
#'
#' @param mu_tau2,mu_S,cov_S,mu_G,mu_alpha Most recent estimates of posterior
#' moments for these variables.
#' @param BOLD (\eqn{T \times V} matrix) preprocessed fMRI data.
#' @param ntime,nICs,nvox Number of timepoints in data, number of ICs, and
#'  the number of data locations.
#'
#' @return List of length two: \code{mu_A} and \code{cov_A}.
#'
#' @keywords internal
update_A <- function(
  mu_tau2, mu_S, cov_S, mu_G, mu_alpha,
  BOLD,
  ntime, nICs, nvox){

  #cov_A (QxQ) -- common across t=1,...,T
  G_inv <- solve(mu_G)
  tmp1 <- mu_S/matrix(sqrt(mu_tau2), nrow=nICs, ncol=nvox, byrow = TRUE)
  tmp1 <- tcrossprod(tmp1)
  tmp2 <- cov_S/array(rep(mu_tau2, each=nICs*nICs), dim=dim(cov_S))
  tmp2 <- apply(tmp2, c(1,2), sum)
  tmp <- tmp1 + tmp2
  cov_A <- solve(tmp + G_inv)

  #mu_A
  tmp1 <- BOLD/matrix(mu_tau2, nrow=ntime, ncol=nvox, byrow = TRUE)
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
#' @return List of length two: \code{mu_S} and \code{cov_S}.
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
    cov_S[,,v] <- solve((1/mu_tau2[v]) * tmp + diag(D_inv[v,]))
  }

  mu_S <- array(0, dim = c(nICs, nvox)) #QxV
  tmp <- t(BOLD) %*% mu_A
  for(v in 1:nvox){
    #sum over t=1...T part
    #D_v_inv <- diag(1/template_var[v,])
    mu_S[,v] <- cov_S[,,v] %*% ( (1/mu_tau2[v]) * tmp[v,] + D_inv_S[v,] )
  }

  list(mu_S=mu_S, cov_S=cov_S)
}

#' Update G for VB FC Template ICA
#'
#' @param mu_A,cov_A,mu_alpha,cov_alpha Most recent estimates of posterior
#' moments for these variables.
#' @param template_FC (list) Parameters of functional connectivity template.
#' @param ntime,nICs Number of timepoints in data and number of ICs.
#'
#' @return List of length two: \code{nu1} and \code{psi1}.
#'
#' @keywords internal
update_G <- function(
  mu_A, cov_A, mu_alpha, cov_alpha, template_FC, ntime, nICs){

  nu0 <- template_FC$nu
  psi0 <- template_FC$psi

  nu1 <- nu0 + ntime
  #sums over t=1,...,T in Psi_G
  tmp1 <- t(mu_A) %*% mu_A + ntime*cov_A
  tmp2 <- colSums(mu_A)
  psi1 <- psi0 + tmp1 -
    2*tcrossprod(mu_alpha, tmp2) +
    ntime * (tcrossprod(mu_alpha) + cov_alpha)

  list(nu1=nu1, psi1=psi1)
}

#' Update tau for VB FC Template ICA
#'
#' @param BOLD (\eqn{T \times V} matrix) preprocessed fMRI data.
#' @param BOLD2_v A precomputed quantity: \code{colSums(BOLD^2)}, a numeric
#'  vector of length \eqn{V}.
#' @param mu_A,mu_S,cov_A,cov_S, Most recent estimates of posterior
#' moments for these variables.
#' @param prior_params Alpha and beta parameters of IG prior on \eqn{\tau^2}
#'  (error variance).
#' @param ntime,nvox Number of timepoints in data and the number of data
#'  locations.
#'
#' @return List of length two: \code{alpha1} and \code{beta1}.
#'
#' @keywords internal
update_tau2 <- function(
  BOLD, BOLD2_v, mu_A, mu_S, cov_A, cov_S, prior_params, ntime, nvox){

  alpha0 <- prior_params[1]
  beta0 <- prior_params[2]

  alpha1 <- alpha0 + ntime/2
  beta1 <- rep(0, nvox)

  #pre-compute sum over t inside of trace (it doesn't involve v)
  tmp2a <- t(mu_A) %*% mu_A + ntime*cov_A

  tmp1 <- t(BOLD) %*% mu_A
  for(v in 1:nvox){
    #sums over t=1,...,T in beta_v
    tmp1v <- crossprod(tmp1[v,], mu_S[,v])
    tmp2b <- tcrossprod(mu_S[,v]) + cov_S[,,v]
    Tr_v <- sum(diag(tmp2a %*% tmp2b)) #TO DO: make trace more efficient with element-wise multiplication and colSums ?
    beta1[v] <- beta0 + (1/2)*BOLD2_v[v] - tmp1v + (1/2)*Tr_v
  }

  list(alpha1=alpha1, beta1=beta1)
}

#' ELBO (Evidence Lower Bound)
#'
#' Computes the Evidence Lower Bound for the VB Algorithm for FC Template ICA
#'
#' @param mu_S (Q x V)
#' @param cov_S (Q x Q x V)
#' @param cov_A (Q x Q) (common across all a_t)
#' @param cov_alpha (Q x Q)
#' @param template_mean (VxQ) The template mean for S
#' @param template_var (VxQ) The template variance for S
#' @param ntime Number of timepoints in data
#'
#' @return ELBO
#'
#' @keywords internal
#' @importFrom mvtnorm dmvnorm
#' @importFrom stats dgamma
ELBO <- function(
  mu_S, cov_S,
  cov_A,
  cov_alpha,
  template_mean, template_var,
  ntime){

  #compute det(X)
  halflogdetX <- function(X){
    cholX <- chol(X)
    sum(log(diag(cholX)))
  }

  nvox <- ncol(mu_S)

  #first part of ELBO
  halflogdet_cov_S_v <- apply(cov_S, 3, halflogdetX)/(nvox*ntime)
  halflogdet_cov_A <- halflogdetX(cov_A)/(nvox*ntime) #scalar
  halflogdet_cov_alpha <- halflogdetX(cov_alpha)/(nvox*ntime) #scalar

  #second part of ELBO
  var_S <- apply(cov_S, 3, diag) #grab the diagonal elements of QxQ cov matrices -- QxV
  part2_v <- 1/2 * colSums((var_S + mu_S^2 - 2*t(template_mean)*mu_S)/t(template_var))/(nvox*ntime)

  ELBO <- sum(-halflogdet_cov_S_v + part2_v) - ntime*halflogdet_cov_A - halflogdet_cov_alpha
  ELBO
}





