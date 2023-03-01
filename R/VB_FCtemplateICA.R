#' @name VB_FCtemplateICA
#' @rdname VB_FCtemplateICA
#'
#' @title VB Algorithm for FC Template ICA Model
#'
#' @param template_mean (\eqn{V \times Q} matrix) mean maps for each IC in template,
#'  where \eqn{Q} is the number of ICs, \eqn{V=nvox} is the number of data locations.
#' @param template_var  (\eqn{V \times Q} matrix) between-subject variance maps for each IC in template
#' @param template_FC (list) Parameters of functional connectivity template
#' @param prior_params Alpha and beta parameters of IG prior on tau^2 (error variance)
#' @param BOLD (\eqn{V \times T} matrix) preprocessed fMRI data
#' @param AS_0 (list) initial guess at latent variables: A (\eqn{TxQ} mixing matrix),
#'  and S (\eqn{QxV} matrix of spatial ICs)
#' @param maxiter Maximum number of VB iterations. Default: 100.
#' @param epsilon Smallest proportion change in parameter estimates between iterations. Default: 0.01.
#' @param verbose If \code{TRUE}, display progress of algorithm. Default: \code{FALSE}.
#'
#'
VB_FCtemplateICA <- function(template_mean, #VxQ
                              template_var, #VxQ
                              template_FC,
                              prior_params=c(0.001, 0.001),
                              BOLD, #VxT
                              A0, S0, S0_var,
                              maxiter=100,
                              epsilon=0.01,
                              verbose=FALSE){

  if(!all.equal(dim(template_var), dim(template_mean))) stop('The dimensions of template_mean and template_var must match.')
  ntime <- ncol(BOLD) #length of timeseries
  nvox <- nrow(BOLD) #number of brain locations
  if(ntime > nvox) warning('More time points than brain locations. Are you sure?')
  if(nrow(template_mean) != nvox) stop('Templates and BOLD must have the same number of brain locations (columns).')

  nICs <- nrow(template_FC$psi)   #number of ICs
  if(ncol(template_mean) != nICs) stop('template_FC is incompatible with template_mean & template_var. Check number of ICs in each.')
  if(nICs > nvox) stop('Cannot estimate more ICs than brain locations.')
  if(nICs > ntime) stop('Cannot estimate more ICs than time points.')

  iter <- 1
  success <- 1
  template_var[template_var < 1e-6] <- 1e-6 #to prevent problems when inverting covariance (this should rarely/never happen with NN variance)

  #1. Compute initial estimates of posterior moments for G, A, S, V(S), tau^2

  mu_A <- A0 #TxQ
  mu_S <- S0 #QxV
  cov_S <- array(0, dim = c(nICs, nICs, nvox)) #QxQxV
  for(v in 1:nvox){cov_S[,,nvox] <- diag(S0_var[,v]) }
  mu_tau2 <- apply(BOLD - t(A %*% S),1,var) #Vx1
  mu_alpha <- colMeans(A0)
  mu_G <- cov(A0)
  cov_alpha <- 1/ntime * mu_G
  BOLD <- t(BOLD) #make the BOLD TxV to match the paper

  #2. Iteratively update approximate posteriors

  err <- 1000 #large initial value for difference between iterations
  while(err > epsilon){

    if(verbose) cat(paste0(' ~~~~~~~~~~~~~~~~~~~~~ VB ITERATION ', iter, ' ~~~~~~~~~~~~~~~~~~~~~ \n'))
    t00 <- Sys.time()

    #a. UPDATE A

    A_new <- update_A(mu_tau2, mu_S, cov_S, mu_G, mu_alpha, BOLD, ntime, nICs, nvox)
    change_A <- mean((A_new$mu_A - mu_A)^2)
    mu_A <- A_new$mu_A
    cov_A <- A_new$cov_A

    #b. UPDATE S
    S_new <- update_S(mu_tau2, mu_A, cov_A, template_var, template_mean, BOLD, ntime, nICs, nvox)
    change_S <- mean((S_new$mu_S - mu_S)^2)
    mu_S <- S_new$mu_S
    cov_S <- S_new$cov_S

    #c1. UPDATE alpha
    mu_alpha <- colMeans(mu_A)
    cov_alpha <- mu_G/ntime

    #c2. UPDATE G
    G_new <- update_G(mu_A, cov_A, mu_alpha, cov_alpha, template_FC, ntime, nICs)
    mu_G_new <- G_new$psi1 / G_new$nu1
    change_G <- mean((mu_G - mu_G_new)^2)
    mu_G <- mu_G_new

    #d. UPDATE tau^2
    tau2_new <- update_tau2(BOLD, mu_A, mu_S, cov_A, cov_S, prior_params, ntime, nvox)
    mu_tau2_new <- tau2_new$beta1/(tau2_new$alpha2 - 1)
    change_tau2 <- mean(mu_tau2 - mu_tau2_new)^2
    mu_tau2 <- mu_tau2_new

    if(verbose) print(Sys.time() - t00)

    change <- c(change_A, change_S, change_G, change_tau2)
    err <- max(change)
    change <- format(change, digits=3, nsmall=3)
    if(verbose) cat(paste0('Iteration ',iter, ': Difference is ',change[1],' for A, ',change[2],' for S, ',change[3],' for G, ',change[4],' for tau2 \n'))

    ### Move to next iteration
    iter <- iter + 1
    if(iter > maxiter){
      success <- 0
      warning(paste0('Failed to converge within ', maxiter,' iterations'))
      break() #exit loop
    }
  } #end iterations

  result <- list(S_mean = mu_S,
                 S_cov = cov_S,
                 A_mean = mu_A,
                 A_cov = cov_A,
                 G_mean = mu_G,
                 success_flag=success,
                 error=err,
                 numiter=iter-1,
                 template_mean,
                 template_var,
                 template_FC)

  return(result)



}

update_A <- function(mu_tau2, mu_S, cov_S, mu_G, mu_alpha, BOLD, ntime, nICs, nvox){

  G_inv <- solve(mu_G)
  cov_A <- array(0, dim = c(ntime, nICs, nICs)) #TxQxQ
  for(t in 1:ntime){
    #sum over v=1...V part
    tmp <- matrix(0, nICs, nICs)
    for(v in 1:nvox){
      tmp <- tmp + (1/mu_tau2[v]) * (tcrossprod(mu_S[,v]) + cov_S[,,v])
    }
    cov_A[t,,] <- solve(tmp + G_inv)
  }

  mu_A <- array(0, dim = c(ntime, nICs)) #TxQ
  for(t in 1:ntime){
    #sum over v=1...V part
    tmp <- rep(0, nICs)
    for(v in 1:nvox){ # TO DO: do this with matrix operations
      tmp <- tmp + (1/mu_tau2[v]) * BOLD[t,v] * mu_S[,v]
    }
    mu_A[t,] <- cov_A[t,,] %*% (tmp + G_inv %*% mu_alpha)
  }

  return(list(mu_A, cov_A))

}


update_S <- function(mu_tau2, mu_A, cov_A, template_var, template_mean, BOLD, ntime, nICs, nvox){
  cov_S <- array(0, dim = c(nICs, nICs, nvox))
  for(v in 1:nvox){
    #sum over t=1...T part
    tmp <- matrix(0, nICs, nICs)
    for(t in 1:ntime){
      tmp <- tmp + (tcrossprod(mu_A[t,]) + cov_A[t,,])
    }
    D_v_inv <- diag(1/template_var[v,])
    cov_S[,,v] <- solve((1/mu_tau2[v]) * tmp + D_v_inv)
  }

  mu_S <- array(0, dim = c(nICs, nvox)) #QxV
  for(v in 1:nvox){
    #sum over t=1...T part
    tmp <- rep(0, nICs)
    for(t in 1:ntime){ # TO DO: do this with matrix operations
      tmp <- tmp + BOLD[t,v] * mu_A[t,]
    }
    mu_S[,v] <- cov_S[,,v] %*% ( (1/mu_tau2[v]) * tmp + D_v_inv %*% template_mean[v,] )
  }

  return(list(mu_S, cov_S))
}

update_G <- function(mu_A, cov_A, mu_alpha, cov_alpha, template_FC, ntime, nICs){

  nu0 <- template_FC$nu
  psi0 <- template_FC$psi

  nu1 <- nu0 + ntime
  #sums over t=1,...,T in Psi_G
  tmp1 <- tmp2 <- matrix(0, nICs, nICs)
  for(t in 1:ntime){
    tmp1 <- tmp1 + tcrossprod(mu_A[t,]) + cov_A[t,,]
    tmp2 <- tmp2 + mu_A[t,] #this is just a colSums or rowSums
  }
  psi1 <- psi0 +
    tmp1 -
    2*tcrossprod(mu_alpha, tmp2) +
    ntime * (tcrossprod(mu_alpha) + cov_alpha)

  return(list(nu1, psi1))

}

update_tau2 <- function(BOLD, mu_A, mu_S, cov_A, cov_S, prior_params, ntime, nvox){

  alpha0 <- prior_params[1]
  beta0 <- prior_params[2]

  alpha1 <- alpha0 + ntime/2
  beta1 <- rep(0, nvox)

  #pre-compute sum over t inside of trace (it doesn't involve v)
  tmp3a <- 0
  for(t in 1:ntime){ tmp3a <- tmp3a + tcrossprod(mu_A[t,]) + cov_A[t,,] }

  for(v in 1:nvox){
    #sums over t=1,...,T in beta_v
    tmp1 <- colSums(BOLD^2) #sum over rows t=1,...,T
    tmp2 <- 0
    for(t in 1:ntime){
      tmp2 <- tmp2 + BOLD[t,v] * mu_A[t,] #TO DO: replace with matrix operations
    }
    tmp3b <- tcrossprod(mu_S[,v]) + cov_S[,,v]
    Tr_v <- sum(diag(tmp3a %*% tmp3b)) #TO DO: make trace more efficient with element-wise multiplication and colSums ?
    beta1[v] <- beta0 + (1/2)*tmp1 - tcrossprod(tmp2, mu_S[,v]) + (1/2)*Tr_v
  }

  return(list(alpha1, beta1))
}


