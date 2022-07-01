#' @name EM_FCtemplateICA
#' @rdname EM_FCtemplateICA
#'
#' @title EM Algorithm for FC Template ICA Model
#'
#' @param template_mean (\eqn{V \times Q} matrix) mean maps for each IC in template,
#'  where \eqn{Q} is the number of ICs, \eqn{V=nvox} is the number of data locations.
#' @param template_var  (\eqn{V \times Q} matrix) between-subject variance maps for each IC in template
#' @param template_FC (list) Parameters of functional connectivity template
#' @param prior_params Alpha and beta parameters of IG prior on tau^2 (error variance)
#' @param BOLD (\eqn{V \times T} matrix) preprocessed fMRI data
#' @param AS_0 (list) initial guess at latent variables: A (\eqn{TxQ} mixing matrix),
#'  and S (\eqn{QxV} matrix of spatial ICs)
#' @param maxiter Maximum number of EM iterations. Default: 100.
#' @param epsilon Smallest proportion change in parameter estimates between iterations. Default: 0.01.
#' @param verbose If \code{TRUE}, display progress of algorithm. Default: \code{FALSE}.
#'
#' @return  A list:
#' theta (list of final parameter estimates),
#' subICmean (estimates of subject-level ICs),
#' subICvar (variance of subject-level ICs),
#' mixing_mean (estimates of subject-level mixing matrix),
#' mixing_var (variance of subject-level mixing matrix),
#' success (flag indicating convergence (\code{TRUE}) or not (\code{FALSE}))
#'
#' @details \code{EM_FCtemplateICA} implements the expectation-maximization
#'  (EM) algorithm for the functional connectivity (FC) template ICA model
#'
#' @keywords internal
#'
EM_FCtemplateICA <- function(template_mean,
                             template_var,
                             template_FC,
                             prior_params=c(0.001, 0.001),
                             BOLD,
                             AS_0,
                             maxiter=100,
                             epsilon=0.01,
                             verbose){

  #get initial values for A and S with dual regression - DONE
  #initialize the posterior expectations of A and S using those values
  #initialize theta using those posterior expectation -- this is the first update of theta using UpdateTheta_FCtemplateICA
  #UpdateTheta_FCtemplateICA will take in as an argument the current posterior moments of A and S
  #within the while loop, first run UpdateTheta_FCtemplateICA, then run Gibbs_AS_posterior to get required moments

  #in Gibbs_AS_posterior, have an argument to determine whether you want the final posterior mean/variance or the elements necessary for the M-step
  #Ani pointed out that we don't want to save all of the samples for computational reasons.  We can compute sums as we go.  Can compute sums in chunks of samples of 100, say, then decide how many chunks to drop before calculating the final means.

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

  #compute initial estimates of posterior moments
  A <- A_init <- AS_0$A #TxQ
  S <- S_init <- AS_0$S #QxV
  AS_init <- A_init %*% S_init #TxV -- for the terms E[a_t's_v]
  AS_sq_init <- (AS_init^2) #TxV -- for the terms E[(a_t's_v)^2]
  # These next are empirical estimates of the hyperparameters that we are estimating with the EM
  tau_v_init <- apply(BOLD - t(A %*% S),1,var)
  alpha_init <- apply(A,2,mean)
  sigma2_alpha <- max(apply(A,2,var))
  G_init <- cov(A)
  theta_new <- list(tau_v_init, alpha_init, G_init)

  #pre-compute sums over t
  A_sum_init <- colSums(A_init) #vector length Q
  AS_sq_sum_init <- colSums(AS_sq_init) #vector length V
  yAS_init <- t(BOLD) * AS_init #element-wise product -- for the terms y_tv * E[a_t's_v]
  yAS_sum_init <- colSums(yAS_init) #vector length V
  AAt_sum_init <- matrix(0, nICs,nICs) #QxQ matrix
  for(t in 1:ntime){ AAt_sum_init <- AAt_sum_init + tcrossprod(A_init[t,]) }
  post_sums <- list(A = A_sum_init, #only actually need sum over t
                    AAt = AAt_sum_init,  #only actually need sum over t
                    yAS = yAS_sum_init,  #only actually need sum over t of Y*AS (element-wise product)
                    AS_sq = AS_sq_sum_init) #only actually need sum over t

  err <- 1000 #large initial value for difference between iterations
  while(err > epsilon){

    if(verbose) cat(paste0(' ~~~~~~~~~~~~~~~~~~~~~ ITERATION ', iter, ' ~~~~~~~~~~~~~~~~~~~~~ \n'))
    theta_old <- theta_new
    t00 <- Sys.time()
    ### TO DO: RUN GIBBS SAMPLER TO SAMPLE FROM (A,S) AND UPDATE POSTERIOR_MOMENTS (RETURN SUMS OVER t=1,...,ntime as above)
    # tricolon <- NULL; Gibbs_AS_posterior <- function(x, ...){NULL} # Damon added this to avoid warnings.
    # post_sums <- Gibbs_AS_posterior(tricolon, final=FALSE)
    post_sums <-
      Gibbs_AS_posteriorCPP(
        nsamp = 1000,
        nburn = 0,
        template_mean = template_mean,
        template_var = template_var,
        S = t(S),
        A = A,
        G = theta_old[[3]],
        tau_v = theta_old[[1]],
        Y = BOLD,
        alpha = theta_old[[2]],
        final = F
      )

    post_sums <-
      Gibbs_AS_posterior(
        nsamp = 10,
        nburn = 0,
        template_mean = template_mean,
        template_var = template_var,
        S = S,
        A = A,
        G = theta_old[[3]],
        tau_v = theta_old[[1]],
        Y = BOLD,
        alpha = theta_old[[2]],
        final = F
      )

    #this function returns a list of tau_sq, alpha, G
    # This is the M-step. It might be better to perform the E-step first, as the
    # M-step assumes that we have a good estimate of the first level of the posterior
    # theta_new <- UpdateTheta_FCtemplateICA(template_mean,
    #                                       template_var,
    #                                       template_FC,
    #                                       prior_params,
    #                                       BOLD,
    #                                       post_sums,
    #                                       verbose=verbose)
    theta_new <-
      UpdateTheta_FCtemplateICAcpp(
        template_mean,
        template_var,
        template_FC,
        theta_old[[3]],
        prior_params,
        BOLD,
        post_sums,
        sigma2_alpha,
        TRUE
      )
    if(verbose) print(Sys.time() - t00)

    ### Compute change in parameters

    G_old <- theta_old[[3]]
    G_new <- theta_new[[3]]
    #2-norm <- largest eigenvalue <- sqrt of largest eigenvalue of AA'
    G_change <- norm(as.vector(G_new - G_old), type="2")/norm(as.vector(G_old), type="2")

    tau_old <- mean(theta_old[[1]])
    tau_new <- mean(theta_new[[1]])
    #2-norm <- largest eigenvalue <- sqrt of largest eigenvalue of AA'
    G_change <- norm(as.vector(G_new - G_old), type="2")/norm(as.vector(G_old), type="2")

    #
    # alpha_old <- theta_old[[2]]
    # alpha_new <- theta_new[[2]]
    # alpha_change <- abs(alpha_new - alpha_old) #avoid dividing by zero
  #
  #   change <- c(A_change, nu0_sq_change)
  #   err <- max(change)
  #   change <- format(change, digits=3, nsmall=3)
  #   if(verbose) cat(paste0('Iteration ',iter, ': Difference is ',change[1],' for A, ',change[2],' for nu0_sq \n'))
  #
  #   ### Move to next iteration
  #   theta <- theta_new
  #   iter <- iter + 1
  #   if(iter > maxiter){
  #     success <- 0
  #     warning(paste0('Failed to converge within ', maxiter,' iterations'))
  #     break() #exit loop
  #   }
  }
  #
  ### Compute final posterior means and variances of A and S

  post_AS <- Gibbs_AS_posterior(tricolon, final=TRUE)

  #
  # result <- list(subjICmean=miu_s,
  #                subjICvar=var_s,
  #                theta_MLE=theta,
  #                success_flag=success,
  #                error=err,
  #                numiter=iter-1,
  #                template_mean = template_mean,
  #                template_var = template_var)
  # return(result)
}


#' @rdname UpdateTheta_FCtemplateICA
#Update theta = (tau_sq, alpha, G)
UpdateTheta_FCtemplateICA <- function(template_mean,
                                      template_var,
                                      template_FC,
                                      prior_params,
                                      BOLD,
                                      post_sums,
                                      verbose){


  nvox <- nrow(BOLD)
  ntime <- ncol(BOLD)
  nICs <- ncol(template_mean)

  ### UPDATE TAU^2 (ERROR VARIANCE)

  tau_sq_new <- rep(NA, nvox)
  alpha_tau <- prior_params[1]
  beta_tau <- prior_params[2]
  #TO DO: Make this more efficient (no loop)
  for(v in 1:nvox){
    tau_sq_new[v] <- 1/(ntime + 2*alpha_tau + 2) * (post_sums$AS_sq[v] - 2*post_sums$yAS[v] + 2*beta_tau)
  }

  ### UPDATE ALPHA (TEMPORAL INTERCEPT)

  alpha_new <- 1/ntime*(post_sums$A)


  ### UPDATE G (TEMPORAL COVARIANCE, I.E. FUNCTIONAL CONNECTIVITY)

  nu0 <- template_FC$nu
  psi0 <- template_FC$psi
  alpha_alpha_t <- tcrossprod(alpha_new)
  tmp <- 2*tcrossprod(post_sums$A, alpha_new) - post_sums$AAt - ntime*alpha_alpha_t - psi0
  G_new <- (ntime + nu0 + nICs + 1)*solve(tmp)

  # RETURN NEW PARAMETER ESTIMATES
  theta_new <- list(
    tau_sq <- tau_sq_new,
    alpha <- alpha_new,
    G <- G_new
  )
  return(theta_new)
}

#' Find expectations from the posterior of A and S using MCMC
#'
#' @param nsamp number of usable samples
#' @param nburn number of burn-in samples
#' @param template_mean V by Q matrix of template mean values
#' @param template_var V by Q matrix of template variance values
#' @param S V by Q matrix of starting values for ICs
#' @param A T by Q matrix of starting values for mixing
#' @param G Q by Q prior covariance of A
#' @param tau_v vector of length V of observational variance for each data location
#' @param Y V by T matrix of observed BOLD values
#' @param alpha length Q vector prior mean of A
#' @param final (logical) should this output samples? Otherwise, summaries are output
#'
#' @return a list of posterior summaries or a list of MCMC samples
#' @export
Gibbs_AS_posterior <- function(nsamp = 1000,
                               nburn = 100,
                               template_mean,
                               template_var,
                               S,
                               A,
                               G,
                               tau_v,
                               Y,
                               alpha,
                               final = FALSE) {
  # number_samples = 1000
  niter <- nsamp + nburn
  Q = ncol(A)
  V = length(tau_v)
  if(ncol(S) == V) S <- t(S) # This is to make sure the dimensions are correct
  ntime <- ncol(Y)
  A_sum = numeric(Q)
  AAt_sum = matrix(0,Q,Q)
  yAS_sum = numeric(V)
  AS_sq = numeric(V)

  G_tau_inv <- Matrix::Diagonal(x = 1/tau_v)
  G_inv <- solve(G)
  alphaGinv <- as.numeric(t(alpha)%*%G_inv)
  YG <- t(Y) %*% G_tau_inv


  ### A is T by Q
  start_time <- proc.time()[3]
  for(i in 1:niter){
    #### update A
    # sigma_a = solve(t(S) %*% G_tau_inv %*% S + solve(G))
    sigma_a_inv <- as.matrix(t(S) %*% G_tau_inv %*% S + G_inv)
    chol_sigma_a_inv <- chol(sigma_a_inv)
    # sigma_a <- solve(sigma_a_inv)
    # cat(sigma_a,'\n')
    YGS <- YG %*% S
    YGSalphaGinv <- apply(YGS,1,function(ygs) ygs + alphaGinv)
    mu_a <- as.matrix(Matrix::solve(sigma_a_inv,YGSalphaGinv))
    cl <- parallel::makeCluster(parallel::detectCores() - 2)
    mu_at <- mu_a[,1]
    A <- t(parallel::parApply(cl,mu_a,2, function(mu_at, chol_sigma_a_inv){
      Qa <- length(mu_at)
      Za <- rnorm(Qa)
      at <- mu_at + as.numeric(Matrix::solve(chol_sigma_a_inv,Za))
      return(at)
    }, chol_sigma_a_inv = chol_sigma_a_inv))
    parallel::stopCluster(cl)
    # cat(YGS[ntime,],"\n")
    # cat(YGS[ntime,] + alphaGinv,"\n")
    # cat("mu_a1:",sigma_a %*% t(YGS[1,] + alphaGinv), "\n")
    # mu_a <- apply(YGS,1,function(ygs) tcrossprod(sigma_a, ygs + alphaGinv))
    # cat(mu_a,"\n")
    # A = mvtnorm::rmvnorm(n = 1, mean = mu_a, sigma = sigma_a) # This doesn't work because mu_a is a matrix
    # A = t(apply(mu_a, 2, function(ma) mvtnorm::rmvnorm(n = 1,mean = ma,sigma = sigma_a)))

    # G_tauv_inv = diag(1/template.var, nrow = T) # Don't actually need to make this
    # G_sv_inv = Matrix::Diagonal(x = 1/sigma_s) ### Need var names for these

    #### update S
    for(v in 1:V){
      sigma_s = solve((1/tau_v[v]) * crossprod(A) + Matrix::Diagonal(x = template_var[v,]))
      AtYvtempVarMean <- (1/tau_v[v]) * crossprod(A,Y[v,]) + Matrix::Diagonal(x = template_var[v,])%*%template_mean[v,]
      mu = sigma_s %*% AtYvtempVarMean
      # if(v == 1) cat("mu_s1:",mu@x,"\n")
      # mu_s = cbind(mu_s, mu)
      S[,v] <- mvtnorm::rmvnorm(n = 1, mean = mu, sigma = sigma_s)
    }

    # S = mvtnorm::rmvnorm(n = 1, mean = mu_s, sigma = sigma_s)

    #### I'm not sure about the definitions of AAt and AS_sq, so the following lines need to be checked.
    if(i > nburn) {
      AtS <- A %*% S
      A_sum = A_sum + colSums(A)
      AAt_sum = AAt_sum + colSums(A %*% t(A))
      yAS_sum = yAS_sum + colSums(t(Y) * AtS)
      AS_sq = AS_sq + apply(AtS^2,2,sum)
    }

    if(niter >= 10 & i %% round(niter / 10) == 0) {
      cat("Posterior sample",i,"of",niter,". Estimated time remaining:", (proc.time()[3] - start_time) / i * (niter - i),"\n")
    }
  }
  #### return estimates of A = A_sum_init, #only actually need sum over t
  #AAt = AAt_sum_init,  #only actually need sum over t
  #yAS = yAS_sum_init,  #only actually need sum over t of Y*AS (element-wise product)
  #AS_sq = AS_sq_sum_init

  A_sum = A_sum/nsamp
  AAt_sum = AAt_sum/nsamp
  yAS_sum = yAS_sum/nsamp
  AS_sq = AS_sq/nsamp
  total_time <- proc.time()[3] - start_time

  return(list(A_sum = A_sum, AAt_sum = AAt_sum, yAS_sum = yAS_sum, AS_sq = AS_sq, total_time = total_time))

}
