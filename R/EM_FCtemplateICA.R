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
#' @importFrom expm sqrtm
#' @importFrom Matrix Diagonal
#' @importFrom matrixStats rowVars
#' @importFrom stats sd cor
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
  S <- t(S)

  #pre-compute sums over t
  A_sum_init <- colSums(A_init) #vector length Q
  AS_sq_sum_init <- colSums(AS_sq_init) #vector length V
  yAS_init <- t(BOLD) * AS_init #element-wise product -- for the terms y_tv * E[a_t's_v]
  yAS_sum_init <- colSums(yAS_init) #vector length V
  AtA_sum_init <- crossprod(A_init) #QxQ matrix
  post_sums <- list(A = A_sum_init,
                    AtA = AtA_sum_init,
                    yAS = yAS_sum_init,
                    AS_sq = AS_sq_sum_init)

  err <- 1000 #large initial value for difference between iterations
  while(err > epsilon){

    if(verbose) cat(paste0(' ~~~~~~~~~~~~~~~~~~~~~ EM ITERATION ', iter, ' ~~~~~~~~~~~~~~~~~~~~~ \n'))
    theta_old <- theta_new
    t00 <- Sys.time()
    ### TO DO: RUN GIBBS SAMPLER TO SAMPLE FROM (A,S) AND UPDATE POSTERIOR_MOMENTS (RETURN SUMS OVER t=1,...,ntime as above)
    # tricolon <- NULL; Gibbs_AS_posterior <- function(x, ...){NULL} # Damon added this to avoid warnings.
    # post_sums <- Gibbs_AS_posterior(tricolon, final=FALSE)
    system.time(post_sums <-
      Gibbs_AS_posteriorCPP(
        nsamp = 1000,
        nburn = 50,
        template_mean = template_mean,
        template_var = template_var,
        S = S,
        G = theta_old[[3]],
        tau_v = theta_old[[1]],
        Y = BOLD,
        alpha = theta_old[[2]],
        final = F,
        return_samp = FALSE
      ))
    S = post_sums$S_post #need to update S because it is used to initialize the Gibbs sampler

    # post_sums <-
    #   Gibbs_AS_posterior(
    #     nsamp = 100,
    #     nburn = 0,
    #     template_mean = template_mean,
    #     template_var = template_var,
    #     S = S,
    #     G = theta_old[[3]],
    #     tau_v = theta_old[[1]],
    #     Y = BOLD,
    #     alpha = theta_old[[2]],
    #     final = F
    #   )
    # S = post_sums$S_post #update S because it is used to start the Gibbs sampler

    #plot.FC_template(cov(A_init), zlim=c(-0.0002, 0.0004), break_by=0.0002, cor=FALSE, title='Initial')

    #this function returns a list of tau_sq, alpha, G
    # This is the M-step. It might be better to perform the E-step first, as the
    # M-step assumes that we have a good estimate of the first level of the posterior

    Y_sq_sum <- rowSums(BOLD^2) # This is a shortcut in computation for the CPP version
    theta_new <-
      UpdateTheta_FCtemplateICAcpp(
        template_mean,
        template_var,
        template_FC,
        theta_old[[3]],
        prior_params,
        BOLD,
        Y_sq_sum,
        post_sums,
        sigma2_alpha,
        TRUE
      )
    # theta_new <- UpdateTheta_FCtemplateICA(
    #     template_mean,
    #     template_var,
    #     template_FC,
    #     prior_params,
    #     BOLD,
    #     post_sums
    #   )

    if(verbose) print(Sys.time() - t00)

    ### Compute change in parameters

    G_old <- theta_old[[3]]
    G_new <- theta_new[[3]]
    #G_change <- max(abs((G_new - G_old)/G_old)) #not good in case there are true zeros
    #Use Frechet distance between distributions, assuming mean zero
    #https://www.sciencedirect.com/science/article/pii/0047259X8290077X
    #Tr(G_old + G_new - 2*(G_old G_new)^(1/2))
    G_orig <- sqrt(sum(diag(G_old)))
    G_diff <- sqrt(sum(diag(G_old + G_new - 2*expm::sqrtm(G_old %*% G_new))))
    G_change <- G_diff / G_orig

    tau_old <- mean(theta_old[[1]])
    tau_new <- mean(theta_new[[1]])
    tau_change <- abs((tau_new - tau_old)/tau_old)

    alpha_old <- theta_old[[2]]
    alpha_new <- theta_new[[2]]
    alpha_change <- sqrt(sum((expm::sqrtm(solve(G_new)) %*% alpha_new - expm::sqrtm(solve(G_old)) %*% alpha_old)^2))

    change <- c(G_change, tau_change, alpha_change)
    err <- max(change)
    change <- format(change, digits=3, nsmall=3)
    if(verbose) cat(paste0('Iteration ',iter, ': Difference is ',change[1],' for G, ',change[2],' for tau^2, ',change[3],' for alpha \n'))

    ### Move to next iteration
    theta_old <- theta_new
    iter <- iter + 1
    if(iter > maxiter){
      success <- 0
      warning(paste0('Failed to converge within ', maxiter,' iterations'))
      break() #exit loop
    }
  }
  #
  ### Compute final posterior means and variances of A and S

  #post_AS <- Gibbs_AS_posterior(tricolon, final=TRUE)

  post_AS <-
    Gibbs_AS_posteriorCPP(
      nsamp = 1000,
      nburn = 50,
      template_mean = template_mean,
      template_var = template_var,
      S = S,
      G = theta_new[[3]],
      tau_v = theta_new[[1]],
      Y = BOLD,
      alpha = theta_new[[2]],
      final = TRUE,
      return_samp = TRUE
    )


  #[TO DO]: Fix indexing issue where first iteration is zero (only affects burn-in)
  S_post <- array(post_AS$S_samp, dim=c(nvox, nICs, 1000-50))
  S_post_mean <- apply(S_post, c(1,2), mean)
  S_post_SE <- apply(S_post, c(1,2), sd)
  A_post <- array(post_AS$A_samp, dim=c(ntime, nICs, 1000-50))
  A_post_mean <- apply(A_post, c(1,2), mean)
  A_post_mean_cor <- cor(A_post_mean)

  #compute cor(A) for each iteration
  corA_post <- apply(A_post, 3, cor)
  corA_post_mean <- matrix(rowMeans(corA_post), nrow=nICs)
  corA_post_SE <- matrix(sqrt(matrixStats::rowVars(corA_post)), nrow=nICs)


  # delta_post <- S_post_mean - template_mean
  # delta_true <- truth_IC - template_mean
  # delta_init <- t(S_init) - template_mean
  # plot(ciftiTools::newdata_xifti(GICA, delta_true), title='true', zlim=c(-0.5,0.5))
  # plot(ciftiTools::newdata_xifti(GICA, delta_post), title='post', zlim=c(-0.5,0.5))
  # plot(ciftiTools::newdata_xifti(GICA, delta_init), title='init', zlim=c(-0.5,0.5))

  # #FC matrices
  # corA_post <- apply(post_AS$covA_final, 3, cov2cor)
  # corA_post_mean <- sapply(corA_post, mean)
  # corA_post_SE <- sapply(corA_post, sd)

  result <- list(S_post,
                 subjICmean = S_post_mean,
                 subjICse = S_post_SE,
                 A_post,
                 A_post_mean_cor,
                 corA_post,
                 subjFCmean = corA_post_mean,
                 subjFCse = corA_post_SE,
                 theta_MLE=theta_new,
                 success_flag=success,
                 error=err,
                 numiter=iter-1,
                 template_mean,
                 template_var,
                 template_FC)

  return(result)
}


#' @rdname UpdateTheta_FCtemplateICA
#' @title Update FC Template ICA parameters (tau_sq, alpha, G)
#' @param template_mean (\eqn{V \times Q} matrix) mean maps for each IC in template
#' @param template_var (\eqn{V \times Q} matrix) between-subject variance maps for each IC in template
#' @param template_FC (list) Parameters of functional connectivity template
#' @param prior_params Alpha and beta parameters of IG prior on tau^2 (error variance)
#' @param BOLD  (\eqn{V \times Q} matrix) dimension-reduced fMRI data
#' @param post_sums TBD
#' @param verbose If \code{TRUE}, display progress of algorithm. Default: \code{FALSE}.
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
  y_sq_sum <- rowSums(BOLD^2)
  for(v in 1:nvox){
    tau_sq_new[v] <- 1/(ntime + 2*alpha_tau + 2) * (post_sums$AS_sq_sum[v] - 2*post_sums$yAS_sum[v] + y_sq_sum[v] + 2*beta_tau)
  }
  tau_sq_new <- rep(mean(tau_sq_new), nvox) #TO DO: FORMALIZE THIS


  ### UPDATE ALPHA (TEMPORAL INTERCEPT)

  alpha_new <- 1/ntime*(post_sums$A_sum)


  ### UPDATE G (TEMPORAL COVARIANCE, I.E. FUNCTIONAL CONNECTIVITY)

  nu0 <- template_FC$nu
  psi0 <- template_FC$psi
  alpha_alpha_t <- tcrossprod(alpha_new)
  a_alpha_t <- tcrossprod(post_sums$A_sum, alpha_new)
  tmp <- post_sums$AtA_sum - 2*a_alpha_t + ntime*alpha_alpha_t + psi0
  G_new <- 1/(ntime + nu0 + nICs + 1)*tmp

  # RETURN NEW PARAMETER ESTIMATES
  theta_new <- list(
    tau_sq <- tau_sq_new,
    alpha <- alpha_new,
    G <- G_new
  )
  return(theta_new)
}

# R mirror of CPP function.  Would need to be updated with return_samp argument
#' #' Find expectations from the posterior of A and S using MCMC
#' #'
#' #' @param nsamp number of usable samples
#' #' @param nburn number of burn-in samples
#' #' @param template_mean V by Q matrix of template mean values
#' #' @param template_var V by Q matrix of template variance values
#' #' @param S V by Q matrix of starting values for ICs
#' #' @param G Q by Q prior covariance of A
#' #' @param tau_v vector of length V of observational variance for each data location
#' #' @param Y V by T matrix of observed BOLD values
#' #' @param alpha length Q vector prior mean of A
#' #' @param final (logical) should this output samples? Otherwise, summaries are output
#' #'
#' #' @return a list of posterior summaries or a list of MCMC samples
#' #' @importFrom Matrix Diagonal solve
#' #' @import parallel
#' #' @export
#' Gibbs_AS_posterior <- function(nsamp = 1000,
#'                                nburn = 100,
#'                                template_mean,
#'                                template_var,
#'                                S,
#'                                G,
#'                                tau_v,
#'                                Y,
#'                                alpha,
#'                                final = FALSE) {
#'   # number_samples = 1000
#'   niter <- nsamp + nburn
#'   Q = ncol(G)
#'   V = length(tau_v)
#'   tau_v <- rep(mean(tau_v), V) #TO DO: FORMALIZE THIS
#'   if(ncol(S) == V) S <- t(S) # This is to make sure the dimensions are correct
#'   ntime <- ncol(Y)
#'
#'   if(!final){
#'     A_sum = numeric(Q)
#'     AtA_sum = matrix(0,Q,Q)
#'     yAS_sum = numeric(V)
#'     AS_sq_sum = numeric(V)
#'     S_post = matrix(0, V, Q)
#'   } else {
#'     A_final = array(0, dim=c(ntime, Q))
#'     covA_final = array(0, dim=c(Q, Q, nsamp)) #save samples of cov(A)
#'     S_final = array(0, dim=c(V, Q, nsamp))
#'   }
#'
#'   G_tau_inv <- Matrix::Diagonal(x = 1/tau_v)
#'   G_inv <- solve(G)
#'   alphaGinv <- as.numeric(t(alpha)%*%G_inv)
#'   YG <- t(Y) %*% G_tau_inv
#'
#'   ### A is T by Q
#'   start_time <- proc.time()[3]
#'   for(i in 1:niter){
#'     #### update A
#'
#'     sigma_a_inv <- as.matrix(t(S) %*% G_tau_inv %*% S + G_inv)
#'     chol_sigma_a_inv <- chol(sigma_a_inv)
#'     YGS <- YG %*% S
#'     YGSalphaGinv <- as.matrix(YGS + matrix(alphaGinv, nrow=ntime, ncol=Q, byrow=TRUE))
#'     mu_a <- as.matrix(Matrix::solve(sigma_a_inv,t(YGSalphaGinv)))
#'
#'     #generate samples from N(mu_at, Sigma_a) for each t=1,...,ntime
#'     cl <- parallel::makeCluster(parallel::detectCores() - 2)
#'     A <- t(parallel::parApply(cl,mu_a,2, function(mu_at, chol_sigma_a_inv){
#'       Qa <- length(mu_at)
#'       Za <- rnorm(Qa)
#'       at <- mu_at + as.numeric(Matrix::solve(chol_sigma_a_inv,Za))
#'       return(at)
#'     }, chol_sigma_a_inv = chol_sigma_a_inv))
#'     parallel::stopCluster(cl)
#'
#'     #### update S
#'     for(v in 1:V){
#'       G_sv_inv = Matrix::Diagonal(x = 1/template_var[v,])
#'       sigma_sv_inv <- (1/tau_v[v]) * crossprod(A) + G_sv_inv
#'       chol_sigma_sv_inv <- chol(sigma_sv_inv)
#'       sigma_sv = solve(sigma_sv_inv)
#'       AtYvtempVarMean <- (1/tau_v[v]) * crossprod(A,Y[v,]) + G_sv_inv %*% template_mean[v,]
#'       mu_sv = as.numeric(Matrix::solve(sigma_sv_inv, AtYvtempVarMean))
#'       Zs <- rnorm(Q)
#'       S[v,] <- mu_sv + as.numeric(Matrix::solve(chol_sigma_sv_inv,Zs))
#'     }
#'
#'     # S = mvtnorm::rmvnorm(n = 1, mean = mu_s, sigma = sigma_s)
#'
#'     if(i > nburn) {
#'
#'       if(!final){
#'         A_sum = A_sum + colSums(A)
#'
#'         #AS_sq = sum_t (a_t's_v)^2
#'         AtS <- A %*% t(S)
#'         AS_sq_sum = AS_sq_sum + apply(AtS^2,2,sum)
#'         #AS_sq = sum_t (a_t's_v)^2
#'         yAS_sum = yAS_sum + colSums(t(Y) * AtS)
#'
#'         #AtA = sum_t a_t a_t' = A'A
#'         AtA_sum = AtA_sum + crossprod(A)
#'
#'         #S
#'         S_post = S_post + S
#'       } else {
#'         print(paste0('Saving sample ',i))
#'         A_final = A_final + A
#'         covA_final[,,i] <- cov(A)
#'         S_final[,,i] = S
#'       }
#'     }
#'
#'     if(niter >= 10 & i %% round(niter / 10) == 0) {
#'       cat("Posterior sample",i,"of",niter,". Estimated time remaining:", (proc.time()[3] - start_time) / i * (niter - i)," sec\n")
#'     }
#'   }
#'   #### return estimates of A = A_sum_init, #only actually need sum over t
#'   #AtA = AtA_sum_init,  #only actually need sum over t
#'   #yAS = yAS_sum_init,  #only actually need sum over t of Y*AS (element-wise product)
#'   #AS_sq = AS_sq_sum_init
#'
#'   if(!final){
#'     A_sum = A_sum/nsamp
#'     AtA_sum = AtA_sum/nsamp
#'     yAS_sum = yAS_sum/nsamp
#'     AS_sq_sum = AS_sq_sum/nsamp
#'     S_post = S_post/nsamp
#'   } else {
#'     A_final = A_final/nsamp
#'   }
#'   total_time <- proc.time()[3] - start_time
#'
#'   if(!final){
#'     return(list(A_sum = A_sum,
#'                 AtA_sum = AtA_sum,
#'                 yAS_sum = yAS_sum,
#'                 AS_sq_sum = AS_sq_sum,
#'                 S_post = S_post,
#'                 total_time = total_time))
#'   } else {
#'     return(list(A_final = A_final,
#'                 covA_final = covA_final,
#'                 S_final = S_final,
#'                 total_time = total_time))
#'   }
#' }
