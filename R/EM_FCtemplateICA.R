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
#' @param AS_init (list) initial guess at latent variables: A (\eqn{TxQ} mixing matrix),
#'  and S (\eqn{QxV} matrix of spatial ICs)
#' @param maxiter Maximum number of EM iterations. Default: 100.
#' @param epsilon Smallest proportion change in parameter estimates between iterations. Default: 0.001.
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
                             prior_params=c(0.001,0.001),
                             BOLD,
                             AS_init,
                             maxiter=100,
                             epsilon=0.001,
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
  template_var[template_var < 1e-6] <- 1e-6 #to prevent problems when inverting covariance

  #compute initial estimates of posterior moments
  A_init <- AS_0$A #TxQ
  S_init <- AS_0$S #QxV
  AS_init <- A_init %*% S_init #TxV -- for the terms E[a_t's_v]
  AS_sq_init <- (AS_init^2) #TxV -- for the terms E[(a_t's_v)^2]

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

    t00 <- Sys.time()
    #this function returns a list of tau_sq, alpha, G
    theta_new <- UpdateTheta_FCtemplateICA(template_mean,
                                          template_var,
                                          template_FC,
                                          prior_params,
                                          BOLD,
                                          post_sums,
                                          verbose=verbose)
    if(verbose) print(Sys.time() - t00)

    ### TO DO: RUN GIBBS SAMPLER TO SAMPLE FROM (A,S) AND UPDATE POSTERIOR_MOMENTS (RETURN SUMS OVER t=1,...,ntime as above)

    post_sums <- Gibbs_AS_posterior(..., final=FALSE)


  #   ### Compute change in parameters
  #
  #   A_old <- theta$A
  #   A_new <- theta_new$A
  #   #2-norm <- largest eigenvalue <- sqrt of largest eigenvalue of AA'
  #   A_change <- norm(as.vector(A_new - A_old), type="2")/norm(as.vector(A_old), type="2")
  #
  #   nu0_sq_old <- theta$nu0_sq
  #   nu0_sq_new <- theta_new$nu0_sq
  #   nu0_sq_change <- abs(nu0_sq_new - nu0_sq_old)/nu0_sq_old
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

  post_AS <- Gibbs_AS_posterior(..., final=TRUE)

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
