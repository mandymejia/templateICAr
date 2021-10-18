#' @name EM_FCtemplateICA
#' @rdname EM_RFtemplateICA
#'
#' @title EM Algorithm for FC Template ICA Model
#'
#' @param template_mean (\eqn{VxQ} matrix) mean maps for each IC in template,
#'  where \eqn{Q} is the number of ICs, \eqn{V=nvox} is the number of data locations.
#' @param template_var  (\eqn{VxQ} matrix) between-subject variance maps for each IC in template
#' @param meshes \code{NULL} for spatial independence model, otherwise a list of
#'  objects of class "templateICA_mesh" containing the triangular mesh (see
#'  \code{\link{make_mesh}}) for each brain structure.
#' @param BOLD  (\eqn{VxQ} matrix) dimension-reduced fMRI data
#' @param theta0 (list) initial guess at parameter values: A (\eqn{QxQ} mixing matrix),
#'  nu0_sq (residual variance from first level) and (for spatial model only)
#'  kappa (SPDE smoothness parameter for each IC map)
#' @param C_diag (\eqn{Qx1}) diagonal elements of matrix proportional to
#'  residual variance.
# @param common_smoothness If \code{TRUE}, use the common smoothness version
#'  of the spatial template ICA model, which assumes that all IC's have the
#'  same smoothness parameter, \eqn{\kappa}
#' @param maxiter Maximum number of EM iterations. Default: 100.
#' @param epsilon Smallest proportion change between iterations. Default: 0.001.
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
EM_FCtemplateICA <- function(template_mean, template_var, BOLD, theta0, C_diag, maxiter=100, epsilon=0.001, verbose){

  #get initial values for A and S with dual regression
  #initialize the posterior expectations of A and S using those values
  #initialize theta using those posterior expectation -- this is the first update of theta using UpdateTheta_FCtemplateICA
  #UpdateTheta_FCtemplateICA will take in as an argument the current posterior moments of A and S
  #within the while loop, first run UpdateTheta_FCtemplateICA, then run Gibbs_AS_posterior to get required moments

  #in Gibbs_AS_posterior, have an argument to determine whether you want the final posterior mean/variance or the elements necessary for the M-step
  #Ani pointed out that we don't want to save all of the samples for computational reasons.  We can compute sums as we go.  Can compute sums in chunks of samples of 100, say, then decide how many chunks to drop before calculating the final means.


  # if(!all.equal(dim(template_var), dim(template_mean))) stop('The dimensions of template_mean and template_var must match.')
  #
  # ntime <- ncol(BOLD) #length of timeseries
  # nvox <- nrow(BOLD) #number of brain locations
  # if(ntime > nvox) warning('More time points than brain locations. Are you sure?')
  # if(nrow(template_mean) != nvox) stop('Templates and BOLD must have the same number of brain locations (columns).')
  #
  # Q <- ncol(template_mean) #number of ICs
  # if(Q > nvox) stop('Cannot estimate more ICs than brain locations.')
  # if(Q > ntime) stop('Cannot estimate more ICs than time points.')
  #
  # iter <- 1
  # theta <- theta0
  # success <- 1
  # template_var[template_var < 1e-6] <- 1e-6 #to prevent problems when inverting covariance
  #
  # err <- 1000 #large initial value for difference between iterations
  # while(err > epsilon){
  #
  #   if(verbose) cat(paste0(' ~~~~~~~~~~~~~~~~~~~~~ ITERATION ', iter, ' ~~~~~~~~~~~~~~~~~~~~~ \n'))
  #
  #   t00 <- Sys.time()
  #   theta_new = UpdateTheta_FCtemplateICA(template_mean, template_var, BOLD, theta, C_diag, verbose=verbose)
  #   if(verbose) print(Sys.time() - t00)
  #
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
  # }
  #
  # ### Compute final posterior mean of subject ICs
  #
  # #A = theta$A
  # At_nu0Cinv <- t(theta$A) %*% diag(1/(C_diag*theta$nu0_sq))
  # At_nu0Cinv_A <- At_nu0Cinv %*% theta$A
  # miu_s <- matrix(NA, nrow=nvox, ncol=Q)
  # var_s <- matrix(NA, nrow=nvox, ncol=Q)
  # for(v in 1:nvox){
  #   y_v <- BOLD[v,]
  #   s0_v <- template_mean[v,]
  #   E_v_inv <- diag(1/template_var[v,])
  #   Sigma_s_v <- solve(E_v_inv + At_nu0Cinv_A)
  #   miu_s[v,] <- Sigma_s_v	%*% (At_nu0Cinv %*% y_v + E_v_inv %*% s0_v) #Qx1
  #   var_s[v,] <- diag(Sigma_s_v)
  # }
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
UpdateTheta_FCtemplateICA <- function(template_mean, template_var, BOLD, theta, C_diag, verbose){

  #Update theta = (tau_sq, alpha, G)

  # Q <- ncol(BOLD)
  # nvox <- nrow(BOLD)
  #
  # #initialize new objects
  # theta_new <- list(A = matrix(NA, Q, Q), nu0_sq = NA)
  # A_part1 <- A_part2 <- matrix(0, Q, Q) #two parts of product for A-hat (construct each looping over voxels)
  #
  #
  # A <- theta$A
  # nu0_sq <- theta$nu0_sq
  # nu0C_inv <- diag(1/(C_diag*nu0_sq)) #Sigma0_inv in matlab code
  # At_nu0Cinv <- t(A) %*% nu0C_inv
  # At_nu0Cinv_A <- At_nu0Cinv %*% A
  #
  # if(verbose) cat('Updating A \n')
  #
  # #store posterior moments for M-step of nu0_sq
  # miu_s <- matrix(NA, nrow=nvox, ncol=Q)
  # miu_ssT <- array(NA, dim=c(nvox, Q, Q))
  #
  # for(v in 1:nvox){
  #
  #   y_v <- BOLD[v,]
  #   s0_v <- template_mean[v,]
  #
  #   ##########################################
  #   ### E-STEP FOR A AND nu0^2: POSTERIOR MOMENTS OF s_i(v)
  #   ##########################################
  #
  #   E_v_inv <- diag(1/template_var[v,])
  #   Sigma_s_v <- solve(E_v_inv + At_nu0Cinv_A)
  #   miu_s_v <- Sigma_s_v	%*% (At_nu0Cinv %*% y_v + E_v_inv %*% s0_v) #Qx1
  #   miu_ssT_v <- (miu_s_v %*% t(miu_s_v)) + Sigma_s_v #QxQ
  #   miu_s[v,] <- miu_s_v #save for M-step of nu0_sq
  #   miu_ssT[v,,] <- miu_ssT_v #save for M-step of nu0_sq
  #
  #   ##########################################
  #   ### M-STEP FOR A: CONSTRUCT PARAMETER ESTIMATES
  #   ##########################################
  #
  #   A_part1 <- A_part1 + y_v %*% t(miu_s_v) #QxQ
  #   A_part2 <- A_part2 + miu_ssT_v #QxQ
  #
  # }
  #
  # #A_hat <- orthonorm(A_part1 %*% solve(A_part2))
  # A_hat <- (A_part1 %*% solve(A_part2))
  #
  # ##########################################
  # ### M-STEP FOR nu0^2: CONSTRUCT PARAMETER ESTIMATES
  # ##########################################
  #
  # # cat('Updating Error Variance nu0_sq \n')
  # #
  # # #use A-hat or A?
  # #
  # # Cinv <- diag(1/C_diag)
  # # Cinv_A <- Cinv %*% A_hat
  # # At_Cinv_A <- t(A_hat) %*% Cinv %*% A_hat
  # # nu0sq_part1 <- nu0sq_part2 <- nu0sq_part3 <- 0
  # #
  # # for(v in 1:nvox){
  # #
  # #   y_v <- BOLD[,v]
  # #   nu0sq_part1 <- nu0sq_part1 + t(y_v) %*% Cinv %*% y_v
  # #   nu0sq_part2 <- nu0sq_part2 + t(y_v) %*% Cinv_A %*% miu_s[,v]
  # #   nu0sq_part3 <- nu0sq_part3 + sum(diag(At_Cinv_A %*% miu_ssT[,,v]))
  # # }
  # #
  # # nu0sq_hat <- 1/(Q*nvox)*(nu0sq_part1 - 2*nu0sq_part2 + nu0sq_part3)
  #
  # nu0sq_hat <- theta$nu0_sq
  #
  #
  # # RETURN NEW PARAMETER ESTIMATES
  #
  # theta_new$A <- A_hat
  # theta_new$nu0_sq <- nu0sq_hat[1]
  # return(theta_new)
}
