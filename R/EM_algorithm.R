#' @name EM_Algorithm
#' @rdname EM_Algorithm
#'
#' @title EM Algorithms for Template ICA Models
#'
#' @param template_mean (QxN matrix) mean maps for each IC in template, where Q is the number of ICs, N is the number of data or mesh locations.
#' @param template_var  (QxN matrix) between-subject variance maps for each IC in template
#' @param mesh NULL for spatial independence model, otherwise an object of class "templateICA_mesh" containing the triangular mesh (see `help(make_mesh)`)
#' @param BOLD  (QxN matrix) dimension-reduced fMRI data
#' @param theta0 (list) initial guess at parameter values: A (QxQ mixing matrix), nu0_sq (residual variance from first level) and (for spatial model only) kappa (SPDE smoothness parameter for each IC map)
#' @param C_diag (Qx1) diagonal elements of matrix proportional to residual variance.
#' @param common_smoothness If TRUE, use the common smoothness version of the spatial template ICA model, which assumes that all IC's have the same smoothness parameter, \eqn{\kappa}
#' @param maxiter maximum number of EM iterations
#' @param epsilon smallest proportion change between iterations (e.g. .001)
#' @param verbose If TRUE, display progress of algorithm
#' @param dim_reduce_flag If FALSE, data is in the original resolution (no dimension reduction).
#'
#' @return  A list with 4 elements: theta (list of final parameter estimates), subICmean (estimates of subject-level ICs), subICvar (variance of subject-level ICs), and success (flag indicating convergence (\code{TRUE}) or not (\code{FALSE}))
#'
#' @details \code{EM_templateICA.spatial} implements the expectation-maximization (EM) algorithm described in Mejia et al. (2019+) for estimating the subject-level ICs and unknown parameters in the template ICA model with spatial priors on subject effects.
#'
#' In both models, if original fMRI timeseries has covariance \eqn{\sigma^2 I_T}, the prewhitened timeseries achieved by premultiplying by (QxT) matrix \eqn{H} from PCA has diagonal covariance \eqn{\sigma^2HH'}, so C_diag is \eqn{diag(HH')}.
#'
#'
NULL

#' @rdname EM_Algorithm
#' @export
#' @importFrom INLA inla.spde2.matern inla.qsolve
#' @importFrom Matrix Diagonal
#' @import SQUAREM
#'
EM_templateICA.spatial = function(template_mean, template_var, mesh, BOLD, theta0, C_diag, common_smoothness=TRUE, maxiter=100, epsilon=0.001, verbose=FALSE, dim_reduce_flag){

  if(!all.equal(dim(template_var), dim(template_mean))) stop('The dimensions of template_mean and template_var must match.')

  ntime <- nrow(BOLD) #length of timeseries
  nvox <- ncol(BOLD) #number of data locations
  if(ntime > nvox) warning('More time points than data locations. Are you sure?')
  if(ncol(template_mean) != nvox) stop('Templates and BOLD must have the same number of data locations (columns).')

  Q <- nrow(template_mean) #number of ICs
  if(Q > nvox) stop('Cannot estimate more ICs than data locations.')
  if(Q > ntime) stop('Cannot estimate more ICs than time points.')

  if(class(mesh) != 'templateICA_mesh') stop('mesh argument should be of class templateICA_mesh. See help(make_mesh).')

	iter = 1
	theta = theta0
	success = 1

	template_var[template_var==0] <- 1e-5

	#pre-compute s0, D and D^{-1}*s0
	V <- ncol(template_mean)
	s0_vec = as.vector(t(template_mean))
	D_vec <- as.vector(sqrt(t(template_var))) #template_var is QxV
	D = Diagonal(V*Q, D_vec)
	Dinv_s0 <- inla.qsolve(Q = D, B=matrix(s0_vec, ncol=1), method='solve')

  ### REFINE STARTING VALUE FOR KAPPA

	if(verbose) cat('Refining starting value for kappa \n')

	# Determine direction of change:
	# Positive change --> search for kappa_max, set kappa_min to kappa1.
	# Negative change --> search for kappa_min, set kappa_max to kappa1.
	kappa_min <- kappa_max <- theta0$kappa[1]
	theta1 <- UpdateTheta.spatial(template_mean, template_var, mesh, BOLD, theta0, C_diag, s0_vec, D, Dinv_s0, common_smoothness=TRUE, verbose=FALSE, dim_reduce_flag=dim_reduce_flag, update='kappa')
	kappa_diff0 <- theta1$kappa[1] - theta0$kappa[1]
	theta <- theta0

	kappa_diff <- kappa_diff0
	if(kappa_diff0 < 0){

	  if(verbose) cat('...Kappa decreasing, finding lower bound for kappa search \n ')

	  kappa_min <- kappa_min/2
	  while(kappa_diff < 0){
	    if(verbose) cat(paste0('... testing kappa = ',round(kappa_min,3),'\n '))
	    theta$kappa <- rep(kappa_min, Q)
	    theta1 <- UpdateTheta.spatial(template_mean, template_var, mesh, BOLD, theta, C_diag, s0_vec, D, Dinv_s0, common_smoothness=TRUE, verbose=FALSE, dim_reduce_flag=dim_reduce_flag, update='kappa')
	    kappa_diff <- theta1$kappa[1] - theta$kappa[1]
	    if(kappa_diff > 0) {
	      #set minimum and stop here
	      kappa_min <- theta1$kappa[1]
	      break
	    } else {
	      #set new value for kappa
	      kappa_min <- kappa_min/2
	    }
	  }
	} else if(kappa_diff0 > 0){

	  if(verbose) cat('...Kappa increasing, finding upper bound for kappa search \n ')

	  kappa_max <- kappa_max*2
	  while(kappa_diff > 0){
	    if(verbose) cat(paste0('... testing kappa = ',round(kappa_max, 3),'\n '))
	    theta$kappa <- rep(kappa_max, Q)
	    theta1 <- UpdateTheta.spatial(template_mean, template_var, mesh, BOLD, theta, C_diag, s0_vec, D, Dinv_s0, common_smoothness=TRUE, verbose=FALSE, dim_reduce_flag=dim_reduce_flag, update='kappa')
	    kappa_diff <- theta1$kappa[1] - theta$kappa[1]
	    if(kappa_diff < 0) {
	      #set maximum and stop here
	      kappa_max <- theta1$kappa[1]
	      break
	    } else {
	      #set new value for kappa
	      kappa_max <- kappa_max*2
	    }
	  }
	}

	#use binary search until convergence
	if(verbose) cat('...Starting binary search for starting value of kappa \n ')
	kappa_test <- (kappa_min + kappa_max)/2
	kappa_change <- 1
	while(kappa_change > epsilon){
	  if(verbose) cat(paste0('... testing kappa = ',round(kappa_test, 3),'\n '))
	  theta$kappa <- rep(kappa_test, Q)
	  theta1 <- UpdateTheta.spatial(template_mean, template_var, mesh, BOLD, theta, C_diag, s0_vec, D, Dinv_s0, common_smoothness=TRUE, verbose=FALSE, dim_reduce_flag=dim_reduce_flag, update='kappa')
	  kappa_diff <- theta1$kappa[1] - theta$kappa[1] #which direction is the estimate of kappa moving in?
	  if(kappa_diff > 0) {
	    kappa_min <- theta1$kappa[1]  #reset minimum to current value
	    kappa_test <- (theta1$kappa[1] + kappa_max)/2 #go halfway to max
	  } else {
	    kappa_max <- theta1$kappa[1]  #reset maximum to current value
	    kappa_test <- (theta1$kappa[1] + kappa_min)/2 #go halfway to min
	  }

	  #how much different is the next value of kappa to be tested versus the current one?
	  kappa_change <- abs((kappa_test - theta1$kappa[1])/theta1$kappa[1])
	}


	### RUN SQUAREM ALGORITHM UNTIL CONVERGENCE

	theta0 <- theta1 #last tested value of kappa0
	theta0$LL <- c(0,0)
	theta0_vec <- unlist(theta0[1:3]) #everything but LL
	names(theta0_vec)[1] <- 0 #store LL value in names of theta0_vec (required for squarem)

	t00000 <- Sys.time()
	result_squarem <- squarem(par=theta0_vec, fixptfn = UpdateThetaSQUAREM, objfn=LL_SQUAREM, control=list(trace=verbose, intermed=TRUE, tol=epsilon, maxiter=maxiter), template_mean, template_var, mesh, BOLD, C_diag, s0_vec, D, Dinv_s0, common_smoothness, verbose=FALSE, dim_reduce_flag)
	if(verbose) print(Sys.time() - t00000)

	path_A <- result_squarem$p.inter[,1:(Q^2)]
	path_kappa <- result_squarem$p.inter[,(Q^2+1)+(1:Q)]
	path_LL <- result_squarem$p.inter[,ncol(result_squarem$p.inter)]
	theta_path <- list(A=path_A, kappa=path_kappa, LL=path_LL)

	theta_MLE <- theta0
	theta_MLE$A <- matrix(result_squarem$par[1:(Q^2)], Q, Q)
	theta_MLE$kappa <- result_squarem$par[(Q^2+1)+(1:Q)]
	theta_MLE$LL <- as.numeric(names(result_squarem$par)[1])

	#success <- (result_squarem$convergence==0) #0 indicates convergence, 1 indicates failure to converge within maxiter
	numiter <- result_squarem$fpevals #number of parameter update steps (approximately 3x the number of SQUAREM iterations)

	### Compute final posterior mean of subject ICs
	if(verbose) cat('Computing final posterior mean of subject ICs \n')
	mu_Omega_s = UpdateTheta.spatial(template_mean, template_var, mesh, BOLD, theta_MLE, C_diag, s0_vec, D, Dinv_s0, common_smoothness=common_smoothness, verbose=verbose, return_MAP=TRUE)

	subjICcov=(D %*% inla.qinv(mu_Omega_s$Omega_s) %*% D)

	result <- list(subjICmean=mu_Omega_s$mu_s, subjICcov=subjICcov, Omega = mu_Omega_s$Omega_s, theta_MLE=theta_MLE, theta_path=theta_path, numiter=numiter, squarem = result_squarem, template = list(mean = t(template_mean), var=t(template_var)))
	return(result)
}

#' @rdname EM_Algorithm
#' @export
EM_templateICA.independent = function(template_mean, template_var, BOLD, theta0, C_diag, maxiter=100, epsilon=0.001){

  if(!all.equal(dim(template_var), dim(template_mean))) stop('The dimensions of template_mean and template_var must match.')

  ntime <- nrow(BOLD) #length of timeseries
  nvox <- ncol(BOLD) #number of brain locations
  if(ntime > nvox) warning('More time points than brain locations. Are you sure?')
  if(ncol(template_mean) != nvox) stop('Templates and BOLD must have the same number of brain locations (columns).')

  Q <- nrow(template_mean) #number of ICs
  if(Q > nvox) stop('Cannot estimate more ICs than brain locations.')
  if(Q > ntime) stop('Cannot estimate more ICs than time points.')

  iter = 1
  theta = theta0
  success = 1
  template_var[template_var < .00001] = .00001 #to prevent problems when inverting covariance

  err = 1000 #large initial value for difference between iterations
  while(err > epsilon){

    print(paste0(' ~~~~~~~~~~~~~~~~~~~~~ ITERATION ', iter, ' ~~~~~~~~~~~~~~~~~~~~~ '))

    t00 <- Sys.time()
    theta_new = UpdateTheta.independent(template_mean, template_var, BOLD, theta, C_diag)
    print(Sys.time() - t00)

    ### Compute change in parameters

    A_old = theta$A
    A_new = theta_new$A
    #2-norm = largest eigenvalue = sqrt of largest eigenvalue of AA'
    A_change = norm(as.vector(A_new - A_old), type="2")/norm(as.vector(A_old), type="2")

    nu0_sq_old = theta$nu0_sq
    nu0_sq_new = theta_new$nu0_sq
    nu0_sq_change = abs(nu0_sq_new - nu0_sq_old)/nu0_sq_old

    change = c(A_change, nu0_sq_change)
    err = max(change)
    change = format(change, digits=3, nsmall=3)
    print(paste0('Iteration ',iter, ': Difference is ',change[1],' for A, ',change[2],' for nu0_sq'))

    ### Move to next iteration
    theta <- theta_new
    iter = iter + 1
    if(iter > maxiter){
      success = 0;
      warning(paste0('Failed to converge within ', maxiter,' iterations'))
      break() #exit loop
    }
  }

  ### Compute final posterior mean of subject ICs

  A = theta$A
  At_nu0Cinv = t(theta$A) %*% diag(1/(C_diag*theta$nu0_sq))
  At_nu0Cinv_A = At_nu0Cinv %*% theta$A
  miu_s = matrix(NA, nrow=Q, ncol=nvox)
  var_s = matrix(NA, nrow=Q, ncol=nvox)
  for(v in 1:nvox){
    y_v <- BOLD[,v]
    s0_v <- template_mean[,v]
    E_v_inv <- diag(1/template_var[,v])
    Sigma_s_v <- solve(E_v_inv + At_nu0Cinv_A)
    miu_s[,v] <- Sigma_s_v	%*% (At_nu0Cinv %*% y_v + E_v_inv %*% s0_v) #Qx1
    var_s[,v] <- diag(Sigma_s_v)
  }

  result <- list(subjICmean=t(miu_s), subjICvar=t(var_s), theta_MLE=theta, success_flag=success, error=err, numiter=iter-1, template = list(mean = t(template_mean), var = t(template_var)))
  #names(result) <- c('subjICmean', 'subjICvar', 'theta_MLE', 'success_flag')
  return(result)
}


# my_squarem <- function (par, fixptfn, objfn, ..., control = list())
# {
#   control.default <- list(K = 1, method = 3, square = TRUE,
#                           step.min0 = 1, step.max0 = 1, mstep = 4, kr = 1, objfn.inc = 1,
#                           tol = 1e-07, maxiter = 1500, trace = FALSE, intermed = FALSE)
#   namc <- names(control)
#   if (!all(namc %in% names(control.default)))
#     stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
#   ctrl <- modifyList(control.default, control)
#   if (ctrl$K > 1 & !(ctrl$method %in% c("rre", "mpe")))
#     ctrl$method <- "rre"
#   if (ctrl$K == 1 & !(ctrl$method %in% c(1, 2, 3)))
#     ctrl$method <- 3
#   if (!missing(objfn)) {
#     if (ctrl$K == 1)
#       sqobj <- my_squarem1(par, fixptfn, objfn, ..., control = ctrl)
#     else if (ctrl$K > 1 | ctrl$method %in% c("rre", "mpe"))
#       sqobj <- cyclem1(par, fixptfn, objfn, ..., control = ctrl)
#   }
#   else {
#     if (ctrl$K == 1)
#       sqobj <- squarem2(par, fixptfn, ..., control = ctrl)
#     else if (ctrl$K > 1 | ctrl$method %in% c("rre", "mpe"))
#       sqobj <- cyclem2(par, fixptfn, ..., control = ctrl)
#   }
#   return(sqobj)
# }
#
#
# my_squarem1 <- function (par, fixptfn, objfn, ..., control = list())
# {
#   control.default <- list(K = 1, square = TRUE, method = 3,
#                           step.min0 = 1, step.max0 = 1, mstep = 4, kr = 1, objfn.inc = 1,
#                           tol = 1e-07, maxiter = 1500, trace = FALSE, intermed = FALSE)
#   namc <- names(control)
#   if (!all(namc %in% names(control.default)))
#     stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
#   ctrl <- modifyList(control.default, control)
#   method <- ctrl$method
#   maxiter <- ctrl$maxiter
#   tol <- ctrl$tol
#   step.min <- ctrl$step.min0
#   step.max0 <- ctrl$step.max0
#   step.max <- ctrl$step.max0
#   mstep <- ctrl$mstep
#   objfn.inc <- ctrl$objfn.inc
#   trace <- ctrl$trace
#   intermed <- ctrl$intermed
#   if (trace)
#     cat("Squarem-1 \n")
#   if (missing(objfn))
#     stop("\n squarem2 should be used if objective function is not available \n\n")
#   iter <- 1
#   objval <- rep(NA, 1)
#   p <- par
#   path_p[[1]] <- p
#   lold <- objfn(p, ...)
#   leval <- 1
#   if (trace)
#     cat("Objective fn: ", lold, "\n")
#   feval <- 0
#   conv <- TRUE
#   p.inter <- c(p, lold)
#   path_p <- vector('list', length=maxiter)
#   while (feval < maxiter) {
#     extrap <- TRUE
#     p1 <- try(fixptfn(p, ...), silent = TRUE)
#     feval <- feval + 1
#     path_p[[feval]] <- p1
#     if (inherits(p1, "try-error") | any(is.nan(unlist(p1))))
#       stop("Error in function evaluation")
#     q1 <- p1 - p
#     sr2 <- crossprod(q1)
#     if (sqrt(sr2) < tol)
#       break
#     p2 <- try(fixptfn(p1, ...), silent = TRUE)
#     feval <- feval + 1
#     path_p[[feval]] <- p2
#     if (inherits(p2, "try-error") | any(is.nan(unlist(p2))))
#       stop("Error in function evaluation")
#     q2 <- p2 - p1
#     sq2 <- sqrt(crossprod(q2))
#     if (sq2 < tol)
#       break
#     sv2 <- crossprod(q2 - q1)
#     srv <- crossprod(q1, q2 - q1)
#     alpha <- switch(method, -srv/sv2, -sr2/srv, sqrt(sr2/sv2))
#     alpha <- max(step.min, min(step.max, alpha))
#     p.new <- p + 2 * alpha * q1 + alpha^2 * (q2 - q1)
#     if (abs(alpha - 1) > 0.01) {
#       p.new <- try(fixptfn(p.new, ...), silent = TRUE)
#       feval <- feval + 1
#     }
#     if (inherits(p.new, "try-error") | any(is.nan(p.new))) {
#       p.new <- p2
#       lnew <- try(objfn(p2, ...), silent = TRUE)
#       leval <- leval + 1
#       if (alpha == step.max)
#         step.max <- max(step.max0, step.max/mstep)
#       alpha <- 1
#       extrap <- FALSE
#     }
#     else {
#       if (is.finite(objfn.inc)) {
#         lnew <- try(objfn(p.new, ...), silent = TRUE)
#         leval <- leval + 1
#       }
#       else lnew <- lold
#       if (inherits(lnew, "try-error") | is.nan(lnew) |
#           (lnew > lold + objfn.inc)) {
#         p.new <- p2
#         lnew <- try(objfn(p2, ...), silent = TRUE)
#         leval <- leval + 1
#         if (alpha == step.max)
#           step.max <- max(step.max0, step.max/mstep)
#         alpha <- 1
#         extrap <- FALSE
#       }
#     }
#     if (alpha == step.max)
#       step.max <- mstep * step.max
#     if (step.min < 0 & alpha == step.min)
#       step.min <- mstep * step.min
#     p <- p.new
#     path_p[[feval]] <- p.new
#     if (!is.nan(lnew))
#       lold <- lnew
#     if (trace)
#       cat("Objective fn: ", lnew, "  Extrapolation: ",
#           extrap, "  Steplength: ", alpha, "\n")
#     if (intermed)
#       p.inter <- rbind(p.inter, c(p, lnew))
#     iter <- iter + 1
#     path_p[[iter]] <-
#   }
#   path_p <- path_p[1:iter]
#   if (feval >= maxiter)
#     conv <- FALSE
#   if (is.infinite(objfn.inc)) {
#     lold <- objfn(p, ...)
#     leval <- leval + 1
#   }
#   rownames(p.inter) <- NULL
#   if (!intermed)
#     return(list(par = p, path = path_p, value.objfn = lold, iter = iter,
#                 fpevals = feval, objfevals = leval, convergence = conv))
#   else return(list(par = p, path = path_p, value.objfn = lold, iter = iter,
#                    fpevals = feval, objfevals = leval, convergence = conv,
#                    p.inter = p.inter))
# }
#
#
#
#
