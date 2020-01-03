#' @name EM_Algorithm
#' @rdname EM_Algorithm
#'
#' @title EM Algorithms for Template ICA Models
#'
#' @param template_mean (QxN matrix) mean maps for each IC in template, where Q is the number of ICs, N is the number of data or mesh locations.
#' @param template_var  (QxN matrix) between-subject variance maps for each IC in template
#' @param BOLD  (QxN matrix) dimension-reduced fMRI data
#' @param mesh NULL for spatial independence model, otherwise an object of class "templateICA_mesh" containing the triangular mesh (see `help(make_mesh)`)
#' @param theta0 (list) initial guess at parameter values: A (QxQ mixing matrix), nu0_sq (residual variance from first level) and (for spatial model only) kappa (SPDE smoothness parameter for each IC map)
#' @param C_diag (Qx1) diagonal elements of matrix proportional to residual variance.
#' @param common_smoothness If TRUE, use the common smoothness version of the spatial template ICA model, which assumes that all IC's have the same smoothness parameter, \eqn{\kappa}
#' @param maxiter maximum number of EM iterations
#' @param epsilon smallest proportion change between iterations (e.g. .001)
#' @param return_kappa_fun If TRUE, return the log likelihood as a function of kappa in a neighborhood of the MLE (common smoothness model only)
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
#' @importFrom INLA inla.spde2.matern
#'
EM_templateICA.spatial = function(template_mean, template_var, mesh, BOLD, theta0, C_diag, common_smoothness=TRUE, maxiter=100, epsilon=0.01, return_kappa_fun=FALSE, verbose=FALSE){

  if(!all.equal(dim(template_var), dim(template_mean))) stop('The dimensions of template_mean and template_var must match.')

  ntime <- nrow(BOLD) #length of timeseries
  nmesh <- ncol(BOLD) #number of mesh locations
  if(ntime > nmesh) warning('More time points than mesh locations. Are you sure?')
  if(ncol(template_mean) != nmesh) stop('Templates and BOLD must have the same number of mesh locations (columns).')

  Q <- nrow(template_mean) #number of ICs
  if(Q > nmesh) stop('Cannot estimate more ICs than mesh locations.')
  if(Q > ntime) stop('Cannot estimate more ICs than time points.')

  if(class(mesh) != 'templateICA_mesh') stop('mesh argument should be of class templateICA_mesh. See help(make_mesh).')
  spde <- inla.spde2.matern(mesh$mesh, alpha=2)
  if(spde$n.spde != nmesh) stop('The mesh does not have the right number of locations.')

	iter = 1
	theta = theta0
	success = 1
	template_var[template_var < .00001] = .00001 #to prevent problems when inverting covariance

	err = 1000 #large initial value for difference between iterations
	theta_path <- vector('list', length=maxiter)
	while(err > epsilon){

	  print(paste0(' ~~~~~~~~~~~~~~~~~~~~~ ITERATION ', iter, ' ~~~~~~~~~~~~~~~~~~~~~ '))

		t00 <- Sys.time()
		theta_new = UpdateTheta.spatial(template_mean, template_var, spde, BOLD, theta, C_diag, common_smoothness=common_smoothness, verbose=verbose, return_kappa_fun=return_kappa_fun)
		print(Sys.time() - t00)

		### Compute change in parameters

		A_old = theta$A
		A_new = theta_new$A
		#2-norm = largest eigenvalue = sqrt of largest eigenvalue of AA'
		A_change = norm(as.vector(A_new - A_old), type="2")/norm(as.vector(A_old), type="2")

		nu0_sq_old = theta$nu0_sq
		nu0_sq_new = theta_new$nu0_sq
		nu0_sq_change = abs(nu0_sq_new - nu0_sq_old)/nu0_sq_old

		kappa_old = theta$kappa
		kappa_new = theta_new$kappa
		kappa_change = norm(kappa_new - kappa_old, type="2")/norm(kappa_old, type="2")

		change = c(A_change, nu0_sq_change, kappa_change)
		err = max(change)
		change = format(change, digits=3, nsmall=3)
		print(paste0('Iteration ',iter, ': Difference is ',change[1],' for A, ',change[2],' for nu0_sq and ',change[3],' for kappa'))
    print(paste0('Current estimate of kappa: ', round(kappa_new[1], 3)))

		### Move to next iteration
		theta <- theta_new
		theta_path[[iter]] <- theta_new
		iter = iter + 1
		if(iter > maxiter){
			success = 0;
			warning(paste0('Failed to converge within ', maxiter,' iterations'))
			break() #exit loop
		}
	}

	### Compute final posterior mean of subject ICs
	print('Computing final posterior mean of subject ICs')
	mu_Omega_s = UpdateTheta.spatial(template_mean, template_var, spde, BOLD, theta, C_diag, common_smoothness=common_smoothness, verbose=verbose, return_kappa_fun=return_kappa_fun, return_MAP=TRUE)

	result <- list(subjICmean=mu_Omega_s$mu_s, subjICprec=mu_Omega_s$Omega_s, theta_MLE=theta, theta_path=theta_path, success_flag=success, error=err, numiter=iter-1)
	return(result)
}

#' @rdname EM_Algorithm
#' @export
EM_templateICA.independent = function(template_mean, template_var, BOLD, theta0, C_diag, maxiter=100, epsilon=0.01){

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

  result <- list(subjICmean=miu_s, subjICvar=var_s, theta_MLE=theta, success_flag=success, error=err, numiter=iter-1)
  #names(result) <- c('subjICmean', 'subjICvar', 'theta_MLE', 'success_flag')
  return(result)
}
















