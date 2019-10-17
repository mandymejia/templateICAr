#' EM Algorithm for Template ICA with Spatial Priors
#'
#' @description Implements the expectation-maximization (EM) algorithm described in Mejia et al. (2019+) for estimating the subject-level ICs and unknown parameters in the template ICA with spatial priors model.
#'
#' @param template_mean (QxN matrix) mean maps for each IC in template. Q is the number of ICs, N is the number of mesh locations.
#' @param template_var  (QxN matrix) between-subject variance maps for each IC in template
#' @param mesh Object of class "templateICA_mesh" containing the triangular mesh (see `help(make_mesh)`)
#' @param BOLD (QxN matrix) dimension-reduced fMRI data at each mesh location.
#' @param theta0 (list) initial guess at parameter values: A (QxQ mixing matrix), nu0_sq (residual variance from first level) and kappa (SPDE smoothness parameter for each IC map)
#' @param C_diag (Qx1) diagonal elements of matrix proportional to residual variance.
#' @param max_iter maximum number of EM iterations
#' @param epsilon smallest proportion change between iterations (e.g. .001)
#'
#' @return  A list with 4 elements: theta (list of final parameter estimates), subICmean (estimates of subject-level ICs), subICvar (variance of subject-level ICs), and success (flag indicating convergence (\code{TRUE}) or not (\code{FALSE}))
#' @export
#' @importFrom INLA inla.spde2.matern
#'
#' @details If original fMRI timeseries has covariance \eqn{\sigma^2 I_T}, the prewhitened timeseries achieved by premultiplying by (QxT) matrix \eqn{H} from PCA has diagonal covariance \eqn{\sigma^2HH'}, so C_diag is \eqn{diag(HH')}.
#'
EM_templateICA.spatial = function(template_mean, template_var, mesh, BOLD, theta0, C_diag, max_iter=100, epsilon=0.01){

  if(!all.equal(dim(template_var), dim(template_mean))) error('The dimensions of template_mean and template_var must match.')

  ntime <- nrow(BOLD) #length of timeseries
  nmesh <- ncol(BOLD) #number of mesh locations
  if(ntime > nmesh) warning('More time points than mesh locations. Are you sure?')
  if(ncol(template_mean) != nmesh) error('Templates and BOLD must have the same number of mesh locations (columns).')

  Q <- nrow(template_mean) #number of ICs
  if(Q > nmesh) error('Cannot estimate more ICs than mesh locations.')
  if(Q > ntime) error('Cannot estimate more ICs than time points.')

  spde <- inla.spde2.matern(mesh$mesh, alpha=2)
  if(spde$n.spde != nmesh) error('The mesh does not have the right number of locations.')


	iter = 1
	theta = theta0
	success = 1
	template_var[template_var < .00001] = .00001 #to prevent problems when inverting covariance

	err = 1000 #large initial value for difference between iterations
	while(err > epsilon){

	  print(paste0('\n ~~~~~~~~~~~~~~~~~~~~~ ITERATION ', iter, ' ~~~~~~~~~~~~~~~~~~~~~ \n'))

		t00 <- Sys.time()
		theta_new = UpdateTheta(template_mean, template_var, spde, BOLD, theta, C_diag)
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

		### Move to next iteration
		theta <- theta_new
		iter = iter + 1
		if(iter > max_iter){
			success = 0;
			warning(paste0('Failed to converge within ', max_iter,' iterations'))
			break() #exit loop
		}
	}

	### Compute final posterior mean of subject ICs

	A = theta$A
	At_nu0Cinv = t(theta$A) %*% diag(1/(C_diag*theta$nu0_sq))
	At_nu0Cinv_A = At_nu0Cinv %*% theta$A
	miu_s = matrix(NA, nrow=Q, ncol=V)
	var_s = matrix(NA, nrow=Q, ncol=V)
	for(v in 1:V){
		y_v <- Y[,v]
		s0_v <- template_mean[,v]
		E_v_inv <- diag(1/template_var[,v])
		Sigma_s_v <- solve(E_v_inv + At_nu0Cinv_A)
		miu_s[,v] <- Sigma_s_v	%*% (At_nu0Cinv %*% y_v + E_v_inv %*% s0_v) #Qx1
		var_s[,v] <- diag(Sigma_s_v)
	}

	result <- list(miu_s, var_s, theta, success)
	names(result) <- c('subjICmean', 'subjICvar', 'theta_MLE', 'success_flag')
	return(result)
}


####################################################################################
####################################################################################
### UpdateTheta() -- spatial ICA EM internal function to update parameters
####################################################################################
####################################################################################

UpdateTheta = function(template_mean, template_var, spde, BOLD, theta, C_diag, verbose=FALSE){

## ###################### OUTPUT ######################
##
##  List with the following elements:
##
##  theta_new (list) - final parameter estimates
##  theta_new$A             : (QxQ) mixing matrix
##  theta_new$nu0_sq        : (1x1) residual variance from first level
##  theta_new$kappa		 : (Qx1) SPDE smoothness parameter for each IC map
##
##  ICmean - estimates of subject-level ICs
##  ICvar - variance of subject-level ICs (for inference)


	Q = nrow(BOLD)
	V = ncol(BOLD)

	#initialize new objects
	theta_new = list(A = matrix(NA, Q, Q), nu0_sq = NA, kappa = rep(NA, Q))
	ICmean = matrix(0, Q, V) #subject IC mean
	ICvar = array(0, dim=c(Q,Q,V)) #subject IC var
	A_part1 = A_part2 = matrix(0, Q, Q) #two parts of product for A-hat (construct each looping over voxels)

	A = theta$A
	nu0_sq = theta$nu0_sq
	nu0C_inv = diag(1/(C_diag*nu0_sq))
	At_nu0Cinv = t(A) %*% nu0C_inv
	At_nu0Cinv_A = At_nu0Cinv %*% A

	print('Updating A')

	#store posterior moments for M-step of nu0_sq
	miu_s = matrix(NA, nrow=Q, ncol=V)
	miu_ssT = array(NA, dim=c(Q, Q, V))

	for(v in 1:V){

		y_v = BOLD[,v]
		s0_v = template_mean[,v]

		##########################################
		### E-STEP FOR A AND nu0^2: POSTERIOR MOMENTS OF s_i(v)
		##########################################

		E_v_inv = diag(1/template_var[,v])
		Sigma_s_v = solve(E_v_inv + At_nu0Cinv_A)
		miu_s_v = Sigma_s_v	%*% (At_nu0Cinv %*% y_v + E_v_inv %*% s0_v) #Qx1
		miu_ssT_v = (miu_s_v %*% t(miu_s_v)) + Sigma_s_v #QxQ
		miu_s[,v] = miu_s_v #save for M-step of nu0_sq
		miu_ssT[,,v] = miu_ssT_v #save for M-step of nu0_sq

		##########################################
		### M-STEP FOR A: CONSTRUCT PARAMETER ESTIMATES
		##########################################

		A_part1 = A_part1 + y_v %*% t(miu_s_v) #QxQ
		A_part2 = A_part2 + miu_ssT_v #QxQ

	}

	A_hat = orthonorm(A_part1 %*% solve(A_part2))

	##########################################
	### M-STEP FOR nu0^2: CONSTRUCT PARAMETER ESTIMATES
	##########################################

	print('Updating Error Variance nu0_sq')

	Cinv = diag(1/C_diag)
	Cinv_A = Cinv %*% A_hat
	At_Cinv_A = t(A_hat) %*% Cinv %*% A_hat
	nu0sq_part1 = nu0sq_part2 = nu0sq_part3 = 0

	for(v in 1:V){

		y_v = BOLD[,v]

		nu0sq_part1 = nu0sq_part1 + t(y_v) %*% Cinv %*% y_v
		nu0sq_part2 = nu0sq_part2 + t(y_v) %*% Cinv_A %*% miu_s[,v]
		nu0sq_part3 = nu0sq_part3 + sum(diag(At_Cinv_A %*% miu_ssT[,,v]))

	}

	nu0sq_hat = 1/(Q*V)*(nu0sq_part1 - 2*nu0sq_part2 + nu0sq_part3)

	##########################################
	### E-STEP for kappa_q: SECOND POSTERIOR MOMENT OF delta_i
	##########################################

	print('Updating SPDE Parameters kappa_q')

	#SPDE matrices, needed to construct R_q_inv
	F = spde$param.inla$M0
	G = 1/2*(spde$param.inla$M1 + t(spde$param.inla$M1))
	GFinvG = spde$param.inla$M2 #confirmed that this equals G %*% solve(F) %*% G

	#set up Sigma (QVxQV) as a sparse block diagonal matrix
	C1 = 1/(4*pi)
	Sigma_inv_list = vector('list', Q)
	for(q in 1:Q){
		D_q_inv = Diagonal(V, 1/sqrt(template_var[q,])) #sparse diagonal matrix
		kappa_q = theta$kappa[q]
		R_q_inv = C1 * (kappa_q^2 * F + 2 * G + kappa_q^(-2) * GFinvG)
		Sigma_q_inv = D_q_inv %*% R_q_inv %*% D_q_inv
		Sigma_inv_list[[q]] = Sigma_q_inv
	}
	Sigma_inv = bdiag(Sigma_inv_list)

	#set up P as a sparse matrix (see OneNote for illustration of this)
	cols = 1:(Q*V)
	rows_P1 = seq(1, (Q-1)*V+1, by=V)
	offset = rep(0:(V-1), each=Q)
	rows = rep(rows_P1, V) + offset
	P = sparseMatrix(i = rows, j = cols)

	#set up A, C and B = A' Cinv A as a sparse matrix
	ones = Diagonal(V)
	bigA = kronecker(ones, A)
	bigAC = kronecker(ones, t(A) %*% nu0C_inv)
	bigB = kronecker(ones, t(A) %*% nu0C_inv %*% A)

	# COMPUTE TERMS INVOLVING POSTERIOR MEAN AND VARIANCE OF DELTA_iq

	cov_delta_inv = t(P) %*% bigB %*% P + Sigma_inv

	# Two parts of Trace() term of MLE for kappa_q's
	# Part1 = Tr(mu_q' * K_q * mu_q), mu_q = T_q * cov_delta * P' * m --> use inla.qsolve to calculate cov_delta * P' * m
	# Part2 = Tr(K_q * T_q * cov_delta * T_q') --> use inla.qsample to draw from N(0, cov_delta), then estimate necessary elements of cov_delta

	### Part 1 (solve system of linear equations)

	print('Computing part 1 of trace terms involving big matrix inverse')

	t0 <- Sys.time()
	yvec <- as.vector(BOLD) # [y(1),...,y(V)]
	s0vec <- as.vector(template_mean) # [s0(1),...,s0(V)]
	Pmvec <- t(P) %*% bigAC %*% (yvec - bigA %*% s0vec)
	inla.setOption(smtp="pardiso")
	#1 minute or less using pardiso!!
	cov_delta_Pm <- inla.qsolve(Q = cov_delta_inv, B=matrix(Pmvec, ncol=1), method='solve')
	Trace1_part1 <- rep(NA, Q) #Trace1 = Tr(Dq_inv F Dq_inv E[delta_delta_q])
	Trace2_part1 <- rep(NA, Q) #Trace2 = Tr(Dq_inv G F_inv G Dq_inv E[delta_delta_q])
	#3 seconds!
	for(q in 1:Q){

		# K_q = sparse matrix appearing in trace, D_q^(-1) * GFinvG * D_q^(-1)
		D_q_inv = Diagonal(V, 1/sqrt(template_var[q,]))
		K_q1 <- D_q_inv %*% F %*% D_q_inv
		K_q2 <- D_q_inv %*% GFinvG %*% D_q_inv

		# T_q = matrix that selects the qth block of size V from a matrix with V*Q rows
		e_q <- Matrix(0, nrow=1, ncol=Q, sparse=TRUE)
		e_q[1,q] <- 1
		T_q <- kronecker(e_q, ones)

		#compute mu_q = T_q * cov_delta * P' * m and mu_q' * K_q * mu_q (part 1 of trace)
		mu_q <- T_q %*% cov_delta_Pm
		Trace1_part1[q] <- t(mu_q) %*% K_q1 %*% mu_q
		Trace2_part1[q] <- t(mu_q) %*% K_q2 %*% mu_q
	}
	print(Sys.time() - t0)

	### Part 2 (Monte Carlo)

	nsamp <- 1000 #number of samples
	print('Estimating part 2 of trace terms involving big matrix inverse')

	print(paste0('Drawing ',nsamp,' Monte Carlo samples'))
	#1-5 min
	print(system.time(musamp <- inla.qsample(nsamp, Q = cov_delta_inv, b=rep(0, V*Q), mu=rep(0, V*Q), num.threads=4))) #1.5 min for 100 samples, 5 min for 1000 samples

	#estimate diagonal blocks of cov_delta
	print('Computing trace of each block diagonal')
	t0 <- Sys.time()
	Trace1_part2 <- rep(NA, Q)
	Trace2_part2 <- rep(NA, Q)
	bigX_left <- KhatriRao(diag(1, V), matrix(1, nrow=nsamp, ncol=V)) # -- 13 SEC (do one time instead of KhatriRao(diag(1, V), Xctr_q) for each q)
	mat_nonzero <- as.matrix(1*(GFinvG != 0))
	diag(mat_nonzero) <- 1 #make sure all diagonal elements estimated (required for Trace 1)
	bigX_right <- KhatriRao(mat_nonzero, matrix(1, nrow=nsamp, ncol=V)) # -- 160 SEC (do one time instead of KhatriRao(as.matrix(1*(K_q != 0)), Xctr_q) for each q)
	bigX_right_cols <- which(GFinvG != 0, arr.ind = TRUE)[,2] #column indices of non-zero locations
	inds_left <- which(bigX_left != 0, arr.ind=TRUE)
	inds_right <- which(bigX_right != 0, arr.ind=TRUE)
	for(q in 1:Q){
		if(verbose) print(paste('Block ',q,' of ',Q))
		inds_q <- (1:V) + (q-1)*V
		D_q_inv = Diagonal(V, 1/sqrt(template_var[q,]))
		K_q1 <- D_q_inv %*% F %*% D_q_inv
		K_q2 <- D_q_inv %*% GFinvG %*% D_q_inv

		Xctr_q <- scale(t(musamp[inds_q,]), scale=FALSE) # < 1 sec

		# #compute cov_delta_qq, set unnecessary terms to zero -- 33 SEC
		# cov_delta_qq <- crossprod(Xctr_q)/(nsamp-1)
		# nonzero_q <- which(K_q != 0)
		# cov_delta_qq_zeros <- K_q
		# cov_delta_qq_zeros[nonzero_q] <- cov_delta_qq[nonzero_q]

		# #compute cov_delta_qq (necessary entries only!) -- 13 SEC

		#left-multiplication matrix
		vals_left <- as.vector(Xctr_q)
		bigX_left_q <- sparseMatrix(i = inds_left[,1], j = inds_left[,2], x = vals_left)
		#right-multiplication matrix
		X_repcols <- Xctr_q[,bigX_right_cols]
		vals_right <- as.vector(X_repcols)
		bigX_right_q <- sparseMatrix(i = inds_right[,1], j = inds_right[,2], x = vals_right)
		#multiply together
		cov_delta_qq_sparse <- t(bigX_left_q) %*% bigX_right_q / (nsamp - 1)

		#compute Trace for this part of the matrix
		Trace1_part2[q] <- sum(K_q1 * cov_delta_qq_sparse) #equivalent to sum(diag(K_q1 %*% cov_delta_qq)) and FAST (note: K_q symmetric)
		Trace2_part2[q] <- sum(K_q2 * cov_delta_qq_sparse) #equivalent to sum(diag(K_q2 %*% cov_delta_qq)) and FAST (note: K_q symmetric)
	}
	print(Sys.time() - t0)

	# NUMERICALLY ESTIMATE MLEs for kappa_q's -- 10 MIN

	print("Performing numerical optimization for kappa_q's")

	kappa_opt <- rep(NA, Q)
	t0 <- Sys.time()
	for(q in 1:Q){
		if(verbose) print(paste('Optimization ',q,' of ',Q))
		kappa_opt_q <- optimize(Q2_kappa_q, lower=-10, upper=10, maximum=TRUE,
								Fmat=F, Gmat=G, GFinvG=GFinvG,
								bigTrace1=Trace1_part1[q] + Trace1_part2[q],
								bigTrace2=Trace2_part1[q] + Trace2_part2[q])
		kappa_opt[q] <- exp(kappa_opt_q$maximum)
	}
	print(Sys.time() - t0)

	# RETURN NEW PARAMETER ESTIMATES

	theta_new$A <- A_hat
	theta_new$nu0_sq <- nu0sq_hat[1]
	theta_new$kappa <- kappa_opt
	return(theta_new)

}

####################################################################################
####################################################################################
### orthonorm() -- orthonormalizes a square, invertible matrix (X_orth = X * (X'X)^(-.5) )
####################################################################################
####################################################################################

orthonorm = function(X){

	#X is a square matrix to be orthonormalized

	#check that X is invertible (required for orthonorm(X) %*% t(orthonorm(X)) = I)
	if(!is.finite(determinant(X)$modulus)) error('X not invertible')

	#compute sqrt of (X'X)
	XtX_sqrt_inv = sqrt_XtX(X, inverse=TRUE)

	#perform orthogonalization
	result = X %*% XtX_sqrt_inv # symmetric orthogonalization
	if(!all.equal(Re(result), result)) error('Complex-valued result')
	return(result)

}

####################################################################################
####################################################################################
### sqrt_XtX() -- computes matrix square root of X'X
####################################################################################
####################################################################################

sqrt_XtX = function(X, inverse=FALSE){

	#X is a square matrix, X=UDV'
	#if inverse=TRUE, compute inverse of square root

	XtX = t(X) %*% X # X'X = V D^2 V'
	e = eigen(XtX)
	Vmat = e$vectors
	d2 = e$values #diagonal elements of D^2

	if(inverse) {
		if(!is.finite(determinant(X)$modulus)) error('X not invertible')
		result = Vmat %*% diag(sqrt(1/d2)) %*% t(Vmat)
	} else {
		result = Vmat %*% diag(sqrt(d2)) %*% t(Vmat)
	}

	return(result)

}



####################################################################################
####################################################################################
### Q2_kappa_q() -- computes part of log-likelihood involving kappa_q, for numerical optimization
####################################################################################
####################################################################################

Q2_kappa_q <- function(logkappa_q, Fmat, Gmat, GFinvG=NULL, bigTrace1, bigTrace2, C1 = 1/(4*pi)){
  # Fmat is a diagonal matrix appearing in SPDE precision
  # Gmat is a sparse nbhd matrix appearing in SPDE precision
  # GFinvG = Gmat %*% (1/Fmat) %*% Gmat (available directly through SPDE object)
  # SD_q is a vector of standard deviations from empirical prior
  # bigTrace1 is Trace(Dq_inv * F * Dq_inv * E[delta_q*delta_q'])
  # bigTrace2 is Trace(Dq_inv * G * F_inv * G * Dq_inv * E[delta_q*delta_q'])
  # C1 : for the unit variance case, tau^2 = C1/kappa^2.  C1 = 1/(4*pi) when alpha=2, nu=1, d=2

  require(Matrix)

  if(is.null(GFinvG)) GFinvG <- Gmat %*% solve(Fmat) %*% Gmat

  kappa_q <- exp(logkappa_q)

  #log determinant part
  part1_mat <- kappa_q^2 * Fmat + 2 * Gmat + kappa_q^(-2) * GFinvG
  logdet <- as.numeric(determinant(part1_mat)) #on log scale
  if(logdet[2] == -1) warning('negative determinant of precision matrix')
  part1 <- logdet[1]

  #first trace term in Q2
  part2a <- C1*(kappa_q^2) * bigTrace1

  #second trace term in Q2
  part2b <- C1*(kappa_q^(-2)) * bigTrace2

  result <- part1 - part2a - part2b
  return(result)

  # part1 decreases with kappa
  # part2=part2a+part2b also decreases with kappa
  # maximum depends on which part dominates
  # return(list(part1, - part2a - part2b, part1 - part2a - part2b))

}

















