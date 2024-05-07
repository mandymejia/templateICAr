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
#' @param Wmat,Gamma_inv Matrices used for modeling temporal correlations in A
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
  Wmat, Gamma_inv,
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
  #mu_tau2 <- apply(BOLD - t(mu_A %*% mu_S),1,var) #Vx1
  mu_tau2 <- mean((BOLD - t(mu_A %*% mu_S))^2) #Vx1
  mu_G <- template_FC$psi / (template_FC$nu - nICs - 1)
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

    #this function will return mu_A and E[A'A], which is needed for tau2 and for S
    #will also return mu_{A-tilde} and E[A-tilde'A-tilde]

    #only ~5s for Q=25, V=90k, T=1200!
    A_new <- update_A(mu_tau2, mu_S, cov_S, mu_G, Gamma_inv, Wmat, #mu_alpha,
                      BOLD, ntime, nICs, nvox)

    mu_A_old <- mu_A
    mu_A <- A_new$mu_A
    mu_A_tilde <- A_new$mu_A_tilde
    change_A <- mean(abs(mu_A - mu_A_old)/(abs(mu_A_old)+0.1)) #add 0.1 to denominator to avoid dividing by zero

    E_AtA <- A_new$E_AtA
    E_AtA_tilde <- A_new$E_AtA_tilde

    #b. UPDATE S
    S_new <- update_S(mu_tau2, mu_A, E_AtA, D_inv, D_inv_S, BOLD, ntime, nICs, nvox)
    change_S <- mean(abs(S_new[[1]] - mu_S)/(abs(mu_S)+0.1)) #add 0.1 to denominator to avoid dividing by zero
    #change_S <- sqrt(crossprod(c(S_new[[1]] - mu_S))) #same as in SQUAREM
    mu_S <- S_new$mu_S
    cov_S <- S_new$cov_S

    #c. UPDATE G
    if(method_FC == 'VB1'){
      G_new <- update_G_IW(mu_A_tilde, E_AtA_tilde, #mu_alpha, cov_alpha,
                        template_FC, ntime, nICs)
      psi_G <- G_new$psi1
      nu_G <- G_new$nu1
      mu_G_new <- psi_G / (nu_G - nICs - 1)
      change_G <- mean(abs(mu_G - mu_G_new)/(abs(mu_G)+0.1)) #add 0.1 to denominator to avoid dividing by zero
      #change_G <- sqrt(crossprod(c(mu_G_new - mu_G))) #same as in SQUAREM
      mu_G <- mu_G_new
      print(diag(mu_G)[1:5])
    } else {

    # [TO DO] implement Cholesky prior-based estimation of G

    }

    # HERE -- pass E_AtA and E_SSt to the next function instead of cov_A, cov_S

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
#' @param Gamma_inv,Wmat matrices used to account for temporal autocorrelation
#' @param BOLD (\eqn{T \times V} matrix) preprocessed fMRI data.
#' @param ntime,nICs,nvox Number of timepoints in data, number of ICs, and
#'  the number of data locations.
#'
#' @return List of length two: \code{mu_A} (TxQ) and \code{cov_A} (QxQ).
#'
#' @keywords internal
update_A <- function(
  mu_tau2, mu_S, cov_S, mu_G, Gamma_inv, Wmat,
  BOLD, ntime, nICs, nvox){

  #cov_A (TQxTQ) -- approximation method 1 (ignore Gamma)
  E_SSt <- apply(cov_S, 1:2, sum) + mu_S %*% t(mu_S) #sum over v to get E[SS']
  E_SSt_big <- bdiag_m2(mat = E_SSt, N = ntime) #I_T \otimes E[SS']
  G_inv <- solve(mu_G)
  Gamma_G_inv <- Matrix::kronecker(Gamma_inv, G_inv) # Gamma^(-1) \otimes G^(-1)
  cov_A_inv <- (1/mu_tau2)*E_SSt_big + Gamma_G_inv

  #mu_A (TQx1)
  tmp <- mu_S %*% t(BOLD) #equivalent to S_otimes %*% y in paper
  tmp <- c(tmp) #vectorize a QxT matrix into tmp_1,...,tmp_T
  tmp <- (1/mu_tau2) * tmp #now we just need to pre-multiply by Cov(a)
  mu_A <- Matrix::solve(a = cov_A_inv, b = tmp) #compute estimate of (a_1,...,a_T)
  mu_A <- matrix(mu_A, nrow=ntime, ncol=nICs, byrow=TRUE) #reshape to TxQ matrix

  ### Constrain each column of A to have var=1 and mean=0
  sd_A <- apply(mu_A, 2, sd)
  D_A <- diag(1/sd_A) #use this below to correct cov(A) for unit-var A
  mu_A <- scale(mu_A)

  ### [TO DO] Appropriately rescale Cov(a_t) and Cov(\tilde{a}_t)
  #cov_A <- diag(1/sd_A) %*% cov_A %*% diag(1/sd_A)

  ### Compute E[A-tilde]
  mu_A_tilde <- Wmat %*% mu_A

  ### APPROXIMATE E[A'A]

  # The idea: invert submatrices at the beginning, middle and end, then extrapolate

  inds_by_t <- matrix(1:(ntime*nICs), nrow=nICs, ncol=ntime) #indices of cov(A) corresponding to all T QxQ blocks
  ar_order <- 5
  half_width <- ar_order*2 #enough to capture effects of autocorrelation and avoid boundary effects (see Appendix of paper)
  #1. middle
  t0 <- round(ntime/2) #this is the block that will serve as proxy for all the middle blocks
  T0 <- (t0 - half_width):(t0 + half_width) #make a window of 10 time points around t0 to capture autocorrelation and avoid spurious boundary effects on t0 itself
  inds_T0 <- inds_by_t[,T0] #which rows/columns of cov(A) to invert
  cov_A0 <- Matrix::solve(cov_A_inv[c(inds_T0),c(inds_T0)]) #invert the middle submatrix
  inds_T0_renum <- (inds_T0 - min(inds_T0) + 1)
  inds_t0 <- inds_T0_renum[,(half_width+1)] #which indices of submatrix correspond to t0
  cov_t0 <- cov_A0[inds_t0,inds_t0] #INGREDIENT 1

  #2. beginning
  T1 <- 1:(2*half_width) #first contiguous segment of time points
  inds_T1 <- inds_by_t[,T1]
  cov_A1 <- Matrix::solve(cov_A_inv[c(inds_T1),c(inds_T1)]) #invert the beginning submatrix
  inds_T1_renum <- (inds_T1 - min(inds_T1) + 1)
  inds_t1 <- inds_T1_renum[,(1:half_width)] #which indices of submatrix to save
  cov_T1 <- array(NA, dim=c(nICs,nICs,half_width)) #only use the first "half_width" estimates
  for(kk in 1:half_width){
    inds_kk <- inds_t1[,kk]
    cov_T1[,,kk] <- as.matrix(cov_A1[inds_kk,inds_kk]) #kth diagonal block of submatrix cov_T1
  } #INGREDIENT 2

  #3. end
  T2 <- tail(1:ntime, half_width*2) #last contiguous segment of time points
  inds_T2 <- inds_by_t[,T2]
  cov_A2 <- Matrix::solve(cov_A_inv[c(inds_T2),c(inds_T2)]) #invert the last submatrix
  inds_T2_renum <- (inds_T2 - min(inds_T2) + 1)
  inds_t2 <- inds_T2_renum[,tail(1:(2*half_width), half_width)] #which indices of submatrix to save
  cov_T2 <- array(NA, dim=c(nICs,nICs,half_width)) #only use the last "half_width" estimates
  for(kk in 1:half_width){
    inds_kk <- inds_t2[,kk]
    cov_T2[,,kk] <- as.matrix(cov_A2[inds_kk,inds_kk]) #kth diagonal block of submatrix cov_T2
  } #INGREDIENT 3

  #sum up to approximate E[A'A]
  cov_A_sum <- apply(cov_T1, 1:2, sum) + apply(cov_T2, 1:2, sum) + (ntime - half_width*2)*cov_t0
  cov_A_sum <- D_A %*% cov_A_sum %*% D_A
  list1 <- list(cov_t0, cov_T1, cov_T2, cov_A_sum) #to return
  E_AtA <- as.matrix(cov_A_sum) + t(mu_A) %*% mu_A


  # APPROXIMATE E[A-tilde'A-tilde]

  Wmat_big <- kronecker(Wmat, Matrix(diag(nICs), sparse=TRUE))
  nA1 <- nrow(cov_A1)
  nA2 <- nrow(cov_A2)
  nWbig <- nrow(Wmat_big)
  Wmat_A1 <- Wmat_big[c(inds_T1),c(inds_T1)]
  Wmat_A2 <- Wmat_big[c(inds_T2),c(inds_T2)]
  Wmat_A0 <- Wmat_big[c(inds_T0),c(inds_T0)]
  sd_A_big <- kronecker(Matrix(diag(length(T1))), D_A, sparse=TRUE) #for standardizing var(A)
  sd_A_big0 <- kronecker(Matrix(diag(length(T0))), D_A, sparse=TRUE) #for standardizing var(A)
  cov_A1_tilde <- Wmat_A1 %*% sd_A_big %*% cov_A1 %*% sd_A_big %*% Wmat_A1
  cov_A2_tilde <- Wmat_A2 %*% sd_A_big %*% cov_A2 %*% sd_A_big %*% Wmat_A2
  cov_A0_tilde <- Wmat_A0 %*% sd_A_big0 %*% cov_A0 %*% sd_A_big0 %*% Wmat_A0
  cov_t0_tilde <- cov_A0_tilde[inds_t0,inds_t0] #INGREDIENT 1
  cov_T1_tilde <- array(NA, dim=c(nICs,nICs,half_width)) #only use the first "half_width" estimates
  for(kk in 1:half_width){
    inds_kk <- inds_t1[,kk]
    cov_T1_tilde[,,kk] <- as.matrix(cov_A1_tilde[inds_kk,inds_kk]) #kth diagonal block of submatrix cov_T1
  } #INGREDIENT 2
  cov_T2_tilde <- array(NA, dim=c(nICs,nICs,half_width)) #only use the last "half_width" estimates
  for(kk in 1:half_width){
    inds_kk <- inds_t2[,kk]
    cov_T2_tilde[,,kk] <- as.matrix(cov_A2_tilde[inds_kk,inds_kk]) #kth diagonal block of submatrix cov_T2
  } #INGREDIENT 3

  #sum up to approximate E[Atilde'Atilde]
  cov_A_tilde_sum <- apply(cov_T1_tilde, 1:2, sum) + apply(cov_T2_tilde, 1:2, sum) + (ntime - half_width*2)*cov_t0_tilde
  E_AtA_tilde <- as.matrix(cov_A_tilde_sum) + t(mu_A_tilde) %*% mu_A_tilde


  # TOO SLOW BUT VERY CLOSE TO THE ABOVE APPROACH
  # # The idea: Approximate V(a) as a block diagonal matrix, then compute V(a-tilde) = W %*% V(a) %*% W
  #
  # mid <- function(x, n){
  #   #grabbing middle n can be achieved by discarding the first 0.5*(length(x)  - n) elements and using head()
  #   discard <- round(0.5*(length(x)  - n))
  #   x <- x[-(1:discard)]
  #   head(x, n)
  # }
  # #a) start and end blocks
  # beginning <- cov_A1[c(inds_t1),c(inds_t1)] #this is for the FIRST "half_width" time points
  # end <- cov_A2[c(inds_t2),c(inds_t2)] #this is for the LAST "half_width" time points
  # #b) main middle blocks
  # t0_plus <- mid(T0, half_width+1) #which time points to save from middle segment
  # inds_t0_plus <- inds_T0_renum[,mid(1:length(T0), half_width+1)] #corresponding indices of submatrix
  # middle <- cov_A0[c(inds_t0_plus),c(inds_t0_plus)] #this is not the whole middle, just a segement for "half_width + 1" time points
  # #c) final "mini" middle block
  # ntime_middle <- ntime - half_width*2 #number of time points we need to fill in-between "beginning" and "end"
  # n_middle <- floor(ntime_middle/(half_width+1)) #the number of full-size "middle" blocks that fit
  # ntime_middle_mini <- ntime_middle - n_middle*(half_width+1) #length of last "mini" middle block
  # if(ntime_middle_mini > 0){
  #   t0_plus_mini <- mid(T0, ntime_middle_mini) #which time points to save from middle segment
  #   inds_t0_plus_mini <- inds_T0_renum[,mid(1:length(T0), ntime_middle_mini)] #corresponding indices of submatrix
  #   middle_mini <- cov_A0[c(inds_t0_plus_mini),c(inds_t0_plus_mini)] #this is not the whole middle, just a segement for "half_width + 1" time points
  # }
  # #d) put them together into a block diagonal matrix
  # #make a list of blocks for block diagonal matrix
  # nK <- 2 + n_middle #number of blocks
  # if(ntime_middle_mini > 0) nK <- nK + 1
  # cov_A_list <- vector('list', length=nK)
  # cov_A_list[[1]] <- beginning
  # cov_A_list[[nK]] <- end
  # for(kk in (1+(1:n_middle))){ cov_A_list[[kk]] <- middle }
  # if(ntime_middle_mini > 0) cov_A_list[[1 + n_middle + 1]] <- middle_mini
  # cov_A <- Matrix::bdiag(cov_A_list)
  # if(nrow(cov_A) != nrow(cov_A_inv)) warning('Dimensions of cov_A are not correct, check the code!')
  # cov_A_tilde <- Wmat_big %*% cov_A %*% Wmat_big
  # cov_A_tilde_sum <- matrix(0, nICs, nICs)
  # for(tt in 1:ntime){
  #   inds_tt <- inds_by_t[,tt]
  #   cov_A_tilde_sum <- cov_A_tilde_sum + as.matrix(cov_A_tilde[c(inds_tt),c(inds_tt)])
  # }
  # E_AtA_tilde <- cov_A_tilde_sum + t(mu_A_tilde) %*% mu_A_tilde

  # #check
  # diff(cov_T1[,,half_width], cov_t0) #SHOULD BE SIMILAR
  # diff(cov_T2[,,1], cov_t0) #SHOULD BE SIMILAR
  # diff(cov_T2[,,1], cov_T1[,,half_width]) #SHOULD BE SIMILAR

  # OLD CODE ASSUMING A TEMPORALLY INDEPENDENT
  # #cov_A (QxQ) -- common across t=1,...,T
  # tmp1 <- mu_S/sqrt(mu_tau2)
  # tmp1 <- tcrossprod(tmp1)
  # tmp2 <- cov_S/mu_tau2
  # tmp2 <- apply(tmp2, c(1,2), sum)
  # tmp <- tmp1 + tmp2
  # cov_A <- solve(tmp + G_inv)
  #
  # #mu_A
  # tmp1 <- BOLD/mu_tau2
  # tmp2 <- tmp1 %*% t(mu_S)
  # mu_A <- array(0, dim = c(ntime, nICs)) #TxQ
  # for(t in 1:ntime) mu_A[t,] <- cov_A %*% tmp2[t,]

  list(mu_A = mu_A,
       mu_A_tilde = mu_A_tilde,
       E_AtA = E_AtA,
       E_AtA_tilde = E_AtA_tilde)
}

#' Update S for VB FC Template ICA
#'
#' @param mu_tau2,mu_A,E_AtA, Most recent posterior estimates
#' @param D_inv,D_inv_S Some pre-computed quantities.
#' #' @param BOLD (\eqn{T \times V} matrix) preprocessed fMRI data.
#' @param ntime,nICs,nvox Number of timepoints in data, number of ICs, and
#'  the number of data locations.
#'
#' @return List of length two: \code{mu_S} (QxV) and \code{cov_S} (QxQxV).
#'
#' @keywords internal
update_S <- function(
  mu_tau2, mu_A, E_AtA,
  D_inv, D_inv_S,
  BOLD,
  ntime, nICs, nvox){

  for(v in 1:nvox){
    cov_S[,,v] <- solve((1/mu_tau2) * E_AtA + diag(D_inv[v,]))
  }

  mu_S <- array(0, dim = c(nICs, nvox)) #QxV
  tmp <- t(mu_A) %*% BOLD
  for(v in 1:nvox){
    #sum over t=1...T part
    #D_v_inv <- diag(1/template_var[v,])
    mu_S[,v] <- cov_S[,,v] %*% ( (1/mu_tau2) * tmp[,v] + D_inv_S[v,] )
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
    mu_A_tilde, E_AtA_tilde, template_FC, ntime, nICs){

  nu0 <- template_FC$nu
  psi0 <- template_FC$psi

  nu1 <- nu0 + ntime
  psi1 <- psi0 + E_AtA_tilde

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
    (1/2)*BOLD2 -
    sum(BOLD * mu_A %*% mu_S) +
    (1/2)*sum(diag(E_aat %*% E_sst))


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

##' Fast version of Matrix :: .bdiag() -- for the case of *many*  (k x k) matrices:
##' @param lmat list(<mat1>, <mat2>, ....., <mat_N>)  where each mat_j is a  k x k 'matrix'
##' @return a sparse (N*k x N*k) matrix of class  \code{"\linkS4class{dgCMatrix}"}.
bdiag_m <- function(lmat) {
  ## Copyright (C) 2016 Martin Maechler, ETH Zurich
  if(!length(lmat)) return(new("dgCMatrix"))
  stopifnot(is.list(lmat), is.matrix(lmat[[1]]),
            (k <- (d <- dim(lmat[[1]]))[1]) == d[2], # k x k
            all(vapply(lmat, dim, integer(2)) == k)) # all of them
  N <- length(lmat)
  if(N * k > .Machine$integer.max)
    stop("resulting matrix too large; would be  M x M, with M=", N*k)
  M <- as.integer(N * k)
  ## result: an   M x M  matrix
  new("dgCMatrix", Dim = c(M,M),
      ## 'i :' maybe there's a faster way (w/o matrix indexing), but elegant?
      i = as.vector(matrix(0L:(M-1L), nrow=k)[, rep(seq_len(N), each=k)]),
      p = k * 0L:M,
      x = as.double(unlist(lmat, recursive=FALSE, use.names=FALSE)))
}

##' Fast version of Matrix :: .bdiag() -- for the case of *many identical*  (k x k) matrices:
##' @param mat a  k x k 'matrix'
##' @param N how many times to repeat \code{mat}
##' @return a sparse (N*k x N*k) matrix of class  \code{"\linkS4class{dgCMatrix}"}.
bdiag_m2 <- function(mat, N) {
  k <- nrow(mat)
  if(N * k > .Machine$integer.max)
    stop("resulting matrix too large; would be  M x M, with M=", N*k)
  M <- as.integer(N * k)
  x <- rep(c(mat), N)
  ## result: an   M x M  matrix
  new("dgCMatrix", Dim = c(M,M),
      ## 'i :' maybe there's a faster way (w/o matrix indexing), but elegant?
      i = as.vector(matrix(0L:(M-1L), nrow=k)[, rep(seq_len(N), each=k)]),
      p = k * 0L:M,
      x = as.double(x))
}

