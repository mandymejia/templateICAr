#' Diagnostic ICA
#'
#' Perform diagnostic independent component analysis (ICA) using expectation-maximization (EM)
#'
#' @param template_mean (A list of G matrices, each \eqn{VxL}) template mean
#'  estimates for each group \eqn{1} to \eqn{G}.
#' @param template_var (A list of G matrices, each \eqn{VxL}) template variance
#'  estimates for each group \eqn{1} to \eqn{G}.
#' @param BOLD (\eqn{VxT} matrix) BOLD fMRI data matrix, where \eqn{T} is the number
#'  of volumes (time points) and \eqn{V} is the number of brain locations
#' @param scale Should BOLD data be scaled by the spatial standard deviation
#'  before model fitting? Default: \code{TRUE}. If done when estimating
#'  templates, should be done here too.
#' @param meshes Either \code{NULL} (assume spatial independence) or a list of objects
#'  of type \code{templateICA_mesh} created by \code{make_mesh} (each list element
#'  corresponds to one brain structure)
#' @param Q2 The number of nuisance ICs to identify. If \code{NULL} (default),
#'  will be estimated. Only provide \eqn{Q2} or \eqn{maxQ} but not both.
#' @param maxQ Maximum number of ICs (template+nuisance) to identify
#'  (\eqn{L <= maxQ <= T}). Only provide \eqn{Q2} or \eqn{maxQ} but not both.
#' @param maxiter Maximum number of EM iterations. Default: 100.
#' @param epsilon Smallest proportion change between iterations. Default: 0.01.
#' @param verbose If \code{TRUE} (default), display progress of algorithm.
#' @param kappa_init Starting value for kappa.  If \code{NULL}, starting value
#'  will be determined automatically. Default: 0.4.
#'
#' @importFrom INLA inla inla.spde.result inla.pardiso.check inla.setOption
#' @importFrom pesel pesel
#' @importFrom stats optim
#' @importFrom abind abind
#' @importFrom matrixStats rowVars
#'
#' @return A list containing the posterior probabilities of group membership,
#'  the estimated independent components \strong{S} (a \eqn{VxL} matrix), their mixing matrix
#'  \strong{A} (a \eqn{TxL} matrix), the number of nuisance ICs estimated (Q_nuis)
#'
#' @export
#'
diagnosticICA <- function(template_mean,
                        template_var,
                        BOLD,
                        scale=TRUE,
                        meshes=NULL,
                        Q2 = NULL,
                        maxQ=NULL,
                        maxiter=100,
                        epsilon=0.01,
                        verbose=TRUE,
                        kappa_init=0.5){

  do_spatial <- !is.null(meshes)

  if(do_spatial){
    flag <- inla.pardiso.check()
    if(grepl('FAILURE',flag)) stop('PARDISO IS NOT INSTALLED OR NOT WORKING. PARDISO for R-INLA is required for computational efficiency. If you already have a PARDISO / R-INLA License, run inla.setOption(pardiso.license = "/path/to/license") and try again.  If not, run inla.pardiso() to obtain a license.')
    inla.setOption(smtp='pardiso')
  }

  if(!is.list(template_mean)) stop('template_mean must be a list')
  if(!is.list(template_var)) stop('template_var must be a list')
  if(length(template_mean) != length(template_var)) stop('template_mean and template_var must have the same length')
  G <- length(template_mean)

  if(!is.matrix(BOLD)) stop('BOLD must be a matrix')
  ntime <- ncol(BOLD) #length of timeseries
  nvox <- nrow(BOLD) #number of data locations
  L <- ncol(template_mean[[1]]) #number of ICs

  #check dimensionality of mesh projection matrices
  if(do_spatial){
    nmesh <- nvox2 <- 0
    for(k in seq_len(length(meshes))){
      Amat <- meshes[[k]]$A # n_orig x n_mesh matrix
      nmesh <- nmesh + ncol(Amat)
      nvox2 <- nvox2 + nrow(Amat)
    }
    if(nvox2 != nvox) stop('Mesh projection matrices (mesh$A) must have a total of nvox rows (nvox is the number of data locations, the columns of BOLD, template_mean and template_var)')
  }


  if(verbose) cat(paste0('Length of timeseries: T = ',ntime,'\n'))
  if(verbose) cat(paste0('Number of voxels/vertices: V = ',nvox,'\n'))
  if(verbose) cat(paste0('Number of ICs: L = ',L,'\n'))
  if(verbose) cat(paste0('Number of groups: G = ',G,'\n'))

  #check that each element of template_mean and template_var is a matrix
  #check that the dimensions of template_mean and template_var elements are ok
  if(any(sapply(template_mean, is.matrix)==FALSE)) stop('Each element of template_mean must be an VxL matrix')
  if(any(sapply(template_var, is.matrix)==FALSE)) stop('Each element of template_var must be an VxL matrix')
  if(any(sapply(template_mean, dim)[2,] != L) | any(sapply(template_mean, dim)[1,] != nvox)) stop('Each element of template_mean must be an VxL matrix')
  if(any(sapply(template_var, dim)[2,] != L) | any(sapply(template_var, dim)[1,] != nvox)) stop('Each element of template_mean must be an VxL matrix')

  #check that the number of data locations (nvox), time points (ntime) and ICs (L) makes sense
  if(ntime > nvox) warning('More time points than voxels. Are you sure?')
  if(L > nvox) stop('The arguments you supplied suggest that you want to estimate more ICs than you have data locations.  Please check the orientation and size of template_mean, template_var and BOLD.')
  if(L > ntime) stop('The arguments you supplied suggest that you want to estimate more ICs than you have time points.  Please check the orientation and size of template_mean, template_var and BOLD.')

  #check that the supplied mesh object is of type templateICA_mesh
  if(do_spatial){
    if(verbose) cat('Fitting a spatial model based on the meshes provided.  Note that computation time and memory demands may be high.')
    if(any(sapply(meshes, class) != 'templateICA_mesh')) stop('Each element of meshes argument should be of class templateICA_mesh. See help(make_mesh).')
  }

  if(round(maxiter) != maxiter | maxiter <= 0) stop('maxiter must be a positive integer')

  #check that maxQ makes sense
  if(!is.null(maxQ)){ if(round(maxQ) != maxQ | maxQ <= 0) stop('maxQ must be NULL or a round positive number') }
  if(is.null(maxQ)) maxQ <- ntime
  if(maxQ < L){
    warning('maxQ must be at least L.  Setting maxQ=L.')
    maxQ <- L
  }
  if(maxQ > ntime){
    warning('maxQ must be no more than T.  Setting maxQ = T.')
    maxQ <- ntime
  }

  if(class(scale) != 'logical' | length(scale) != 1) stop('scale must be a logical value')


  ### IDENTIFY AND REMOVE ANY BAD VOXELS/VERTICES

  keep <- rep(TRUE, nvox)
  for(g in 1:G){
    keep[rowSums(is.nan(template_mean[[g]])) > 0] <- FALSE
    keep[rowSums(is.na(template_mean[[g]])) > 0] <- FALSE
    keep[rowSums(is.nan(template_var[[g]])) > 0] <- FALSE
    keep[rowSums(is.na(template_var[[g]])) > 0] <- FALSE
  }
  keep[rowSums(is.nan(BOLD)) > 0] <- FALSE
  keep[rowSums(is.na(BOLD)) > 0] <- FALSE
  keep[rowVars(BOLD) == 0] <- FALSE

  if(sum(!keep) > 0){
    stop('Bad voxels detected. Check templates and BOLD for NA or NaN values, and check BOLD for zero-variance voxel timeseries. Set inds_data in make_mesh to exclude bad voxels.')
    # template_mean_orig <- template_mean
    # template_var_orig <- template_var
    # nvox <- sum(keep)
    # if(verbose) cat(paste0('Excluding ',sum(!keep),' bad (NA, NaN or flat) voxels/vertices from analysis.\n'))
    # for(g in 1:G){
    #   template_mean[[g]] <- template_mean[[g]][keep,]
    #   template_var[[g]] <- template_var[[g]][keep,]
    # }
    # BOLD <- BOLD[keep,]
  }


  ### 1. ESTIMATE AND DEAL WITH NUISANCE ICS (unless maxQ = L)

  template_mean_avg <- apply(abind(template_mean, along=3), c(1,2), mean)

  if(maxQ > L){

    #i. PERFORM DUAL REGRESSION TO GET INITIAL ESTIMATE OF TEMPLATE ICS
    BOLD1 <- scale_BOLD(BOLD, scale=scale)
    DR1 <- dual_reg(BOLD1, template_mean_avg)

    #ii. SUBTRACT THOSE ESTIMATES FROM THE ORIGINAL DATA --> BOLD2
    fit <- t(DR1$S) %*% t(DR1$A)
    BOLD2 <- BOLD1 - fit #data without template ICs

    #iii. ESTIMATE THE NUMBER OF REMAINING ICS
    #pesel function expects nxp data and will determine asymptotic framework
    #here, we consider n=T (volumes) and p=V (vertices), and will use p-asymptotic framework
    if(is.null(Q2)){
      if(verbose) cat(paste0('DETERMINING NUMBER OF NUISANCE COMPONENTS.... '))
      pesel_BOLD2 <- pesel(BOLD2, npc.max=maxQ-L, method='homo')
      Q2 <- pesel_BOLD2$nPCs #estimated number of nuisance ICs
      if(verbose) cat(paste0(Q2,'\n'))
    }

    #iv. ESTIMATE THE NUISANCE ICS USING GIFT/INFOMAX
    #if(verbose) cat(paste0('ESTIMATING AND REMOVING ',Q2,' NUISANCE COMPONENTS\n'))
    #ICA_BOLD2 <- icaimax(BOLD2, nc=Q2, center=TRUE)
    #fit <- ICA_BOLD2$M %*% t(ICA_BOLD2$S)

    #iv. INSTEAD OF ESTIMATING ICS, JUST ESTIMATE PCS!
    #THE RESIDUAL (BOLD3) IS THE EXACT SAME BECAUSE THE ICS ARE JUST A ROTATION OF THE PCS
    #IF THE NUISANCE ICS ARE NOT OF INTEREST, CAN TAKE THIS APPROACH
    svd_BOLD2 <- svd(t(BOLD2) %*% BOLD2, nu=Q2, nv=0)
    vmat <- diag(1/svd_BOLD2$d[1:Q2]) %*% t(svd_BOLD2$u) %*% t(BOLD2)
    fit <- svd_BOLD2$u %*% diag(svd_BOLD2$d[1:Q2]) %*% vmat

    #v. SUBTRACT THOSE ESTIMATES FROM THE ORIGINAL DATA --> BOLD3
    BOLD3 <- BOLD1 - t(fit) #original data without nuisance ICs

  } else {

    # USE ORIGINAL DATA, SCALED, SINCE WE ARE ASSUMING NO NUISANCE COMPONENTS
    BOLD3 <- scale_BOLD(BOLD, scale=scale) #center, and if scale=TRUE, scale

  }

  ### 2. PERFORM DIMENSION REDUCTION --> BOLD4

  dat_list <- dim_reduce(BOLD3, L)

  BOLD4 <- dat_list$data_reduced
  H <- dat_list$H
  Hinv <- dat_list$H_inv
  C_diag <- dat_list$C_diag
  C_diag <- C_diag * (dat_list$sigma_sq) #(nu^2)HH' in paper

  ### 3. SET INITIAL VALUES FOR PARAMETERS

  #initialize mixing matrix (use dual regression-based estimate for starting value)
  dat_DR <- dual_reg(BOLD3, template_mean_avg)
  HA <- H %*% dat_DR$A #apply dimension reduction
  theta0 <- list(A = HA)
  theta00 <- theta0
  theta00$nu0_sq = dat_list$sigma_sq #required for template ICA

  #temporary
  S_DR <- t(dat_DR$S)
  for(g in 1:G) {
    num_smallvar <- sum(template_var[[g]] < 1e-6)
    if(num_smallvar>0){
      #if(verbose) cat(paste0('Setting ',num_smallvar,' (',round(100*num_smallvar/length(template_var[[g]]),1),'%) very small variance values in group ',g,' template to ',1e-6,'.\n'))
      template_var[[g]][template_var[[g]] < 1e-6] = 1e-6 #to prevent problems when inverting covariance
    }
  }
  template_var_max <- template_var[[1]]
  for(g in 2:G){
    template_var_diff <- template_var[[g]] - template_var_max
    template_var_diff[template_var_diff < 0] <- 0 #places where max is already greater, do not change
    template_var_max <- template_var_max + template_var_diff #places where max is not greater, add to max
  }
  template_var_max[template_var_max < 1e-6] <- 1e-6


  #print('Distance using Group Templates')
  #dist1 <- colSums((S_DR - template_mean[[1]])^2/(template_var[[1]]))
  #dist2 <- colSums((S_DR - template_mean[[2]])^2/(template_var[[2]]))
  #print(cbind(dist1, dist2, dist1<dist2))
  #print(mean(dist1<dist2))
  #print(sum(dist1) < sum(dist2))

  dist1 <- colSums((S_DR - template_mean[[1]])^2/(template_var_max))
  dist2 <- colSums((S_DR - template_mean[[2]])^2/(template_var_max))
  #print(cbind(dist1, dist2, dist1<dist2))
  #print(mean(dist1<dist2))
  #print(sum(dist1) < sum(dist2))
  #for 2 groups only:
  pr1 <- mean(dist1<dist2)
  pr2 <- 1-pr1
  pr_z = c(pr1, pr2)
  print(paste0('Initial Group Probabilities: ',paste(round(pr_z,3), collapse=', ')))

  #pr_z <- c(0.5, 0.5)

  #use pr_z as the starting value to get the first set of estimates of s
  #to do that, could run the Update function with returnMAP=TRUE
  #once we have a set of s estimates, use those to update the probability in the same way
  #the probability is one of the parameters, so we should run the model till it and other parameters converge

  theta0$pr_z <- pr_z

  #end temporary code

  ### 4. RUN EM ALGORITHM!

  #NON-SPATIAL DIAGNOSTIC ICA

  if(verbose) cat('INITIATING WITH STANDARD DIAGNOSTIC ICA\n')
  resultEM <- EM_diagnosticICA.independent(template_mean,
                                           template_var,
                                           BOLD=BOLD4,
                                           theta0=theta0,
                                           C_diag=C_diag,
                                           maxiter=maxiter,
                                           epsilon=epsilon,
                                           verbose=verbose)
  resultEM$A <- Hinv %*% resultEM$theta_MLE$A
  class(resultEM) <- 'dICA'

  #SPATIAL DIAGNOSTIC ICA

  if(do_spatial){

    # template_var1_avg <- apply(abind(template_var, along=3), c(1,2), mean) #avg within-group var
    # template_var2_avg <- apply(abind(template_mean, along=3), c(1,2), var) #between-group var (the unbiased N-1 version), equals 1/2*(y1-y2)^2 for N=2
    # template_var_avg <- template_var1_avg + template_var2_avg #between-subject var (pooled groups)
    # resultEM <- EM_templateICA.independent(template_mean_avg, template_var_avg, BOLD4, theta00, C_diag, maxiter=maxiter, epsilon=epsilon, verbose=verbose)
    # class(resultEM) <- 'tICA'

    resultEM_dICA <- resultEM
    # which_group <- which.max(resultEM$group_probs) # never used

    ## TO DO: Delete this or generalize to multiple meshes

    # #starting value for kappas
    # if(is.null(kappa_init)){
    #   if(verbose) cat('Using ML on initial estimates to determine starting value for kappa\n')
    #
    #   #organize data and replicates
    #   for(q in 1:L){
    #     #print(paste0('IC ',q,' of ',L))
    #     d_q <- resultEM$subjICmean[,q] - template_mean[[which_group]][,q]
    #     rep_q <- rep(q, length(d_q))
    #     D_diag_q <- sqrt(template_var[[which_group]][,q])
    #     if(q==1) {
    #       dev <- d_q
    #       rep <- rep_q
    #       D_diag <- D_diag_q
    #     } else {
    #       dev <- c(dev, d_q)
    #       rep <- c(rep, rep_q)
    #       D_diag <- c(D_diag, D_diag_q)
    #     }
    #   }
    #
    #   #determine MLE of kappa
    #   # print(system.time(opt <- optim(par=c(0,-20),
    #   #                fn=loglik_kappa_est,
    #   #                method='L-BFGS-B',
    #   #                lower=c(-5,-20),
    #   #                upper=c(1,Inf), #kappa usually less than 1, log(1)=0
    #   #                delta=dev,
    #   #                D_diag=D_diag,
    #   #                mesh=mesh,
    #   #                Q=L)))
    #   print(system.time(opt <- optim(par=0,
    #                  fn=loglik_kappa_est,
    #                  method='L-BFGS-B',
    #                  lower=-5,
    #                  upper=1, #kappa usually less than 1, log(1)=0. exp(1) = 2.7, which is much larger than kappa in most cases
    #                  log_var = -15, #optimizing over kappa and residual var typically results in very small estimates of the latter
    #                  delta=dev,
    #                  D_diag=D_diag,
    #                  mesh=mesh,
    #                  Q=L)))
    #   kappa_init <- exp(opt$par[1])
    #
    #   if(verbose) print(paste0('Starting value for kappa = ',paste(round(kappa_init,3), collapse=', ')))
    # }

    if(kappa_init <= 0 | !is.numeric(kappa_init)) stop('Value of kappa_init must be greater than zero')
    theta0$kappa <- rep(kappa_init, L)


 # resultEM$A <- Hinv %*% resultEM$theta_MLE$A

    if(verbose) cat('ESTIMATING SPATIAL MODEL ... ')
    t000 <- Sys.time()
    resultEM <- EM_diagnosticICA.spatial(template_mean=template_mean,
                                         template_var=template_var,
                                         meshes=meshes,
                                         BOLD=BOLD4,
                                         theta0=theta0,
                                         C_diag=C_diag,
                                         maxiter=maxiter,
                                         epsilon=epsilon,
                                         verbose=verbose)
    if(verbose) cat(paste0(round((Sys.time() - t000)/60,2),' minutes\n'))

    #organize estimates and variances in matrix form
    resultEM$subjICmean <- matrix(resultEM$subjICmean, ncol=L)
    resultEM$subjICvar <- matrix(diag(resultEM$subjICcov), ncol=L)
  }

  resultEM$A <- Hinv %*% resultEM$theta_MLE$A

  #for stICA, return tICA estimates also
  if(do_spatial){
    resultEM$result_dICA <- resultEM_dICA
    class(resultEM) <- 'sdICA'
  }

  return(resultEM)

}
