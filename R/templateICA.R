#' Template ICA
#'
#' Perform template independent component analysis (ICA) using expectation-maximization (EM)
#'
#' @param template_mean (VxL matrix) template mean estimates, i.e. mean of empirical population prior for each of L independent components
#' @param template_var (VxL matrix) template variance estimates, i.e. between-subject variance of empirical population prior for each of L ICs
#' @param BOLD (VxT matrix) BOLD fMRI data matrix, where T is the number of volumes (time points) and V is the number of brain locations. Or, a list of such data matrices. 
#' @param scale Logical indicating whether BOLD data should be scaled by the spatial standard deviation before model fitting. If done when estimating templates, should be done here too.
#' @param normA Normalize the A matrix (spatial maps)?
#' @param meshes Either NULL (assume spatial independence) or a list of objects of type \code{templateICA_mesh}
#' created by \code{make_mesh} (spatial priors are assumed on each independent component).
#' Each list element represents a brain structure, between which spatial independence is assumed (e.g. left and right hemispheres)
#' @param Q2 The number of nuisance ICs to identify. If \code{NULL}, will be estimated.
#'  Only provide \code{Q2} or \code{maxQ} but not both.
#' @param maxQ Maximum number of ICs (template+nuisance) to identify 
#'  (L <= maxQ <= T). If \code{maxQ == L}, then do not remove any nuisance regressors.
#'  Only provide \code{Q2} or \code{maxQ} but not both.
#' @param maxiter Maximum number of EM iterations
#' @param epsilon Smallest proportion change between iterations (e.g. .01)
#' @param verbose If \code{TRUE}. display progress of algorithm
# @param common_smoothness If \code{TRUE}. use the common smoothness version of the spatial template ICA model, which assumes that all IC's have the same smoothness parameter, \eqn{\kappa}
#' @param kappa_init Starting value for kappa.  Default: \code{0.2}.
#'
#' @return A list containing the estimated independent components S (a VxL matrix), their mixing matrix A (a TxL matrix), and the number of nuisance ICs estimated (Q_nuis)
#'
#' @export
#'
# @importFrom INLA inla inla.spde.result inla.pardiso.check inla.setOption
#' @importFrom stats optim
#' @importFrom matrixStats rowVars
#'
templateICA <- function(template_mean,
                        template_var,
                        BOLD,
                        scale=TRUE,
                        normA=FALSE,
                        meshes=NULL,
                        Q2=NULL,
                        maxQ=NULL,
                        maxiter=100,
                        epsilon=0.001,
                        verbose=TRUE,
                        #common_smoothness=TRUE,
                        kappa_init=0.2){

  #TO DO: Define do_FC, pass in FC template
  do_FC <- FALSE

  if (!is.null(meshes)) {
    INLA_check()
    flag <- INLA::inla.pardiso.check()
    if (grepl('FAILURE',flag)) {
      stop(
        'PARDISO IS NOT INSTALLED OR NOT WORKING. ', 
        'PARDISO for R-INLA is required for computational efficiency. ',
        'If you already have a PARDISO / R-INLA License, run inla.setOption(pardiso.license = "/path/to/license") and try again. ',
        'If not, run inla.pardiso() to obtain a license.'
      )
    }
    INLA::inla.setOption(smtp='pardiso')
  }

  # Handle if the inputs are xiftis
  template_mean <- as.matrix(template_mean)
  template_var <- as.matrix(template_var)

  # If `BOLD` is a list, ensure that all dimensions are the same.
  if (is.list(BOLD)) {
    if (length(BOLD) == 1) {
      BOLD <- BOLD[[1]]
    } else {
      stopifnot(all(vapply(BOLD, is.matrix, FALSE)))
      stopifnot(length(unique(vapply(BOLD, nrow, 0))) == 1)
    }
  }
  multi_scans <- is.list(BOLD)

  # Get data dimensions and the number of ICs
  if (multi_scans) {
    ntime <- vapply(BOLD, ncol, 0); nvox <- nrow(BOLD[[1]])
  } else {
    ntime <- ncol(BOLD); nvox <- nrow(BOLD)
  }
  L <- ncol(template_mean) #number of ICs

  #check that the number of data locations (nvox), time points (ntime) and ICs (L) makes sense
  if (sum(ntime) > nvox) warning('More time points than voxels. Are you sure?')
  if (L > nvox) stop('The arguments you supplied suggest that you want to estimate more ICs than you have data locations.  Please check the orientation and size of template_mean, template_var and BOLD.')
  if (L > sum(ntime)) stop('The arguments you supplied suggest that you want to estimate more ICs than you have time points.  Please check the orientation and size of template_mean, template_var and BOLD.')

  #check that all arguments have consistent number of data locations (nvox) and ICs (L)
  if (nrow(template_mean) != nvox | nrow(template_var) != nvox) stop('template_mean, template_var and BOLD must have same number of data locations, but they do not.')
  if (ncol(template_var) != L) stop('template_mean and template_var must have the same number of ICs (rows), but they do not.')

  #check that the supplied mesh object is of type templateICA_mesh
  do_spatial <- !is.null(meshes)
  if (do_spatial) {
    if (verbose) cat('Fitting a spatial model based on the mesh provided.  Note that computation time and memory demands may be high.\n')
    if (!is.list(meshes)) stop('meshes argument must be a list.')
    mesh_classes <- sapply(meshes, 'class')
    if (any(mesh_classes != 'templateICA_mesh')) stop('Each element of meshes argument should be of class templateICA_mesh. See help(make_mesh).')
    #if(class(common_smoothness) != 'logical' | length(common_smoothness) != 1) stop('common_smoothness must be a logical value')
  }
  #if(!do_spatial & !is.null(kappa_init)) stop('kappa_init should only be provided if mesh also provided for spatial modeling')

  if (round(maxiter) != maxiter | maxiter <= 0) stop('maxiter must be a positive integer')

  if (!is.null(kappa_init)) {
    if(length(kappa_init) != 1 | kappa_init <= 0) stop('kappa_init must be a positive scalar or NULL')
  }

  #check that maxQ makes sense
  if (!is.null(maxQ)) { 
    if(round(maxQ) != maxQ || maxQ <= 0) stop('maxQ must be NULL or a round positive number.')
  } else {
    maxQ <- round(sum(ntime)/2)
  }
  if (maxQ < L) {
    warning('maxQ must be at least L.  Setting maxQ=L.')
    maxQ <- L
  }
  # This is to avoid the area of the pesel objective function that spikes close 
  #   to rank(X), which often leads to nPC close to rank(X)
  if (maxQ > sum(ntime)*0.75) {
    warning('maxQ too high, setting to 75% of T.')
    maxQ <- round(sum(ntime)*0.75)
  }

  if(!is.logical(scale) | length(scale) != 1) stop('scale must be a logical value')

  if (multi_scans) {
    BOLD <- lapply(BOLD, scale_BOLD, scale=scale)
  } else {
    BOLD <- scale_BOLD(BOLD, scale=scale)
  }

  ### IDENTIFY AND REMOVE ANY BAD VOXELS/VERTICES
  keep <- rep(TRUE, nvox)
  keep[rowSums(is.nan(template_mean)) > 0] <- FALSE
  keep[rowSums(is.na(template_mean)) > 0] <- FALSE
  keep[rowSums(is.nan(template_var)) > 0] <- FALSE
  keep[rowSums(is.na(template_var)) > 0] <- FALSE
  keep[rowSums(is.nan(BOLD)) > 0] <- FALSE
  keep[rowSums(is.na(BOLD)) > 0] <- FALSE
  keep[rowVars(BOLD) == 0] <- FALSE
  if(sum(!keep) > 0){
    stop('flat or NA voxels detected in data or templates')
    # For this part, would need to also update "A" matrix (projection from mesh to data locations)
    # template_mean_orig <- template_mean
    # template_var_orig <- template_var
    # nvox <- sum(keep)
    # if(verbose) cat(paste0('Excluding ',sum(!keep),' bad (NA, NaN or flat) voxels/vertices from analysis.\n'))
    # template_mean <- template_mean[keep,]
    # template_var <- template_var[keep,]
    # BOLD <- BOLD[keep,]
  }

  ### 1. ESTIMATE AND DEAL WITH NUISANCE ICS (unless maxQ = L)

  if(maxQ > L){
    if (multi_scans) {
      BOLD <- lapply(BOLD, rm_nuisIC, template_mean=template_mean, Q2=Q2, Q2_max=maxQ-L, verbose=verbose)
      BOLD <- do.call(cbind, BOLD)
    } else {
      BOLD <- rm_nuisIC(BOLD, template_mean=template_mean, Q2=Q2, Q2_max=maxQ-L, verbose=verbose)
    }
  } 

  # Concatenate if multiple sessions exist.
  if (multi_scans) { BOLD <- do.call(cbind, BOLD) }
  ntime <- sum(ntime)

  #initialize mixing matrix (use dual regression-based estimate for starting value)
  dat_DR <- dual_reg(BOLD, template_mean, normA=normA)

  ### 4. RUN EM ALGORITHM!

  #Three algorithms to choose from:
  #1) FC Template ICA (new)
  #2) Template ICA
  #3) Spatial Template ICA (initialize with standard Template ICA)

  if(do_FC) {
    ### FC TEMPLATE ICA

    ## TO DO: FILL IN HERE
    prior_params <- c(0.001, 0.001) #alpha, beta (uninformative) -- note that beta is scale parameter in IG but rate parameter in the Gamma
    template_FC <- NULL
    EM_FCtemplateICA <- function(){NULL}
    resultEM <- EM_FCtemplateICA(template_mean,
                                 template_var,
                                 template_FC,
                                 prior_params, #for prior on tau^2
                                 BOLD=BOLD,
                                 AS_init = dat_DR, #initial values for A and S
                                 maxiter=maxiter,
                                 epsilon=epsilon,
                                 verbose=verbose)

    ## TO DO: FILL IN HERE

  } else {
    ### TEMPLATE ICA AND SPATIAL TEMPLATE ICA
    BOLD <- dim_reduce(BOLD, L)
    err_var <- BOLD$sigma_sq
    BOLD2 <- BOLD$data_reduced
    H <- BOLD$H
    Hinv <- BOLD$H_inv
    C_diag <- BOLD$C_diag #in original template ICA model nu^2 is separate, for spatial template ICA it is part of C
    if(do_spatial) C_diag <- C_diag * (BOLD$sigma_sq) #(nu^2)HH' in paper



    rm(BOLD)
    HA <- H %*% dat_DR$A #apply dimension reduction
    theta0 <- list(A = HA)
    # #initialize residual variance --- no longer do this, because we use dimension reduction-based estimate
    # theta0$nu0_sq = dat_list$sigma_sq
    # if(verbose) print(paste0('nu0_sq = ',round(theta0$nu0_sq,1)))
    
    #TEMPLATE ICA
    if(do_spatial) if(verbose) cat('INITIATING WITH STANDARD TEMPLATE ICA\n')
    theta00 <- theta0
    theta00$nu0_sq <- err_var
    resultEM <- EM_templateICA.independent(template_mean,
                                           template_var,
                                           BOLD=BOLD,
                                           theta0=theta00,
                                           C_diag=C_diag,
                                           maxiter=maxiter,
                                           epsilon=epsilon,
                                           verbose=verbose)
    resultEM$A <- Hinv %*% resultEM$theta_MLE$A
    class(resultEM) <- 'tICA'

    #SPATIAL TEMPLATE ICA
    if(do_spatial){

      resultEM_tICA <- resultEM
      theta0$kappa <- rep(kappa_init, L)
      if(verbose) cat('ESTIMATING SPATIAL MODEL\n')
      t000 <- Sys.time()
      resultEM <- EM_templateICA.spatial(template_mean,
                                         template_var,
                                         meshes,
                                         BOLD=BOLD,
                                         theta0,
                                         C_diag,
                                         maxiter=maxiter,
                                         epsilon=epsilon,
                                         verbose=verbose)
      #common_smoothness=common_smoothness)
      print(Sys.time() - t000)

      #organize estimates and variances in matrix form
      resultEM$subjICmean <- matrix(resultEM$subjICmean, ncol=L)
      resultEM$subjICvar <- matrix(diag(resultEM$subjICcov), ncol=L)
    }

    resultEM$A <- Hinv %*% resultEM$theta_MLE$A

    #for stICA, return tICA estimates also
    if(do_spatial){
      resultEM$result_tICA <- resultEM_tICA
      class(resultEM) <- 'stICA'
    }
  }

  #return DR estimates
  resultEM$result_DR <- dat_DR

  #This part problematic for spatial template ICA, but can bring back
  #for template ICA and FC template ICA.  When we check for bad locations,
  #can return an error only for spatial template ICA.

  #resultEM$keep <- keep
  # #map estimates & templates back to original locations
  # if(sum(!keep)>0){
  #   #estimates
  #   subjICmean <- subjICse <- matrix(nrow=length(keep), ncol=L)
  #   subjICmean[keep,] <- resultEM$subjICmean
  #   subjICse[keep,] <- resultEM$subjICse
  #   resultEM$subjICmean <- subjICmean
  #   resultEM$subjICse <- subjICse
  #   #templates
  #   resultEM$template_mean <- template_mean_orig
  #   resultEM$template_var <- template_var_orig
  # }

  return(resultEM)

}