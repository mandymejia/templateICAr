#' Perform template independent component analysis (ICA) using expectation-maximization (EM)
#'
#' @param template_mean (LxV matrix) template mean estimates, i.e. mean of empirical population prior for each of L independent components
#' @param template_var (LxV matrix) template variance estimates, i.e. between-subject variance of empirical population prior for each of L ICs
#' @param BOLD (TxV matrix) BOLD fMRI data matrix, where T is the number of volumes (time points) and V is the number of brain locations
#' @param scale_BOLD Logical indicating whether BOLD data should be scaled by the spatial standard deviation before model fitting. If done when estimating templates, should be done here too.
#' @param mesh Either NULL (assume spatial independence) or an object of type \code{templateICA_mesh} created by \code{make_mesh()} (spatial priors are assumed on each independent component)
#' @param maxQ Maximum number of ICs (template+nuisance) to identify (L <= maxQ <= T)
#' @param common_smoothness If TRUE, use the common smoothness version of the spatial template ICA model, which assumes that all IC's have the same smoothness parameter, \eqn{\kappa}
#' @param maxiter Maximum number of EM iterations
#' @param epsilon Smallest proportion change between iterations (e.g. .01)
#' @param verbose If TRUE, display progress of algorithm
#' @param return_kappa_fun If TRUE, return the log likelihood as a function of kappa in a neighborhood of the MLE (slow) (common smoothness model only)
#' @param kappa_init (Optional) Starting value for kappa.  If NULL, starting value will be determined automatically.
#' @param dim_reduce_flag If FALSE, do not perform dimension reduction in spatial template ICA. Not applicable for standard template ICA.
#'
#' @return A list containing the estimated independent components S (a LxV matrix), their mixing matrix A (a TxL matrix), and the number of nuisance ICs estimated (Q_nuis)
#' @export
#' @importFrom INLA inla inla.spde.result
#'
templateICA <- function(template_mean, template_var, BOLD, scale_BOLD=FALSE, mesh=NULL, maxQ=NULL, common_smoothness=TRUE, maxiter=100, epsilon=0.01, verbose=TRUE, return_kappa_fun=FALSE, kappa_init=NULL, dim_reduce_flag=FALSE){

  ntime <- nrow(BOLD) #length of timeseries
  nvox <- ncol(BOLD) #number of data locations
  L <- nrow(template_mean) #number of ICs

  #check that the number of data locations (nvox), time points (ntime) and ICs (L) makes sense
  if(ntime > nvox) warning('More time points than voxels. Are you sure?')
  if(L > nvox) stop('The arguments you supplied suggest that you want to estimate more ICs than you have data locations.  Please check the orientation and size of template_mean, template_var and BOLD.')
  if(L > ntime) stop('The arguments you supplied suggest that you want to estimate more ICs than you have time points.  Please check the orientation and size of template_mean, template_var and BOLD.')

  #check that all arguments have consistent number of data locations (nvox) and ICs (L)
  if(ncol(template_mean) != nvox | ncol(template_var) != nvox) stop('template_mean, template_var and BOLD must have same number of data locations (columns), but they do not.')
  if(nrow(template_var) != L) stop('template_mean and template_var must have the same number of ICs (rows), but they do not.')

  #check that the supplied mesh object is of type templateICA_mesh
  if(is.null(mesh)){
    message('No mesh supplied: Using standard template ICA model, which assumes spatial independence. If this is not what you want, stop and supply a valid mesh. See help(make_mesh).')
  } else if(class(mesh) != 'templateICA_mesh'){
    stop('mesh argument should be of class templateICA_mesh. See help(make_mesh).')
  }

  if(!is.null(maxQ)){
    if(round(maxQ) != maxQ | maxQ <= 0) stop('maxQ must be NULL or a round positive number')
  }

  if(round(maxiter) != maxiter | maxiter <= 0) stop('maxiter must be a positive integer')

  if(!is.null(kappa_init)){
    if(length(kappa_init) != 1 | kappa_init <= 0) stop('kappa_init must be a positive scalar or NULL')
  }

  #check that maxQ makes sense
  if(is.null(maxQ)) maxQ <- ntime
  if(maxQ < L){
    warning('maxQ must be at least L.  Setting maxQ=L.')
    maxQ <- L
  }
  if(maxQ > ntime){
    warning('maxQ must be no more than T.  Setting maxQ = T.')
    maxQ <- ntime
  }

  if(class(scale_BOLD) != 'logical' | length(scale_BOLD) != 1) stop('scale_BOLD must be a logical value')
  if(class(common_smoothness) != 'logical' | length(common_smoothness) != 1) stop('common_smoothness must be a logical value')
  if(class(return_kappa_fun) != 'logical' | length(return_kappa_fun) != 1) stop('return_kappa_fun must be a logical value')
  if(return_kappa_fun==TRUE & common_smoothness==FALSE){
    warning('I can only return the kappa function when common_smoothness==TRUE.  Setting return_kappa_fun to FALSE.')
    return_kappa_fun <- FALSE
  }


  ### 1. ESTIMATE AND DEAL WITH NUISANCE ICS (unless maxQ = L)

  if(maxQ > L){

    ## use BOLD_mesh! if(!is.null(mesh))

    #i. PERFORM DUAL REGRESSION TO GET INITIAL ESTIMATE OF TEMPLATE ICS
    #ii. SUBTRACT THOSE ESTIMATES FROM THE ORIGINAL DATA --> BOLD2
    #iii. ESTIMATE THE NUMBER OF REMAINING ICS USING THE MINKA METHOD (TO DO)
    #iv. ESTIMATE THE NUISANCE ICS USING GIFT/INFOMAX (TO DO)
    #v. SUBTRACT THOSE ESTIMATES FROM THE ORIGINAL DATA --> BOLD3

  } else {

    # USE ORIGINAL DATA, SINCE WE ARE ASSUMING NO NUISANCE COMPONENTS
    BOLD3 <- BOLD

  }

  ### 2. PERFORM DIMENSION REDUCTION --> BOLD4

  if(dim_reduce_flag) if(verbose) print('PERFORMING DIMENSION REDUCTION')

  BOLD3 <- scale_BOLD(BOLD3, scale=FALSE)
  dat_list <- dim_reduce(BOLD3, L)

  BOLD4 <- dat_list$data_reduced
  H <- dat_list$H
  Hinv <- dat_list$H_inv
  C_diag <- dat_list$C_diag

  ### 3. SET INITIAL VALUES FOR PARAMETERS

  #initialize mixing matrix (use dual regression-based estimate for starting value)
  dat_DR <- dual_reg(BOLD3, template_mean)
  #dat_DR$A <- scale(dat_DR$A) #would this have the same effect as the code below?

  # Keep?
  HA <- H %*% dat_DR$A #apply dimension reduction
  # HA <- orthonorm(HA)  #orthogonalize
  # sd_A <- sqrt(colVars(Hinv %*% HA)) #get scale of A (after reverse-prewhitening)
  # HA <- HA %*% diag(1/sd_A) #standardize scale of A
  theta0 <- list(A = HA)


  #initialize residual variance
  theta0$nu0_sq = dat_list$sigma_sq
  if(verbose) print(paste0('nu0_sq = ',round(theta0$nu0_sq,1)))


  ### 4. RUN EM ALGORITHM!

  template_mean <- t(scale(t(template_mean), scale=FALSE)) #center template IC means

  if(!is.null(mesh)) print('INITIATING WITH STANDARD TEMPLATE ICA')
  resultEM <- EM_templateICA.independent(template_mean, template_var, BOLD4, theta0, C_diag, maxiter=maxiter)

  if(!is.null(mesh)){

    if(dim_reduce_flag == FALSE){
      BOLD4 <- BOLD3
      C_diag <- rep(1, ntime)
      theta0$A <- dat_DR$A
    }

    # theta0$A <- resultEM$theta_MLE$A
    # theta0$nu0_sq <- resultEM$theta_MLE$nu0_sq
    # if(verbose) print(paste0('Starting value for nu0_sq = ',round(theta0$nu0_sq,1)))

    #starting value for kappas
    if(is.null(kappa_init)){
      if(verbose) print('Running INLA on tICA estimates to determine starting value(s) for kappa')
      theta0$kappa <- rep(NA, L)
      for(q in 1:L){
        print(paste0('IC ',q,' of ',L))
        d_q <- resultEM$subjICmean[q,] - template_mean[q,]
        data_inla_q <- list(y = d_q, x = mesh$mesh$idx$loc)
        formula_q <- y ~ -1 + f(x, model = mesh$spde)
        result_q <- inla(formula_q, data = data_inla_q, verbose = FALSE)
        result_spde_q <- inla.spde.result(result_q, name='x', spde=mesh$spde)
        theta0$kappa[q] <- exp(result_spde_q$summary.log.kappa$mean)
      }
      if(common_smoothness) theta0$kappa <- rep(mean(theta0$kappa), L)
      if(verbose) print(paste0('Starting value for kappa = ',paste(round(theta0$kappa,3), collapse=', ')))
    } else {
      theta0$kappa <- rep(kappa_init, L)
    }

    #project BOLD and templates to mesh locations
    Amat <- mesh$A # n_orig x n_mesh matrix
    nmesh <- ncol(Amat)
    if(nrow(Amat) != nvox) stop('Mesh projection matrix (mesh$A) must have nvox rows (nvox is the number of data locations, the columns of BOLD, template_mean and template_var)')
    template_mean <- template_mean %*% Amat
    template_var <- template_var %*% Amat

    BOLD4_orig <- BOLD4
    BOLD4 <- BOLD4 %*% Amat
    #BOLD4 <- BOLD3 %*% Amat  #Keep (no dim reduction)?

    print('RUNNING SPATIAL TEMPLATE ICA')
    resultEM <- EM_templateICA.spatial(template_mean, template_var, mesh, BOLD=BOLD4, theta0, C_diag, Hinv, common_smoothness=common_smoothness, maxiter=maxiter, return_kappa_fun=return_kappa_fun, verbose=verbose, dim_reduce_flag=dim_reduce_flag)

    #project estimates back to data locations
    resultEM$subjICmean2 <- t(matrix(resultEM$subjICmean, ncol=L)) %*% t(Amat)

  }

  if(dim_reduce_flag==TRUE) {
    A <- Hinv %*% resultEM$theta_MLE$A
    resultEM$A <- A
  }

  #regression-based final estimate of A
  if(is.null(mesh)){ S <- resultEM$subjICmean } else { S <- resultEM$subjICmean2 }
  tmp <- dual_reg(BOLD3, S)
  resultEM$A_reg <- tmp$A
  return(resultEM)
}
