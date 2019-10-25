#' Perform template independent component analysis (ICA) using expectation-maximization (EM)
#'
#' @param template_mean (LxV matrix) template mean estimates, i.e. mean of empirical population prior for each of L independent components
#' @param template_var (LxV matrix) template variance estimates, i.e. between-subject variance of empirical population prior for each of L ICs
#' @param BOLD (TxV matrix) BOLD fMRI data matrix, where T is the number of volumes (time points) and V is the number of brain locations
#' @param mesh Either NULL (assume spatial independence) or an object of type \code{templateICA_mesh} created by \code{make_mesh()} (spatial priors are assumed on each independent component)
#' @param maxQ Maximum number of ICs (template+nuisance) to identify (L <= maxQ <= T)
#' @param maxiter Maximum number of EM iterations
#' @param epsilon Smallest proportion change between iterations (e.g. .01 or 1%)
#'
#' @return A list containing the estimated independent components S (a LxV matrix), their mixing matrix A (a TxL matrix), and the number of nuisance ICs estimated (Q_nuis)
#' @export
#'
#' @examples
templateICA <- function(template_mean, template_var, BOLD, mesh=NULL, maxQ=NULL, maxiter=100, epsilon=0.01){

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

  #project data to the mesh locations if using spatial priors
  if(!is.null(mesh)){
    Amat <- mesh$A # n_orig x n_mesh matrix
    nmesh <- ncol(Amat)
    if(nrow(Amat) != nvox) error('Mesh projection matrix (mesh$A) must have nvox rows (nvox is the number of data locations, the columns of BOLD, template_mean and template_var)')
    template_mean_mesh <- template_mean %*% Amat
    template_var_mesh <- template_var %*% Amat
    BOLD_mesh <- BOLD %*% Amat
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
    if(!is.null(mesh)) BOLD3 <- BOLD_mesh else BOLD3 <- BOLD

  }

  ### 2. PERFORM DIMENSION REDUCTION --> BOLD4

  BOLD3 <- scale_BOLD(BOLD3)
  dat_list <- dim_reduce(BOLD3, L)
  BOLD4 <- dat_list$data_reduced
  H <- dat_list$H
  Hinv <- dat_list$H_inv
  C_diag <- dat_list$C_diag

  ### 3. SET INITIAL VALUES FOR PARAMETERS

  #initialize mixing matrix (use dual regression-based estimate for starting value)
  dat_DR <- dual_reg(BOLD3, template_mean_mesh)
  #dat_DR$A <- scale(dat_DR$A) #would this have the same effect as the code below?
  HA <- H %*% dat_DR$A #apply dimension reduction
  HA <- orthonorm(HA)  #orthogonalize
  sd_A <- sqrt(colVars(Hinv %*% HA)) #get scale of A (after reverse-prewhitening)
  HA <- HA %*% diag(1/sd_A) #standardize scale of A
  theta0 <- list(A = HA)

  #initialize residual variance
  theta0$nu0_sq = dat_list$sigma_sq #this value is huge.  what if we switch the observations and variables in the dimension reduction step?

  #initialize kappa parameters (spatial model only)
  if(!is.null(mesh)) theta0$kappa <- rep(1, L)


  ### 4. RUN EM ALGORITHM!

  if(is.null(mesh)) {
    resultEM <- EM_templateICA.independent(template_mean, template_var, BOLD4, theta0, C_diag)
  } else {
    resultEM <- EM_templateICA.spatial(template_mean_mesh, template_var_mesh, mesh, BOLD_mesh=BOLD4, theta0, C_diag)

#TO DO: revise EM_algorithm function to take mesh (made from make_mesh) instead of spde
#TO DO: make a wrapper EM_templateICA function to call the right algorithm

  }



}
