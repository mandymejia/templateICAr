#' Run template ICA based on CIFTI-format BOLD data and CIFTI-based template
#'
#' @param cifti_fname File path of CIFTI-format timeseries data (ending in .dtseries.nii).
#' @param template Result of call to \code{estimate_template.cifti()}, either (1) an object of class
#' \code{template.cifti} OR (2) the file path and basename prefix (the part before "_mean.dscalar.nii"
#' and "_var.dscalar.nii") of cifti files written out by \code{estimate_template.cifti()}. For
#' spatial template ICA (\code{spatial_model=TRUE}), must be the latter.
#' @param spatial_model Should spatial modeling be performed? (Default FALSE) If TRUE, surface
#' geometries will be used to fit a spatial Bayesian model. Computational demands (time and memory)
#' are much higher with spatial_model=TRUE than spatial_model=FALSE, but results are more accurate.
#' @param brainstructures Character vector indicating which brain structure(s)
#'  to obtain: \code{"left"} (left cortical surface), \code{"right"} (right
#'  cortical surface) and/or \code{"subcortical"} (subcortical and cerebellar
#'  gray matter). Can also be \code{"all"} (obtain all three brain structures).
#'  Default: \code{c("left","right")} (cortical surface only).  For spatial
#'  modeling (performed if `surfL_fname` and/or `surfR_fname` provided), subcortical
#'  locations will be ignored.
#' @param surfL_fname (Required for spatial modeling) File path of GIFTI surface geometry
#'  file representing the left cortex. Only required if \code{brainstructures} includes \code{"right"}.
#' @param surfR_fname (Required for spatial modeling) File path of GIFTI surface geometry
#'  file representing the right cortex. Only required if \code{brainstructures} includes \code{"left"}.
#' @param resamp_res (Only recommended for spatial modeling) Target resolution for resampling (number of
#'  cortical surface vertices per hemisphere). A value less than 10000 is recommended for computational
#'  feasibility. If \code{NULL} (default) or \code{FALSE}, do not perform resampling.
#' @param scale Logical indicating whether BOLD data should be scaled by the spatial
#' standard deviation before model fitting. If done when estimating templates, should be done here too.
#' @param Q2 The number of nuisance ICs to identify. If NULL, will be estimated. Only provide Q2 or maxQ but not both.
#' @param maxQ Maximum number of ICs (template+nuisance) to identify (L <= maxQ <= T). Only provide Q2 or maxQ but not both.
#' @param maxiter Maximum number of EM iterations
#' @param epsilon Smallest proportion change between iterations (e.g. .01)
#' @param verbose If TRUE, display progress of algorithm
#' @param kappa_init Starting value for kappa.  If NULL, starting value will be determined automatically.
#' @param common_smoothness If TRUE, use the common smoothness version of the spatial template ICA model, which assumes that all IC's have the same smoothness parameter, \eqn{\kappa}
#' @param write_dir Where should any output files be written? \code{NULL} (default) will write them to the current working directory.
#'
#' @importFrom INLA inla.pardiso.check inla.setOption
#' @importFrom ciftiTools read_cifti
#'
#' @return A list containing the subject IC estimates (class 'xifti'), the subject IC variance estimates (class 'xifti'), and the result of the model call to \code{templateICA} (class 'dICA')
#' @export
#'
templateICA.cifti <- function(cifti_fname,
                              template,
                              brainstructures=c('left','right'),
                              spatial_model=FALSE,
                              surfL_fname=NULL,
                              surfR_fname=NULL,
                              resamp_res=NULL,
                              scale=TRUE,
                              Q2=NULL,
                              maxQ=NULL,
                              maxiter=100,
                              epsilon=0.001,
                              verbose=TRUE,
                              common_smoothness=TRUE,
                              kappa_init=NULL,
                              write_dir=NULL){

  if (is.null(write_dir)) { write_dir <- getwd() }


  brainstructures <- match_input(
    brainstructures, c("left","right","subcortical","all"),
    user_value_label="brainstructures"
  )
  if ("all" %in% brainstructures) {
    brainstructures <- c("left","right","subcortical")
  }

  do_left <- do_right <- do_sub <- FALSE
  if('left' %in% brainstructures) do_left <- TRUE
  if('right' %in% brainstructures) do_right <- TRUE
  if('subcortical' %in% brainstructures) do_sub <- TRUE

  if(spatial_model){
    if(do_sub) stop('If spatial_model=TRUE, only applicable to "left" and/or "right" brainstructures. Check brainstructures argument and try again.')
    if(!is.character(template)) stop('If spatial_model=TRUE, template argument must be file path prefix to cifti files written by estimate_template.cifti().')
    flag <- inla.pardiso.check()
    if(grepl('FAILURE',flag)) stop('PARDISO IS NOT INSTALLED OR NOT WORKING. PARDISO for R-INLA is required for computational efficiency. If you already have a PARDISO / R-INLA License, run inla.setOption(pardiso.license = "/path/to/license") and try again.  If not, run inla.pardiso() to obtain a license.')
    inla.setOption(smtp='pardiso')
  }

  if(!spatial_model){
    if(!is.null(resamp_res)) warning('Resampling typically only required for spatial modeling.  Recommend setting resamp_res to NULL.')
  }

  # GET TEMPLATE MEAN AND VARIANCE (xifti objects)
  template_class <- class(template)
  if(template_class=='template.cifti'){
    template_mean <- template$template_mean
    template_var <- template$template_var
  } else if(is.character(template)){
    fname_mean <- paste0(template,'_mean.dscalar.nii')
    fname_var <- paste0(template,'_var.dscalar.nii')
    if(!file.exists(fname_mean)) stop(paste0('The file ', fname_mean, ' does not exist.'))
    if(!file.exists(fname_var)) stop(paste0('The file ', fname_var, ' does not exist.'))
    template_mean <- read_cifti(fname_mean, brainstructures=brainstructures, resamp_res=resamp_res, surfL_fname = surfL_fname, surfR_fname = surfR_fname)
    template_var <- read_cifti(fname_var, brainstructures=brainstructures, resamp_res=resamp_res, surfL_fname = surfL_fname, surfR_fname = surfR_fname)
  } else {
    stop('template argument must be an object of class template.cifti or file path prefix to result of estimate_template.cifti() (same as out_fname argument passed to estimate_template.cifti().')
  }


  # READ IN BOLD TIMESERIES DATA
  if(!file.exists(cifti_fname)) stop(paste0('The BOLD timeseries file ',cifti_fname,' does not exist.'))
  BOLD_cifti <- read_cifti(cifti_fname,
                     surfL_fname = surfL_fname,
                     surfR_fname = surfR_fname,
                     brainstructures = brainstructures,
                     resamp_res=resamp_res)

  #TO DO: SINGLE SPATIAL MODEL

  # IF SPATIAL MODELING, LOOP OVER HEMISPHERES
  if(spatial_model) {
    if(do_left & !do_right) models <- c('lh')
    if(do_right & !do_left) models <- c('rh')
    if(do_left & do_right) models <- c('lh','rh')
  } else {
    models <- 'single'
  }

  #set up xifti objects for IC mean and variance estimates
  clear_data <- function(x){
    if(!is.null(x$data$cortex_left)) x$data$cortex_left <- matrix(0, nrow(x$data$cortex_left), 1)
    if(!is.null(x$data$cortex_right)) x$data$cortex_right <- matrix(0, nrow(x$data$cortex_right), 1)
    if(!is.null(x$data$subcort)) x$data$subcort <- matrix(0, nrow(x$data$subcort), 1)
    return(x)
  }
  subjICmean_xifti <- subjICvar_xifti <- clear_data(BOLD_cifti)

  models_list <- vector('list', length=length(models))
  names(models_list) <- models
  for(mod in models){

    # IF SPATIAL MODELING, CONSTRUCT MESH
    if(spatial_model){
      if(mod=='lh') {
        surf <- BOLD_cifti$surf$cortex_left
        wall_mask <- which(BOLD_cifti$meta$cortex$medial_wall_mask$left)
      }
      if(mod=='rh') {
        surf <- BOLD_cifti$surf$cortex_right
        wall_mask <- which(BOLD_cifti$meta$cortex$medial_wall_mask$right)
      }
      mesh <- make_mesh(surf = surf, inds_mesh = wall_mask) #remove wall for greater computational efficiency
    } else {
      mesh <- NULL
    }

    #TO DO: Test spatial model through this function.  Are the medial wall locations dealt with appropriately?

    # FORM DATA MATRIX AND TEMPLATE MATRICES
    if(mod=='lh') {
      BOLD_mat <- BOLD_cifti$data$cortex_left
      template_mean_mat <- template_mean$data$cortex_left
      template_var_mat <- template_var$data$cortex_left
    } else if(mod=='rh') {
      BOLD_mat <- BOLD_cifti$data$cortex_right
      template_mean_mat <- template_mean$data$cortex_right
      template_var_mat <- template_var$data$cortex_right
    } else {
      BOLD_mat <- rbind(BOLD_cifti$data$cortex_left, BOLD_cifti$data$cortex_right, BOLD_cifti$data$subcort)
      template_mean_mat <- rbind(template_mean$data$cortex_left, template_mean$data$cortex_right, template_mean$data$subcort)
      template_var_mat <- rbind(template_var$data$cortex_left, template_var$data$cortex_right, template_var$data$subcort)
    }

    # CALL TEMPLATE ICA FUNCTION

    result_mod <- templateICA(template_mean = template_mean_mat,
                              template_var = template_var_mat,
                              BOLD = BOLD_mat,
                              scale = scale,
                              mesh = mesh,
                              Q2 = Q2,
                              maxQ=maxQ,
                              maxiter=maxiter,
                              epsilon=epsilon,
                              verbose=verbose,
                              common_smoothness=common_smoothness,
                              kappa_init=kappa_init)

    models_list[[which(models==mod)]] <- result_mod

    # MAP ESTIMATES AND VARIANCE TO XIFTI FORMAT
    if(mod=='lh') { #left hemisphere spatial model
      subjICmean_xifti$data$cortex_left <- result_mod$subjICmean
      subjICvar_xifti$data$cortex_left <- result_mod$subjICvar
    } else if(mod=='rh') { #right hemisphere spatial model
      subjICmean_xifti$data$cortex_right <- result_mod$subjICmean
      subjICvar_xifti$data$cortex_left <- result_mod$subjICvar
    } else { #single non-spatial model
      n_left <- n_right <- n_sub <- 0
      if(do_left) {
        n_left <- nrow(subjICmean_xifti$data$cortex_left)
        subjICmean_xifti$data$cortex_left <- result_mod$subjICmean[1:n_left,]
        subjICvar_xifti$data$cortex_left <- result_mod$subjICvar[1:n_left,]
      }
      if(do_right) {
        n_right <- nrow(subjICmean_xifti$data$cortex_right)
        subjICmean_xifti$data$cortex_right <- result_mod$subjICmean[n_left+(1:n_right),]
        subjICvar_xifti$data$cortex_right <- result_mod$subjICvar[n_left+(1:n_right),]
      }
      if(do_sub) {
        n_sub <- nrow(subjICmean_xifti$data$subcort)
        subjICmean_xifti$data$subcort <- result_mod$subjICmean[n_left + n_right + (1:n_sub),]
        subjICvar_xifti$data$subcort <- result_mod$subjICvar[n_left + n_right + (1:n_sub),]
      }
    }
  }

    list(subjICmean_xifti = subjICmean_xifti,
         subjICvar_xifti = subjICvar_xifti,
         model_result = models_list)

}

