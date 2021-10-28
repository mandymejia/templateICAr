#' Template ICA for CIFTI
#'
#' Run template ICA based on CIFTI-format BOLD data and CIFTI-based template
#'
#' @param cifti_fname File path of CIFTI-format timeseries data (ending in .dtseries.nii).
#' @param template Result of call to \code{estimate_template.cifti}, either (1) an object of class
#' \code{template_cifti} OR (2) the file path and basename prefix (the part before "_mean.dscalar.nii"
#' and "_var.dscalar.nii") of cifti files written out by \code{estimate_template.cifti}. For
#' spatial template ICA (\code{spatial_model=TRUE}), must be the latter.
#' @param spatial_model Should spatial modeling be performed? (Default \eqn{FALSE}) If \code{TRUE}. surface
#' geometries will be used to fit a spatial Bayesian model. Computational demands (time and memory)
#' are much higher with spatial_model=TRUE than spatial_model=FALSE, but results are more accurate.
#' @param brainstructures Character vector indicating which brain structure(s)
#'  to obtain: \code{"left"} (left cortical surface), \code{"right"} (right
#'  cortical surface) and/or \code{"subcortical"} (subcortical and cerebellar
#'  gray matter). Can also be \code{"all"} (obtain all three brain structures).
#'  Default: \code{c("left","right")} (cortical surface only).  For spatial
#'  modeling (performed if \code{surfL_fname} and/or \code{surfR_fname} provided), subcortical
#'  locations will be ignored.
#' @param surfL_fname (Required for spatial modeling) File path of GIFTI surface geometry
#'  file representing the left cortex. Only required if \code{brainstructures} includes \code{"left"}.
#' @param surfR_fname (Required for spatial modeling) File path of GIFTI surface geometry
#'  file representing the right cortex. Only required if \code{brainstructures} includes \code{"right"}.
#' @param resamp_res (Only recommended for spatial modeling) Target resolution for resampling (number of
#'  cortical surface vertices per hemisphere). A value less than 10000 is recommended for computational
#'  feasibility. If \code{NULL} (default) or \code{FALSE}, do not perform resampling.
#' @param scale Logical indicating whether BOLD data should be scaled by the spatial
#' standard deviation before model fitting. If done when estimating templates, should be done here too.
#' @param Q2 The number of nuisance ICs to identify. If NULL, will be estimated. Only provide \eqn{Q2} or \eqn{maxQ} but not both.
#' @param maxQ Maximum number of ICs (template+nuisance) to identify (L <= maxQ <= T). Only provide \eqn{Q2} or \eqn{maxQ} but not both.
#' @param maxiter Maximum number of EM iterations. Default: 100.
#' @param epsilon Smallest proportion change between iterations. Default: 0.001.
#' @param verbose If \code{TRUE} (default), display progress of algorithm.
#' @param kappa_init Starting value for kappa.  Default: \code{0.2}.
# @param common_smoothness If \code{TRUE}. use the common smoothness version of the spatial template ICA model, which assumes that all IC's have the same smoothness parameter, \eqn{\kappa}
#' @param write_dir Where should any output files be written? \code{NULL} (default) will write them to the current working directory.
#' @param rm_mwall Should medial wall (missing data) locations be removed from mesh?  If TRUE, faster estimation but less accurate estimates on boundary of wall.
#' @param time_inds Subset of fMRI BOLD volumes to include in analysis. If NULL (default), use all volumes.
#'
# @importFrom INLA inla.pardiso.check inla.setOption
#' @importFrom ciftiTools read_cifti newdata_xifti
#'
#' @return A list containing the subject IC estimates (class 'xifti'), the subject IC variance estimates (class 'xifti'), and the result of the model call to \code{templateICA} (class 'dICA')
#'
#' @export
#'
templateICA.cifti <- function(cifti_fname,
                              template,
                              brainstructures=c("left","right"),
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
                              #common_smoothness=TRUE,
                              kappa_init=0.2,
                              write_dir=NULL,
                              rm_mwall=TRUE,
                              time_inds=NULL){

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
    INLA_check()
    if(do_sub) stop('If spatial_model=TRUE, only applicable to "left" and/or "right" brainstructures. Check brainstructures argument and try again.')
    if(!is.character(template)) stop('If spatial_model=TRUE, template argument must be file path prefix to cifti files written by estimate_template.cifti().')
    flag <- INLA::inla.pardiso.check()
    if(grepl('FAILURE',flag)) stop('PARDISO IS NOT INSTALLED OR NOT WORKING. PARDISO for R-INLA is required for computational efficiency. If you already have a PARDISO / R-INLA License, run inla.setOption(pardiso.license = "/path/to/license") and try again.  If not, run inla.pardiso() to obtain a license.')
    INLA::inla.setOption(smtp='pardiso')
  }

  if(!spatial_model){
    if(!is.null(resamp_res)) warning('Resampling typically only required for spatial modeling.  Recommend setting resamp_res to NULL.')
  }

  # READ IN BOLD TIMESERIES DATA
  if (!file.exists(cifti_fname)) stop(paste0('The BOLD timeseries file ',cifti_fname,' does not exist.'))
  if (verbose) cat('Reading in BOLD timeseries data.\n')
  BOLD <- read_cifti(
    cifti_fname, surfL_fname = surfL_fname, surfR_fname = surfR_fname,
    brainstructures = brainstructures, resamp_res=resamp_res
  )
  # Separate data and metadata
  BOLD_meta <- BOLD; BOLD_meta["data"] <- list(NULL)
  BOLD <- as.matrix(BOLD)

  # GET TEMPLATE MEAN AND VARIANCE (xifti objects)
  template_class <- class(template)
  if (template_class != 'template_cifti') {
    if (verbose) cat('Reading in templates.\n')
    fname_mean <- paste0(template,'_mean.dscalar.nii')
    fname_var <- paste0(template,'_var.dscalar.nii')
    if(!file.exists(fname_mean)) stop(paste0('The file ', fname_mean, ' does not exist.'))
    if(!file.exists(fname_var)) stop(paste0('The file ', fname_var, ' does not exist.'))
    template <- list(mean=NULL, var=NULL)
    template$mean <- read_cifti(
      fname_mean, brainstructures=brainstructures, resamp_res=resamp_res, 
      surfL_fname = surfL_fname, surfR_fname = surfR_fname
    )
    template$var <- read_cifti(
      fname_var, brainstructures=brainstructures, resamp_res=resamp_res, 
      surfL_fname = surfL_fname, surfR_fname = surfR_fname
    )
  } else {
    stop(
      'template argument must be an object of class template_cifti or file path',
      'prefix to result of estimate_template.cifti() (same as out_fname argument',
      'passed to estimate_template.cifti().'
    )
  }

  #models_list <- vector('list', length=length(models))
  #names(models_list) <- models
  #for(mod in models){

  # IF SPATIAL MODELING, CONSTRUCT MESH
  if(spatial_model){
    meshes <- NULL
    if(do_left) {
      surf <- BOLD_meta$surf$cortex_left
      wall_mask <- which(BOLD_meta$meta$cortex$medial_wall_mask$left)
      if (rm_mwall) mesh <- make_mesh(surf = surf, inds_mesh = wall_mask) #remove wall for greater computational efficiency
      if(!rm_mwall) mesh <- make_mesh(surf = surf, inds_data = wall_mask) #retain wall in mesh for more accuracy along boundary with wall
      meshes <- c(meshes, list(mesh))
    }
    if(do_right) {
      surf <- BOLD_meta$surf$cortex_right
      wall_mask <- which(BOLD_meta$meta$cortex$medial_wall_mask$right)
      if(rm_mwall) mesh <- make_mesh(surf = surf, inds_mesh = wall_mask) #remove wall for greater computational efficiency
      if(!rm_mwall) mesh <- make_mesh(surf = surf, inds_data = wall_mask) #retain wall in mesh for more accuracy along boundary with wall
      meshes <- c(meshes, list(mesh))
    }
  } else {
    meshes <- NULL
  }

  if(!is.null(time_inds)){
    all_inds <- seq(ncol(BOLD))
    if (any(!(time_inds %in% all_inds))) {
      stop(paste0('time_inds contains indices outside of the range of 1 to ', max(all_inds)))
    }
    BOLD <- BOLD[,time_inds,drop=FALSE]
  }

  # CALL TEMPLATE ICA FUNCTION
  result <- templateICA(
    template_mean = as.matrix(template$mean), 
    template_var = as.matrix(template$var),
    BOLD = BOLD,
    scale = scale,
    meshes = meshes,
    Q2 = Q2,
    maxQ=maxQ,
    maxiter=maxiter,
    epsilon=epsilon,
    verbose=verbose,
    #common_smoothness=common_smoothness,
    kappa_init=kappa_init
  )
  rm(BOLD)

  result$subjICmean <- newdata_xifti(template$mean, result$subjICmean)
  result$subjICse <- newdata_xifti(template$var, result$subjICse)

  if (spatial_model) {
    result$result_tICA$subjICmean <- newdata_xifti(template$mean, result$result_tICA$subjICmean)
    result$result_tICA$subjICse <- newdata_xifti(template$var, result$result_tICA$subjICse)
  }

  return(result)
}