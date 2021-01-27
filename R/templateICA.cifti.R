#' Template ICA for CIFTI
#'
#' Run template ICA based on CIFTI-format BOLD data and CIFTI-based template
#'
#' @param cifti_fname File path of CIFTI-format timeseries data (ending in .dtseries.nii).
#' @param template Result of call to \code{estimate_template.cifti}, either (1) an object of class
#' \code{template.cifti} OR (2) the file path and basename prefix (the part before "_mean.dscalar.nii"
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
#' @param surfL_fname_vis (Optional) File path of GIFTI surface geometry file representing the left cortex
#' to be used for visualization. Only provide if \code{brainstructures} includes \code{"left"}.
#' @param surfR_fname_vis (Optional) File path of GIFTI surface geometry file representing the right cortex
#' to be used for visualization. Only provide if \code{brainstructures} includes \code{"right"}.
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
#'
# @importFrom INLA inla.pardiso.check inla.setOption
#' @importFrom ciftiTools read_cifti
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
                              surfL_fname_vis=NULL,
                              surfR_fname_vis=NULL,
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
                              rm_mwall=TRUE){

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
    if (!requireNamespace("INLA", quietly = TRUE)) {
      stop(
        paste0(
          "Package \"INLA\" needed to for spatial modeling.",
          "Please install it at http://www.r-inla.org/download.",
        ), call. = FALSE
      )
    }
    if(do_sub) stop('If spatial_model=TRUE, only applicable to "left" and/or "right" brainstructures. Check brainstructures argument and try again.')
    if(!is.character(template)) stop('If spatial_model=TRUE, template argument must be file path prefix to cifti files written by estimate_template.cifti().')
    flag <- INLA::inla.pardiso.check()
    if(grepl('FAILURE',flag)) stop('PARDISO IS NOT INSTALLED OR NOT WORKING. PARDISO for R-INLA is required for computational efficiency. If you already have a PARDISO / R-INLA License, run inla.setOption(pardiso.license = "/path/to/license") and try again.  If not, run inla.pardiso() to obtain a license.')
    INLA::inla.setOption(smtp='pardiso')
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
    if(verbose) cat('Reading in templates.\n')
    fname_mean <- paste0(template,'_mean.dscalar.nii')
    fname_var <- paste0(template,'_var.dscalar.nii')
    if(!file.exists(fname_mean)) stop(paste0('The file ', fname_mean, ' does not exist.'))
    if(!file.exists(fname_var)) stop(paste0('The file ', fname_var, ' does not exist.'))
    if(!is.null(surfL_fname_vis)) surfL_fname_template <- surfL_fname_vis else surfL_fname_template <- surfL_fname
    if(!is.null(surfR_fname_vis)) surfR_fname_template <- surfR_fname_vis else surfR_fname_template <- surfR_fname
    template_mean <- read_cifti(fname_mean, brainstructures=brainstructures, resamp_res=resamp_res, surfL_fname = surfL_fname_template, surfR_fname = surfR_fname_template)
    template_var <- read_cifti(fname_var, brainstructures=brainstructures, resamp_res=resamp_res, surfL_fname = surfL_fname_template, surfR_fname = surfR_fname_template)
  } else {
    stop('template argument must be an object of class template.cifti or file path prefix to result of estimate_template.cifti() (same as out_fname argument passed to estimate_template.cifti().')
  }


  # READ IN BOLD TIMESERIES DATA
  if(!file.exists(cifti_fname)) stop(paste0('The BOLD timeseries file ',cifti_fname,' does not exist.'))
  if(verbose) cat('Reading in BOLD timeseries data.\n')
  BOLD_cifti <- read_cifti(cifti_fname,
                     surfL_fname = surfL_fname,
                     surfR_fname = surfR_fname,
                     brainstructures = brainstructures,
                     resamp_res=resamp_res)

  #set up xifti objects for IC mean and variance estimates
  subjICmean_xifti <- subjICvar_xifti <- clear_data(template_mean)

  #models_list <- vector('list', length=length(models))
  #names(models_list) <- models
  #for(mod in models){

  # IF SPATIAL MODELING, CONSTRUCT MESH
  if(spatial_model){
    meshes <- NULL
    if(do_left) {
      surf <- BOLD_cifti$surf$cortex_left
      wall_mask <- which(BOLD_cifti$meta$cortex$medial_wall_mask$left)
      if(rm_mwall) mesh <- make_mesh(surf = surf, inds_mesh = wall_mask) #remove wall for greater computational efficiency
      if(!rm_mwall) mesh <- make_mesh(surf = surf, inds_data = wall_mask) #retain wall in mesh for more accuracy along boundary with wall
      meshes <- c(meshes, list(mesh))
    }
    if(do_right) {
      surf <- BOLD_cifti$surf$cortex_right
      wall_mask <- which(BOLD_cifti$meta$cortex$medial_wall_mask$right)
      if(rm_mwall) mesh <- make_mesh(surf = surf, inds_mesh = wall_mask) #remove wall for greater computational efficiency
      if(!rm_mwall) mesh <- make_mesh(surf = surf, inds_data = wall_mask) #retain wall in mesh for more accuracy along boundary with wall
      meshes <- c(meshes, list(mesh))
    }
  } else {
    meshes <- NULL
  }

  # FORM DATA MATRIX AND TEMPLATE MATRICES
  BOLD_mat <- do.call(rbind, BOLD_cifti$data)
  template_mean_mat <- do.call(rbind, template_mean$data)
  template_var_mat <- do.call(rbind, template_var$data)

  # if(do_left) {
  #     BOLD_mat <- rbind(BOLD_mat, BOLD_cifti$data$cortex_left)
  #     template_mean_mat <- rbind(template_mean_mat, template_mean$data$cortex_left)
  #     template_var_mat <- rbind(template_var_mat, template_var$data$cortex_left)
  #   }
  # if(mod=='rh') {
  #     BOLD_mat <- BOLD_cifti$data$cortex_right
  #     template_mean_mat <- template_mean$data$cortex_right
  #     template_var_mat <- template_var$data$cortex_right
  #   } else {
  #     BOLD_mat <- rbind(BOLD_cifti$data$cortex_left, BOLD_cifti$data$cortex_right, BOLD_cifti$data$subcort)
  #     template_mean_mat <- rbind(template_mean$data$cortex_left, template_mean$data$cortex_right, template_mean$data$subcort)
  #     template_var_mat <- rbind(template_var$data$cortex_left, template_var$data$cortex_right, template_var$data$subcort)
  #   }

  # CALL TEMPLATE ICA FUNCTION

  result <- templateICA(template_mean = template_mean_mat,
                            template_var = template_var_mat,
                            BOLD = BOLD_mat,
                            scale = scale,
                            meshes = meshes,
                            Q2 = Q2,
                            maxQ=maxQ,
                            maxiter=maxiter,
                            epsilon=epsilon,
                            verbose=verbose,
                            #common_smoothness=common_smoothness,
                            kappa_init=kappa_init)


  # MAP ESTIMATES AND VARIANCE TO XIFTI FORMAT
  n_left <- n_right <- n_sub <- 0
  if(do_left) {
    n_left <- nrow(subjICmean_xifti$data$cortex_left)
    subjICmean_xifti$data$cortex_left <- result$subjICmean[1:n_left,]
    subjICvar_xifti$data$cortex_left <- result$subjICvar[1:n_left,]
  }
  if(do_right) {
    n_right <- nrow(subjICmean_xifti$data$cortex_right)
    subjICmean_xifti$data$cortex_right <- result$subjICmean[n_left+(1:n_right),]
    subjICvar_xifti$data$cortex_right <- result$subjICvar[n_left+(1:n_right),]
  }
  if(do_sub) {
    n_sub <- nrow(subjICmean_xifti$data$subcort)
    subjICmean_xifti$data$subcort <- result$subjICmean[n_left + n_right + (1:n_sub),]
    subjICvar_xifti$data$subcort <- result$subjICvar[n_left + n_right + (1:n_sub),]
  }

  subjICmean_xifti$meta$cifti$names <- paste0('IC ', seq_len(length(subjICmean_xifti$meta$cifti$names)))
  subjICvar_xifti$meta$cifti$names <- paste0('IC ', seq_len(length(subjICmean_xifti$meta$cifti$names)))

  RESULT <- list(
    subjICmean_xifti = subjICmean_xifti,
    subjICvar_xifti = subjICvar_xifti,
    model_result = result
  )
  class(RESULT) <- 'templateICA.cifti'
  return(RESULT)
}


#' Activations of (s)tICA
#'
#' Identify areas of activation in each independent component map
#'
#' @param result Result of templateICA.cifti model call
#' @param spatial_model Should spatial model result be used, if available?  If FALSE, will use standard template ICA result. If NULL, use spatial model result if available.
#' @param u Activation threshold, default = 0
#' @param alpha Significance level for joint PPM, default = 0.1
#' @param type Type of region.  Default is '>' (positive excursion region).
#' @param method_p Type of multiple comparisons correction to use for p-values for standard template ICA, or NULL for no correction.  See \code{help(p.adjust)}.
#' @param verbose If \code{TRUE}, display progress of algorithm. Default: \code{FALSE}.
#' @param which.ICs Indices of ICs for which to identify activations.  If NULL, use all ICs.
#' @param deviation If \code{TRUE}. identify significant deviations from the template mean, rather than significant areas of engagement
#'
#' @return A list containing activation maps for each IC and the joint and marginal PPMs for each IC.
#'
#' @export
#'
#'
activations.cifti <- function(result, spatial_model=NULL, u=0, alpha=0.01, type=">", method_p='BH', verbose=FALSE, which.ICs=NULL, deviation=FALSE){

  if(class(result) != 'templateICA.cifti') stop("result argument must be of class 'templateICA.cifti', the result of a call to function templateICA.cifti()")

  model_result <- result$model_result
  active_xifti <- clear_data(result$subjICmean_xifti)
  nleft <- nrow(result$subjICmean_xifti$data$cortex_left)
  nright <- nrow(result$subjICmean_xifti$data$cortex_right)

  if(class(model_result)=='stICA' & is.null(spatial_model)) spatial_model <- TRUE
  if(class(model_result)=='tICA' & is.null(spatial_model)) spatial_model <- FALSE
  if((spatial_model==TRUE) & (class(model_result) == 'tICA')) {
    warning('spatial_model set to TRUE but class of model result is tICA. Setting spatial_model = FALSE, performing inference using standard template ICA.')
    spatial_model <- FALSE
  }

  #if spatial model available but spatial_model set to FALSE, grab standard template ICA model result
  if(class(model_result)=='stICA' & spatial_model==FALSE){
    model_result <- model_result$result_tICA
  }

  #run activations function
  activations_result <- activations(model_result, u=u, alpha=alpha, type=type, method_p=method_p, verbose=verbose, which.ICs=which.ICs, deviation=deviation)

  #construct xifti object for activation maps
  act_viz <- activations_result$active*1
  act_viz[act_viz==0] <- NA
  active_xifti$data$cortex_left <- act_viz[1:nleft,]
  active_xifti$data$cortex_right <- act_viz[nleft+(1:nright),]

  #include convert_to_dlabel function

  return(list(activations = activations_result,
                 active_xifti = active_xifti))
}


