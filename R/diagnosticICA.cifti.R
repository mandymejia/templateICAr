#' Diagnostic ICA for CIFTI
#' 
#' Run diagnostic ICA based on CIFTI-format BOLD data and CIFTI-based template
#'
#' @param cifti_fname File path of CIFTI-format timeseries data (ending in 
#'  ".dtseries.nii").
#' @param templates Set of templates, each the result of call to 
#'  \code{\link{estimate_template.cifti}} (one template for each group).
#' Either a list of objects of class \code{"template.cifti"}, OR a vector of 
#'  file path prefixes (the part before "_mean.dscalar.nii"
#'  and "_var.dscalar.nii") of cifti files written out by 
#'  \code{\link{estimate_template.cifti}}.
#' @param spatial_model Should spatial modeling be performed? Default: 
#'  \code{FALSE}. If \code{TRUE}, surface geometries will be used to fit a 
#'  spatial Bayesian model. Computational demands (time and memory) are much 
#'  higher with \code{spatial_model=TRUE} than \code{spatial_model=FALSE}, but 
#'  results are more accurate.
#' @param brainstructures Character vector indicating which brain structure(s)
#'  to obtain: \code{"left"} (left cortical surface), \code{"right"} (right
#'  cortical surface) and/or \code{"subcortical"} (subcortical and cerebellar
#'  gray matter). Can also be \code{"all"} (obtain all three brain structures).
#'  Default: \code{c("all")}.  For spatial
#'  modeling (performed if \code{surfL_fname} and/or \code{surfR_fname} provided), 
#'  subcortical locations will be ignored.
#' @param surfL_fname (Required for spatial modeling) File path of GIFTI surface geometry
#'  file representing the left cortex. Only required if \code{brainstructures} includes \code{"right"}.
#' @param surfR_fname (Required for spatial modeling) File path of GIFTI surface geometry
#'  file representing the right cortex. Only required if \code{brainstructures} includes \code{"left"}.
#' @param resamp_res (Only recommended for spatial modeling) Target resolution for resampling (number of
#'  cortical surface vertices per hemisphere). A value less than 10000 is recommended for computational
#'  feasibility. If \code{NULL} (default) or \code{FALSE}, do not perform resampling.
#' @param scale Should BOLD data be scaled by the spatial standard deviation 
#'  before model fitting? Default: \code{TRUE}. If done when estimating 
#'  templates, should be done here too.
#' @param Q2 The number of nuisance ICs to identify. If \code{NULL} (default), 
#'  will be estimated. Only provide \eqn{Q2} or \eqn{Q2_max} but not both.
#' @param Q2_max Maximum number of ICs (template+nuisance) to identify 
#'  (\eqn{L <= Q2_max <= T}). Only provide \eqn{Q2} or \eqn{Q2_max} but not both.
#' @param maxiter Maximum number of EM iterations. Default: 100.
#' @param epsilon Smallest proportion change between iterations. Default: 0.01.
#' @param verbose If \code{TRUE} (default), display progress of algorithm.
#' @param kappa_init Starting value for kappa.  If \code{NULL}, starting value 
#'  will be determined automatically. Default: 0.4.
#' @param write_dir Where should any output files be written? \code{NULL} 
#'  (default) will write them to the current working directory.
#'
# @importFrom INLA inla.pardiso.check inla.setOption
#'
#' @return A list containing the subject IC estimates (class \code{"xifti"}), the 
#'  subject IC variance estimates (class \code{"xifti"}), and the result of the model 
#'  call to \code{diagnosticICA} (class \code{"dICA"}).
#' 
#' @importFrom fMRItools match_input
#' @keywords internal
#'
diagnosticICA.cifti <- function(cifti_fname,
                              templates,
                              brainstructures="all",
                              spatial_model=FALSE,
                              surfL_fname=NULL,
                              surfR_fname=NULL,
                              resamp_res=NULL,
                              scale=TRUE,
                              Q2=NULL,
                              Q2_max=NULL,
                              maxiter=100,
                              epsilon=0.01,
                              verbose=TRUE,
                              kappa_init=0.2,
                              write_dir=NULL){

  if (!requireNamespace("ciftiTools", quietly = TRUE)) {
    stop("Package \"ciftiTools\" needed to read NIFTI data. Please install it.", call. = FALSE)
  }

  if (is.null(write_dir)) { write_dir <- getwd() }

  brainstructures <- fMRItools::match_input(
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
    if(!is.character(templates[[1]])) stop('If spatial_model=TRUE, template argument must be file path prefix to cifti files written by estimate_template.cifti().')
    flag <- INLA::inla.pardiso.check()
    if (any(grepl('FAILURE',flag))) stop('PARDISO IS NOT INSTALLED OR NOT WORKING. PARDISO for R-INLA is required for computational efficiency. If you already have a PARDISO / R-INLA License, run inla.setOption(pardiso.license = "/path/to/license") and try again.  If not, run inla.pardiso() to obtain a license.')
    INLA::inla.setOption(smtp='pardiso')
  }

  if(!spatial_model){
    if(!is.null(resamp_res)) warning('Resampling only necessary for spatial modeling.  Setting resamp_res to NULL.')
    resamp_res <- NULL
  }

  G <- length(templates)
  if(G>5) stop(paste0('Length of templates is ',G,' which is a large number of groups. Check that the templates argument is formatted correctly.'))
  cat(paste0('Number of group templates provided: ',G,'\n'))


  # READ IN BOLD TIMESERIES DATA
  if(!file.exists(cifti_fname)) stop(paste0('The BOLD timeseries file ',cifti_fname,' does not exist.'))
  BOLD <- ciftiTools::read_cifti(cifti_fname,
                    surfL_fname = surfL_fname,
                    surfR_fname = surfR_fname,
                    brainstructures = brainstructures,
                    resamp_res=resamp_res)
  locs_left <- BOLD_cifti$meta$cortex$medial_wall_mask$left
  locs_right <- BOLD_cifti$meta$cortex$medial_wall_mask$right
  # Separate data and metadata
  BOLD_meta <- BOLD; BOLD_meta["data"] <- list(NULL)
  BOLD <- as.matrix(BOLD)

  # GET TEMPLATE MEAN AND VARIANCE FOR EACH GROUP
  template_class <- sapply(templates, class)
  if(length(unique(template_class))>1) stop('All elements of templates argument must be of the same class')
  template_mean_cifti <- template_var_cifti <- vector('list', length=G)

  for(g in 1:G){

    #Obtain template_mean_g and template_var_g
    if (template_class[1]=='template.cifti') {
      template_mean_cifti[[g]] <- templates[[g]]$template_mean #class xifti
      template_var_cifti[[g]] <- templates[[g]]$template_var #class xifti
    } else if (template_class[1]=='character') {
      fname_mean <- paste0(templates[[g]],'_mean.dscalar.nii')
      fname_var <- paste0(templates[[g]],'_var.dscalar.nii')
      if(!file.exists(fname_mean)) stop(paste0('The file ', fname_mean, ' does not exist.'))
      if(!file.exists(fname_var)) stop(paste0('The file ', fname_var, ' does not exist.'))
      template_mean_cifti[[g]] <- ciftiTools::read_cifti(fname_mean, surfL_fname = surfL_fname, surfR_fname = surfR_fname, brainstructures=brainstructures, resamp_res=resamp_res)
      template_var_cifti[[g]] <- ciftiTools::read_cifti(fname_var, surfL_fname = surfL_fname, surfR_fname = surfR_fname, brainstructures=brainstructures, resamp_res=resamp_res)
    } else {
      stop('template argument must be an object of class template.cifti or file path prefix to result of estimate_template.cifti() (same as out_fname argument passed to estimate_template.cifti().')
    }

    #Extract data matrix from template_mean_g and template_var_g
    #template_mean_mat[[g]] <- rbind(template_mean_g$data$cortex_left, template_mean_g$data$cortex_right, template_mean_g$data$subcort)
    #template_var_mat[[g]] <- rbind(template_var_g$data$cortex_left, template_var_g$data$cortex_right, template_var_g$data$subcort)

    # DELETE THIS WITH ciftiTools version 1.5
    #locs_left <- (!is.nan(template_mean_cifti[[g]]$data$cortex_left[,1]) & !is.nan(template_var_cifti[[g]]$data$cortex_left[,1]))
    #locs_right <- (!is.nan(template_mean_cifti[[g]]$data$cortex_right[,1]) & !is.nan(template_var_cifti[[g]]$data$cortex_right[,1]))
    #template_mean_cifti[[g]]$data$cortex_left <- template_mean_cifti[[g]]$data$cortex_left[locs_left,]
    #template_var_cifti[[g]]$data$cortex_left <- template_var_cifti[[g]]$data$cortex_left[locs_left,]
    #template_mean_cifti[[g]]$data$cortex_right <- template_mean_cifti[[g]]$data$cortex_right[locs_right,]
    #template_var_cifti[[g]]$data$cortex_right <- template_var_cifti[[g]]$data$cortex_right[locs_right,]
    #template_mean_cifti[[g]]$meta$cortex$medial_wall_mask$left <- template_var_cifti[[g]]$meta$cortex$medial_wall_mask$left <- locs_left
    #template_mean_cifti[[g]]$meta$cortex$medial_wall_mask$right <- template_var_cifti[[g]]$meta$cortex$medial_wall_mask$right <- locs_right

    locs_left <- locs_left & template_mean_cifti[[g]]$meta$cortex$medial_wall_mask$left & template_var_cifti[[g]]$meta$cortex$medial_wall_mask$left
    locs_right <- locs_right & template_mean_cifti[[g]]$meta$cortex$medial_wall_mask$right & template_var_cifti[[g]]$meta$cortex$medial_wall_mask$right

  }

  # DELETE THIS WITH ciftiTools version 1.5, which will hopefully result in the same medial wall for BOLD and the templates (only a problem with resampling)
  ntime <- ncol(BOLD)
  #left cortex
  dat_all <- matrix(NA, nrow=length(locs_left), ncol=ntime)
  dat_all[BOLD_meta$meta$cortex$medial_wall_mask$left,] <- BOLD$data$cortex_left
  BOLD$data$cortex_left <- dat_all[locs_left,,drop=FALSE]
  BOLD_cifti$meta$cortex$medial_wall_mask$left <- locs_left
  #right cortex
  dat_all <- matrix(NA, nrow=length(locs_right), ncol=ntime)
  dat_all[BOLD_cifti$meta$cortex$medial_wall_mask$right,] <- BOLD$data$cortex_right
  BOLD$data$cortex_right <- dat_all[locs_right,,drop=FALSE]
  BOLD_cifti$meta$cortex$medial_wall_mask$right <- locs_right

  #TO DO: SINGLE SPATIAL MODEL

  # # IF SPATIAL MODELING, LOOP OVER HEMISPHERES
  # if(spatial_model) {
  #   if(do_left & !do_right) models <- c('lh')
  #   if(do_right & !do_left) models <- c('rh')
  #   if(do_left & do_right) models <- c('lh','rh')
  # } else {
  #   models <- 'single'
  # }

  #models_list <- vector('list', length=length(models))
  #names(models_list) <- models
  #for(mod in models){

  # IF SPATIAL MODELING, CONSTRUCT MESH
  if(spatial_model){
    #num_surfs <- do_left + do_right # never used
    meshes <- vector('list', 2)
    # ind <- 1 # never used
    if(do_left) {
      surf <- BOLD_meta$surf$cortex_left;
      wall_mask <- which(BOLD_meta$meta$cortex$medial_wall_mask$left)
      meshes[[1]] <- make_mesh(surf = surf, inds_mesh = wall_mask) #remove wall for greater computational efficiency
    }
    if(do_right) {
      surf <- BOLD_meta$surf$cortex_right;
      wall_mask <- which(BOLD_meta$meta$cortex$medial_wall_mask$right)
      meshes[[2]] <- make_mesh(surf = surf, inds_mesh = wall_mask) #remove wall for greater computational efficiency
    }
  } else {
    meshes <- NULL
  }

  # CALL DIAGNOSTIC ICA FUNCTION
  result <- diagnosticICA(template_mean = lapply(template_mean_cifti, as.matrix),
                            template_var = lapply(template_var_cifti, as.matrix),
                            BOLD = BOLD,
                            scale = scale,
                            meshes = meshes,
                            Q2 = Q2,
                            Q2_max = Q2_max,
                            maxiter = maxiter,
                            epsilon = epsilon,
                            verbose = verbose,
                            kappa_init = kappa_init)

  #HERE (BELOW IS FOR SINGLE MODEL.  ACTUALLY FOR SPATIAL MODEL NEED TO GENERALIZE ABOVE TO SINGLE MODEL.)

  result$subjICmean <- ciftiTools::newdata_xifti(
    ciftiTools::select_xifti(template_mean_cifti[[1]], rep(1, ncol(result$subjICmean))),
    result$subjICmean
  )
  result$subjICse <- ciftiTools::newdata_xifti(
    ciftiTools::select_xifti(template_var_cifti[[1]], rep(1, ncol(result$subjICse))),
    result$subjICse
  )

  # RETURN XIFTI RESULTS AND MODEL RESULT
  result
}
