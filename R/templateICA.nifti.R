#' Template ICA for NIFTI-format data
#'
#' Run template ICA based on NIFTI-format BOLD data and NIFTI-based template
#'
#' @param BOLD_fname File path of NIFTI-format timeseries data (ending in .nii or .nii.gz).
#' @param template Result of call to \code{estimate_template.nifti}, either (1) an object of class
#' \code{template.nifti} OR (2) the file path and basename prefix (e.g., the part before "_mean.nii"
#' of NIFTI files written out by \code{estimate_template.nifti}.
#' @param mask_fname (Optional) File path to brain mask (0/1 or logical). Will also try to use the brain mask
#' output by \code{estimate_template.nifti} if available.
#' @param scale Logical indicating whether BOLD data should be scaled by the spatial
#' standard deviation before model fitting. If done when estimating templates, should be done here too.
#' @param Q2 The number of nuisance ICs to identify. If NULL, will be estimated. Only provide \eqn{Q2} or \eqn{maxQ} but not both.
#' @param maxQ Maximum number of ICs (template+nuisance) to identify (L <= maxQ <= T). Only provide \eqn{Q2} or \eqn{maxQ} but not both.
#' @param maxiter Maximum number of EM iterations. Default: 100.
#' @param epsilon Smallest proportion change between iterations. Default: 0.001.
#' @param verbose If \code{TRUE} (default), display progress of algorithm.
#' @param write_dir Where should any output files be written? \code{NULL} (default) will write them to the current working directory.
#' @param time_inds Subset of fMRI BOLD volumes to include in analysis. If NULL (default), use all volumes.
#'
#' @importFrom oro.nifti readNIfTI writeNIfTI
#' @importFrom matrixStats rowVars
#'
#' @return A list containing the subject IC estimates (class 'xifti'), the subject IC variance estimates (class 'xifti'), and the result of the model call to \code{templateICA} (class 'dICA')
#'
#' @export
#'
templateICA.nifti <- function(BOLD_fname,
                              template,
                              mask_fname=NULL,
                              scale=TRUE,
                              Q2=NULL,
                              maxQ=NULL,
                              maxiter=100,
                              epsilon=0.001,
                              verbose=TRUE,
                              write_dir=NULL,
                              time_inds=NULL){

  if (is.null(write_dir)) { write_dir <- getwd() }

  # READ IN BOLD TIMESERIES DATA
  if(!file.exists(BOLD_fname)) stop(paste0('The BOLD timeseries file ',BOLD_fname,' does not exist.'))
  BOLD_nifti <- readNIfTI(BOLD_fname, reorient=FALSE)

  # GET TEMPLATE MEAN AND VARIANCE
  template_class <- class(template)
  if(template_class=='template.nifti'){
    template_mean <- template$template_mean
    template_var <- template$template_var
    template_mask <- template$mask2
  } else if(is.character(template)){
    if(verbose) cat('Reading in templates.\n')
    fname_mean <- paste0(template,'_mean.nii')
    fname_var <- paste0(template,'_var.nii')
    fname_mask <- paste0(template,'_mask.nii')
    if(!file.exists(fname_mean)) stop(paste0('The file ', fname_mean, ' does not exist.'))
    if(!file.exists(fname_var)) stop(paste0('The file ', fname_var, ' does not exist.'))
    if(!file.exists(fname_mask)) stop(paste0('The file ', fname_mask, ' does not exist.'))
    if(verbose) cat(paste0('Reading in template mean file: ',fname_mean, '\n'))
    template_mean <- readNIfTI(fname_mean, reorient = FALSE)
    if(verbose) cat(paste0('Reading in template variance file: ',fname_var, '\n'))
    template_var <- readNIfTI(fname_var, reorient = FALSE)
    if(verbose) cat(paste0('Reading in template mask: ',fname_mask, '\n'))
    template_mask <- readNIfTI(fname_mask, reorient = FALSE)
  } else {
    stop('template argument must be an object of class template.nifti or file path prefix to result of estimate_template.nifti() (matches out_fname argument passed to estimate_template.nifti().')
  }

  # READ IN MASK IF PROVIDED
  if(!is.null(mask_fname)){
    mask <- readNIfTI(mask_fname, reorient=FALSE)
    if(sum(mask==TRUE, na.rm=TRUE) + sum(mask==FALSE, na.rm=TRUE) != length(mask)) stop('Input mask contains values that cannot be coerced to logical.')
    if(!all.equal(dim(mask), dim(template_mask))) stop('Dimensions of template mask and provided mask do not match.')
    } else {
    mask <- template_mask
  }
  if(sum(template_mask==TRUE, na.rm=TRUE) + sum(template_mask==FALSE, na.rm=TRUE) != length(template_mask)) stop('Template mask contains values that cannot be coerced to logical.')
  mask <- as.logical(mask*template_mask) #intersection mask
  V <- sum(mask==TRUE)
  if(verbose) cat(paste0('Mask contains ', V, ' locations to be included in analysis.'))

  # CHECK DIMENSIONS
  dim3d <- dim(mask)[1:3]
  Q <- dim(template_mean)[4]
  ntime <- dim(BOLD_nifti)[4]
  if(!all.equal(dim3d, dim(template_mean)[1:3])) stop('Dimensions of template mean image do not match mask dimensions')
  if(!all.equal(dim3d, dim(template_var)[1:3])) stop('Dimensions of template var image do not match mask dimensions')
  if(!all.equal(dim3d, dim(BOLD_nifti)[1:3])) stop('Dimensions of BOLD image do not match mask dimensions')
  if(dim(template_var)[4] != Q) stop('Number of ICs in template mean and template var images do not match.')

  if(verbose) cat(paste0('Number of template ICs: ', Q, '\n'))
  if(verbose) cat(paste0('Length of BOLD timeseries: ', ntime, '\n'))

  # FORM DATA MATRIX AND TEMPLATE MATRICES
  vectorize <- function(img3d, mask){ return(img3d[mask=TRUE]) }
  BOLD_mat <- apply(BOLD_nifti, 4, vectorize)
  template_mean_mat <- apply(template_mean, 4, vectorize)
  template_var_mat <- apply(template_var, 4, vectorize)

  if(!is.null(time_inds)){
    all_inds <- 1:ncol(BOLD_mat)
    if(any(!(time_inds %in% all_inds))) stop(paste0('time_inds contains indices outside of the range of 1 to ', max(all_inds)))
    BOLD_mat <- BOLD_mat[,time_inds]
  }

  # CHECK FOR FLAT TIMESERIES
  flat_vox <- (rowVars(BOLD_mat)==0)
  if(sum(flat_vox)>0) {
    warning(paste0(sum(flat_vox), ' flat voxels detected. Removing these from the mask. Updated mask will be returned.'))
    mask[mask==1] <- (!flat_vox)
    template_mean_mat <- template_mean_mat[!flat_vox,]
    template_var_mat <- template_var_mat[!flat_vox,]
    BOLD_mat <- BOLD_mat[!flat_vox,]
    V <- sum(mask)
  }



  # CALL TEMPLATE ICA FUNCTION

  result <- templateICA(template_mean = template_mean_mat,
                            template_var = template_var_mat,
                            BOLD = BOLD_mat,
                            scale = scale,
                            Q2 = Q2,
                            maxQ=maxQ,
                            maxiter=maxiter,
                            epsilon=epsilon,
                            verbose=verbose)

  # HERE -----

  # MAP ESTIMATES AND VARIANCE TO NIFTI FORMAT

  #reformat column names to IC 1, IC 2, etc.
  if(!grepl('IC',template_mean$meta$cifti$names[1])){
    template_mean$meta$cifti$names <- paste0('IC ',template_mean$meta$cifti$names)
  }

  subjICmean_xifti <- subjICvar_xifti <- template_mean
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

#return mask!!

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


