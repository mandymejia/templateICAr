#' Diagnostic ICA for NIFTI
#' 
#' Run diagnostic ICA based on NIFTI-format BOLD data and NIFTI-based template
#'
#' @param nifti_fname File path of NIFTI BOLD timeseries data
#' @param templates Set of templates, each the result of call to 
#'  \code{\link{estimate_template.nifti}}.
#' (one template for each group). Either a list of objects of class 
#'  \code{"template.nifti"}.
#' @param scale Should BOLD data be scaled by the spatial standard deviation 
#'  before model fitting? Default: \code{TRUE}. If done when estimating 
#'  templates, should be done here too.
#' @param Q2 The number of nuisance ICs to identify. If \code{NULL} (default), 
#'  will be estimated. Only provide \eqn{Q2} or \eqn{maxQ} but not both.
#' @param maxQ Maximum number of ICs (template+nuisance) to identify 
#'  (\eqn{L <= maxQ <= T}). Only provide \eqn{Q2} or \eqn{maxQ} but not both.
#' @param maxiter Maximum number of EM iterations. Default: 100.
#' @param epsilon Smallest proportion change between iterations. Default: 0.01.
#' @param verbose If \code{TRUE} (default), display progress of algorithm.
#' @param out_fname The path and base name prefix of the NIFTI files to write.
#'  Will be appended with "_subjICmean.nii" for IC mean maps and 
#'  "_subjICse.nii" for IC variance maps.
#'
#' @importFrom oro.nifti readNIfTI writeNIfTI
#'
#' @return A list containing the subject IC estimates (class 'xifti'), the 
#'  subject IC variance estimates (class 'xifti'), and the result of the model
#'  call to \code{diagnosticICA} (class 'dICA')
#' 
#' @export
#'
diagnosticICA.nifti <- function(nifti_fname,
                              templates,
                              scale=TRUE,
                              Q2=NULL,
                              maxQ=100,
                              maxiter=100,
                              epsilon=0.01,
                              verbose=TRUE,
                              out_fname=NULL){

  G <- length(templates)
  if(G>5) stop(paste0('Length of templates is ',G,' which is a large number of groups. Check that the templates argument is formatted correctly.'))
  cat(paste0('Number of group templates provided: ',G,'\n'))


  # READ IN BOLD TIMESERIES DATA
  if(!file.exists(nifti_fname)) stop(paste0('The BOLD timeseries file ',nifti_fname,' does not exist.'))
  BOLD_nifti <- readNIfTI(nifti_fname, reorient=FALSE)

  # GET TEMPLATE MEAN AND VARIANCE FOR EACH GROUP
  template_class <- sapply(templates, class)
  if(any(template_class != 'template.nifti')) stop('The templates argument must be a list of template.nifti objects. See estimate_template.nifti().')
  template_mean_mat <- template_var_mat <- vector('list', length=G)

  for(g in 1:G){

    #Obtain template_mean_g and template_var_g
    template_mean_mat[[g]] <- templates[[g]]$template_mean
    template_var_mat[[g]] <- templates[[g]]$template_var

   }

  #set up mask
  masks <- vector('list', length=G)
  mask_all <- templates[[1]]$mask2
  for(g in 1:G){
    masks[[g]] <- templates[[g]]$mask2
    mask_all <- mask_all*masks[[g]]
  }
  #if masks differ over groups, re-mask templates
  for(g in 1:G){
    if(verbose) cat(paste0('Group ',g,': '))
    remask_g <- mask_all[masks[[g]]==1] #vectorized re-mask
    if(verbose) cat(paste0('removing ',sum(!remask_g),' voxels not present in all groups \n'))
    template_mean_mat[[g]] <- template_mean_mat[[g]][remask_g==1,]
    template_var_mat[[g]] <- template_var_mat[[g]][remask_g==1,]
  }

  #form BOLD data matrix
  V <- sum(mask_all)
  ntime <- dim(BOLD_nifti)[4]
  BOLD_mat <- matrix(NA, V, ntime)
  for(t in 1:ntime){
    BOLD_mat[,t] <- BOLD_nifti[,,,t][mask_all==1]
  }

  # CALL DIAGNOSTIC ICA FUNCTION
  result <- diagnosticICA(template_mean = template_mean_mat,
                            template_var = template_var_mat,
                            BOLD = BOLD_mat,
                            scale = scale,
                            Q2 = Q2,
                            maxQ = maxQ,
                            maxiter = maxiter,
                            epsilon = epsilon,
                            verbose = verbose)

  L <- ncol(result$subjICmean)
  BOLD_nifti@.Data <- BOLD_nifti@.Data[,,,1:L] #remove non-template ICs
  BOLD_nifti@dim_[5] <- L
  subjICmean_nifti <- subjICse_nifti <- BOLD_nifti #copy over header information from GICA
  img_tmp <- mask_all
  for(l in 1:L){
    img_tmp[mask_all==1] <- result$subjICmean[,l]
    subjICmean_nifti@.Data[,,,l] <- img_tmp
    img_tmp[mask_all==1] <- result$subjICse[,l]
    subjICse_nifti@.Data[,,,l] <- img_tmp
  }

  if(!is.null(out_fname)){
    out_fname_mean <- paste0(out_fname, '_subjICmean')
    out_fname_var <- paste0(out_fname, '_subjICse')
    writeNIfTI(subjICmean_nifti, out_fname_mean)
    writeNIfTI(subjICse_nifti, out_fname_var)
    writeNIfTI(mask_all, 'mask_all')
  }

  return(result)

}
