#' Estimate NIFTI template
#'
#' Estimate template for Template or Diagnostic ICA based on NIFTI-format data
#'
#' @param nifti_fnames Vector of file paths of NIFTI-format fMRI timeseries for template estimation.
#' @param nifti_fnames2 (Optional) Vector of file paths of "retest" NIFTI-format fMRI
#'  timeseries for template estimation.  Must be from the same subjects and in the same
#'  order as nifti_fnames.  Should only be provided if nifti_fnames provided, but not required.
#'  If none specified, will create pseudo test-retest data from single session.
#' @param GICA_fname File path of NIFTI-format group ICA maps (Q IC's)
#' @param mask_fname File path of NIFTI-format binary brain map
#' @param inds Indicators of which L <= Q group ICs to include in template. If NULL,
#'  use all Q original group ICs.
#' @param scale Logical indicating whether BOLD data should be scaled by the
#'  spatial standard deviation before template estimation.
#' @param verbose If \code{TRUE}. display progress updates
#' @param out_fname The path and base name prefix of the NIFTI files to write.
#' Will be appended with "_mean.nii" for template mean maps and "_var.nii" for
#' template variance maps.
#'
#' @return List of two elements: template mean of class nifti and template variance of class nifti
#'
#' @importFrom oro.nifti readNIfTI writeNIfTI
#' @importFrom matrixStats rowVars
#'
#' @export
#'
estimate_template.nifti <- function(
  nifti_fnames,
  nifti_fnames2=NULL,
  GICA_fname,
  mask_fname,
  inds=NULL,
  scale=TRUE,
  verbose=TRUE,
  out_fname=NULL){

  # Check arguments.
  if (!is.logical(scale) || length(scale) != 1) { stop('scale must be a logical value') }

  # Read GICA result
  if(verbose) cat('\n Reading in GICA result')
  GICA <- readNIfTI(GICA_fname, reorient = FALSE)
  mask2 <- mask <- readNIfTI(mask_fname, reorient = FALSE)
  dims <- dim(mask)
  vals <- sort(unique(as.vector(mask)))
  if(any(!(vals %in% c(0,1)))) stop('Mask must be binary.')
  V <- sum(mask)
  Q <- dim(GICA)[4]
  if(any(dim(GICA)[-4] != dims)) stop('GICA dims and mask dims do not match')
  GICA_mat <- matrix(NA, V, Q)
  for(q in 1:Q){
    GICA_mat[,q] <- GICA[,,,q][mask==1]
  }

  # Center each IC map.
  GICA_mat <- scale(GICA_mat, scale=FALSE)

  if(verbose){
    cat(paste0('\n Number of data locations: ',V))
    cat(paste0('\n Number of original group ICs: ',Q))
  }

  L <- Q
  if(!is.null(inds)){
    if(any(!(inds %in% 1:Q))) stop('Invalid entries in inds argument.')
    L <- length(inds)
  } else {
    inds <- 1:Q
  }

  N <- length(nifti_fnames)

  if(verbose){
    cat(paste0('\n Number of template ICs: ',L))
    cat(paste0('\n Number of training subjects: ',N))
  }

  if(!is.null(nifti_fnames2)) retest <- TRUE else retest <- FALSE
  if(retest){
    if(length(nifti_fnames) != length(nifti_fnames2)) stop('If provided, nifti_fnames2 must have same length as nifti_fnames and be in the same subject order.')
  }

  # PERFORM DUAL REGRESSION ON (PSEUDO) TEST-RETEST DATA
  DR1 <- DR2 <- array(NA, dim=c(N, L, V))
  missing_data <- NULL
  ntime_vec <- c()
  for(ii in 1:N){

    ### READ IN BOLD DATA AND PERFORM DUAL REGRESSION

    if(verbose) cat(paste0('\n Reading in data for subject ',ii,' of ',N))

    #read in BOLD
    fname_ii <- nifti_fnames[ii]
    if(!file.exists(fname_ii)) {
      missing_data <- c(missing_data, fname_ii)
      if(verbose) cat(paste0('\n Data not available'))
      next
    }
    BOLD1_ii <- readNIfTI(fname_ii, reorient = TRUE)
    ntime <- dim(BOLD1_ii)[4]
    if(any(dim(BOLD1_ii)[-4] != dims)) stop('BOLD dims and mask dims do not match')
    ntime_vec <- c(ntime_vec, ntime)
    if(length(ntime_vec)>2){
      if(var(ntime_vec) > 0) stop('All BOLD timeseries should have the same duration')
    }

    BOLD1_ii_mat <- matrix(NA, V, ntime)
    for(t in 1:ntime){
      BOLD1_ii_mat[,t] <- BOLD1_ii[,,,t][mask2==1]
    }

    if(nrow(BOLD1_ii_mat) != nrow(GICA_mat)) stop(paste0('The number of data locations in GICA and BOLD timeseries data from subject ',ii,' do not match.'))
    rm(BOLD1_ii)

    #read in BOLD retest data OR create pseudo test-retest data
    if(!retest){
      part1 <- 1:round(ntime/2)
      part2 <- setdiff(1:ntime, part1)
      BOLD2_ii_mat <- BOLD1_ii_mat[,part2]
      BOLD1_ii_mat <- BOLD1_ii_mat[,part1]
    } else {
      #read in BOLD from retest
      fname_ii <- nifti_fnames2[ii]
      if(!file.exists(fname_ii)) {
        missing_data <- c(missing_data, fname_ii)
        if(verbose) cat(paste0('\n Data not available'))
        next
      }
      BOLD2_ii <- readNIfTI(fname_ii, reorient = TRUE)
      if(any(dim(BOLD2_ii)[-4] != dims)) stop('Retest BOLD dims and mask dims do not match')
      if(dim(BOLD2_ii)[4] != ntime) stop('Retest BOLD data has different duration from first session.')
      BOLD2_ii_mat <- matrix(NA, V, ntime)
      for(t in 1:ntime){
        BOLD2_ii_mat[,t] <- BOLD2_ii[,,,t][mask2==1]
      }
      rm(BOLD2_ii)
      if(nrow(BOLD2_ii_mat) != nrow(GICA_mat)) stop(paste0('The number of data locations in GICA and BOLD retest timeseries data from subject ',ii,' do not match.'))
    }

    flat_vox <- ((rowVars(BOLD1_ii_mat)==0) | (rowVars(BOLD2_ii_mat)==0))
    if(sum(flat_vox)>0) {
      warning(paste0(sum(flat_vox), ' flat voxels detected. Removing these from the mask for this and future subjects. Updated mask will be returned with estimated templates.'))
      mask2[mask2==1] <- (!flat_vox)
      GICA_mat <- GICA_mat[!flat_vox,]
      BOLD1_ii_mat <- BOLD1_ii_mat[!flat_vox,]
      BOLD2_ii_mat <- BOLD2_ii_mat[!flat_vox,]
      DR1 <- DR1[,,!flat_vox]
      DR2 <- DR2[,,!flat_vox]
      V <- sum(mask2)
    }

    #perform dual regression on test and retest data
    DR1_ii <- dual_reg(BOLD1_ii_mat, GICA_mat, scale=scale)$S
    DR2_ii <- dual_reg(BOLD2_ii_mat, GICA_mat, scale=scale)$S
    DR1[ii,,] <- DR1_ii[inds,]
    DR2[ii,,] <- DR2_ii[inds,]
  }

  cat(paste0('Total number of voxels in updated mask: ', V, '\n'))

  # Estimate template
  if (verbose) { cat("\nEstimating template.\n") }
  template <- estimate_template_from_DR(DR1, DR2, var_method="unbiased")
  rm(DR1, DR2)

  if(!is.null(out_fname)){
    out_fname_mean <- paste0(out_fname, '_mean')
    out_fname_var <- paste0(out_fname, '_var')
    GICA@.Data <- GICA@.Data[,,,inds] #remove non-template ICs
    GICA@dim_[5] <- length(inds)
    template_mean_nifti <- template_var_nifti <- GICA #copy over header information from GICA
    img_tmp <- mask2
    for(l in 1:L){
      img_tmp[mask2==1] <- template_mean_nifti[,l]
      template_mean_nifti@.Data[,,,l] <- img_tmp
      img_tmp[mask2==1] <- template_var[,l]
      template_var_nifti@.Data[,,,l] <- img_tmp
    }
    writeNIfTI(template_mean_nifti, out_fname_mean)
    writeNIfTI(template_var_nifti, out_fname_var)
    writeNIfTI(mask2, 'mask2')
  }

  result <- list(template_mean=template_mean_nifti, template_var=template_var_nifti, scale=scale, mask=mask, mask2=mask2, inds=inds)
  class(result) <- 'template.nifti'
  return(result)
}
