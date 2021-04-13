#' Estimate CIFTI template
#'
#' Estimate template for Template or Diagnostic ICA based on CIFTI-format data
#'
#' @param cifti_fnames Vector of file paths of CIFTI-format fMRI timeseries
#'  (*.dtseries.nii) for template estimation
#' @param cifti_fnames2 Vector of file paths of "retest" CIFTI-format fMRI
#'  timeseries (*.dtseries.nii) for template estimation.  Must be from the same
#'  subjects and in the same order as cifti_fnames.  If none specified, will
#'  create pseudo test-retest data from single session.
#' @param GICA_fname File path of CIFTI-format group ICA maps (ending in .d*.nii)
#' @param inds Indicators of which group ICs to include in template. If NULL,
#'  use all group ICs.
#' @param scale Logical indicating whether BOLD data should be scaled by the
#'  spatial standard deviation before template estimation.
#' @param brainstructures Character vector indicating which brain structure(s)
#'  to obtain: \code{"left"} (left cortical surface), \code{"right"} (right
#'  cortical surface) and/or \code{"subcortical"} (subcortical and cerebellar
#'  gray matter). Can also be \code{"all"} (obtain all three brain structures).
#'  Default: \code{c("left","right")} (cortical surface only).
#' @param verbose If \code{TRUE}. display progress updates
#' @param out_fname (Required if templates are to be resampled to a lower spatial
#' resolution, usually necessary for spatial template ICA.) The path and base name
#' prefix of the CIFTI files to write. Will be appended with "_mean.dscalar.nii" for
#' template mean maps and "_var.dscalar.nii" for template variance maps.
#' @param FC If TRUE, include template for functional connectivity
#'
#' @importFrom ciftiTools read_cifti write_cifti
#'
#' @return List of two elements: template mean of class xifti and template
#'  variance of class xifti
#'
#' @export
#'
estimate_template.cifti <- function(
  cifti_fnames,
  cifti_fnames2=NULL,
  GICA_fname,
  inds=NULL,
  scale=TRUE,
  brainstructures=c("left","right"),
  verbose=TRUE,
  out_fname=NULL,
  FC = FALSE){

  #TO DOs:
  # Create function to print and check template, template.cifti and template.nifti objects

  #IDEA:
  # Denoise the data by estimating and removing noise IC's for each scan, as we do in template ICA
  # Have an argument denoise (logical).  For small datasets, denoising may result in cleaner templates.
  # For larger datasets, denoising may be slow and unnecessary.

  if(!is.null(out_fname)){
    if(!dir.exists(dirname(out_fname))) stop('directory part of out_fname does not exist')
  }

  retest <- !is.null(cifti_fnames2)
  if(retest){
    if(length(cifti_fnames) != length(cifti_fnames2)) stop('If provided, cifti_fnames2 must have same length as cifti_fnames and be in the same subject order.')
  }

  notthere <- sum(!file.exists(cifti_fnames))
  if(notthere == length(cifti_fnames)) stop('The files in cifti_fnames do not exist.')
  if(notthere > 0) warning(paste0('There are ', notthere, ' files in cifti_fnames that do not exist. These scans will be excluded from template estimation.'))
  if(retest) {
    notthere2 <- sum(!file.exists(cifti_fnames2))
    if(notthere2 == length(cifti_fnames2)) stop('The files in cifti_fnames2 do not exist.')
    if(notthere2 > 0) warning(paste0('There are ', notthere2, ' files in cifti_fnames2 that do not exist. These scans will be excluded from template estimation.'))
  }

  # Check arguments.
  if (!is.logical(scale) || length(scale) != 1) { stop('scale must be a logical value') }
  brainstructures <- match_input(
    brainstructures, c("left","right","subcortical","all"),
    user_value_label="brainstructures"
  )
  if ("all" %in% brainstructures) {
    brainstructures <- c("left","right","subcortical")
  }

  # Read GICA result
  if(verbose) { cat('\n Reading in GICA result') }
  GICA <- read_cifti(GICA_fname, brainstructures=brainstructures)
  GICA_flat <- do.call(rbind, GICA$data)
  V <- nrow(GICA_flat); Q <- ncol(GICA_flat)
  # Center each IC map.
  GICA_flat <- scale(GICA_flat, scale=FALSE)
  if(verbose) {
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

  N <- length(cifti_fnames)

  if(verbose){
    cat(paste0('\n Number of template ICs: ',L))
    cat(paste0('\n Number of training subjects: ',N))
  }

  # # Obtain the brainstructure mask for the flattened CIFTIs.
  # flat_bs_mask <- vector("logical", 0)
  # if ("left" %in% brainstructures) {
  #   left_mwall <- GICA$meta$cortex$medial_wall_mask$left
  #   if (is.null(left_mwall)) {
  #     left_mwall <- rep(TRUE, nrow(GICA$data$cortex_left))
  #   }
  #   flat_bs_mask <- c(flat_bs_mask, ifelse(left_mwall, "left", "mwall"))
  # }
  # if ("right" %in% brainstructures) {
  #   right_mwall <- GICA$meta$cortex$medial_wall_mask$right
  #   if (is.null(right_mwall)) {
  #     right_mwall <- rep(TRUE, nrow(GICA$data$cortex_right))
  #   }
  #   flat_bs_mask <- c(flat_bs_mask, ifelse(right_mwall, "right", "mwall"))
  # }
  # if ("subcortical" %in% brainstructures) {
  #   flat_bs_mask <- c(flat_bs_mask, rep("subcortical", nrow(GICA$data$subcort)))
  # }


  # PERFORM DUAL REGRESSION ON (PSEUDO) TEST-RETEST DATA
  DR1 <- DR2 <- array(NA, dim=c(N, L, V))
  if(FC) FC1 <- FC2 <- array(NA, dim=c(N, L, L)) #for functional connectivity template
  missing_data <- NULL
  for(ii in 1:N){

    if(verbose) cat(paste0('\n Reading and analyzing data for subject ',ii,' of ',N))

    #read in BOLD
    fname_ii <- cifti_fnames[ii]
    if(!file.exists(fname_ii)) {
      missing_data <- c(missing_data, fname_ii)
      if(verbose) cat(paste0('\n Data not available for file:', fname_ii))
      next
    }

    BOLD1_ii <- do.call(rbind, read_cifti(fname_ii, brainstructures=brainstructures)$data)
    if(nrow(BOLD1_ii) != nrow(GICA_flat)) {
      stop(paste0(
        'The number of data locations in GICA and',
        'BOLD data from subject ', ii,' do not match.'
      ))
    }
    ntime <- ncol(BOLD1_ii)

    #read in BOLD retest data OR create pseudo test-retest data
    if(!retest){
      part1 <- 1:round(ntime/2)
      part2 <- setdiff(1:ntime, part1)
      BOLD2_ii <- BOLD1_ii[,part2]
      BOLD1_ii <- BOLD1_ii[,part1]
    } else {
      #read in BOLD from retest
      fname_ii <- cifti_fnames2[ii]
      if(!file.exists(fname_ii)) {
        missing_data <- c(missing_data, fname_ii)
        if(verbose) cat(paste0('\n Data not available for this file:', fname_ii))
        next
      }
      BOLD2_ii <- do.call(rbind, read_cifti(fname_ii, brainstructures=brainstructures)$data)
      if (nrow(BOLD2_ii) != nrow(GICA_flat)) {
        stop(paste0(
          'The number of data locations in GICA and',
          'BOLD data from subject ', ii,' do not match.'
        ))
      }
    }

    #perform dual regression on test and retest data
    DR1_ii <- dual_reg(BOLD1_ii, GICA_flat, scale=scale)
    DR2_ii <- dual_reg(BOLD2_ii, GICA_flat, scale=scale)
    DR1[ii,,] <- DR1_ii$S[inds,]
    DR2[ii,,] <- DR2_ii$S[inds,]
    if(FC) FC1[ii,,] <- cov(DR1_ii$A[,inds])
    if(FC) FC2[ii,,] <- cov(DR2_ii$A[,inds])
  }
  #save(FC1, FC2, file='FC_temp.RData')

  # ESTIMATE MEAN

  if(verbose) cat('\n Estimating Template Mean')
  mean1 <- apply(DR1, c(2,3), mean, na.rm=TRUE)
  mean2 <- apply(DR2, c(2,3), mean, na.rm=TRUE)
  template_mean <- t((mean1 + mean2)/2)

  # ESTIMATE SIGNAL (BETWEEN-SUBJECT) VARIANCE

  # mean_lmer <- var_noise_lmer <- var_signal_lmer <- matrix(NA, V, L)
  #
  # ids <- c(1:N, 1:N)
  # for(v in 1:V){
  #   print(v)
  #   for(l in 1:L){
  #     #fit an lmer
  #     DR_vl <- c(DR1[,l,v], DR2[,l,v])
  #     #bad <- is.na(DR1[,l,v]) | is.na(DR2[,l,v]) #exclude any subjects with any missing data
  #     lmer_vl <- lmer(DR_vl ~ (1|ids))#, subset = c(!bad, !bad))
  #     mean_lmer[v,l] <- summary(lmer_vl)$coefficients[1,1]
  #     var_noise_lmer[v,l] <- (summary(lmer_vl)$sigma)^2
  #     var_signal_lmer[v,l] <- VarCorr(lmer_vl)$ids[1]
  #   }
  # }

  # total variance
  if(verbose) cat('\n Estimating Total Variance')
  var_tot1 <- apply(DR1, c(2,3), var, na.rm=TRUE)
  var_tot2 <- apply(DR2, c(2,3), var, na.rm=TRUE)
  var_tot <- t((var_tot1 + var_tot2)/2)

  # noise (within-subject) variance
  if(verbose) cat('\n Estimating Within-Subject Variance')
  DR_diff <- DR1 - DR2;
  var_noise <- t((1/2)*apply(DR_diff, c(2,3), var, na.rm=TRUE))

  # signal (between-subject) variance
  if(verbose) cat('\n Estimating Template (Between-Subject) Variance \n')
  template_var <- var_tot - var_noise
  template_var[template_var < 0] <- 0

  rm(DR1, DR2, mean1, mean2, var_tot1, var_tot2, var_tot, DR_diff)

  # Format template as "xifti"s

  xifti_mean <- xifti_var <- GICA
  nleft <- nrow(GICA$data$cortex_left)
  nright <- nrow(GICA$data$cortex_right)
  nsub <- nrow(GICA$data$subcort)
  if ("left" %in% brainstructures) {
    xifti_mean$data$cortex_left <- template_mean[1:nleft,]  #[flat_bs_mask == "left",, drop=FALSE]
    xifti_var$data$cortex_left <- template_var[1:nleft,] #[flat_bs_mask == "left",, drop=FALSE]
    #xifti_noisevar$data$cortex_left <- var_noise[flat_bs_mask == "left",, drop=FALSE]
  }
  if ("right" %in% brainstructures) {
    xifti_mean$data$cortex_right <- template_mean[nleft+(1:nright),] #[flat_bs_mask == "right",, drop=FALSE]
    xifti_var$data$cortex_right <- template_var[nleft+(1:nright),] #[flat_bs_mask == "right",, drop=FALSE]
    #xifti_noisevar$data$cortex_right <- var_noise[flat_bs_mask == "right",, drop=FALSE]
  }
  if ("subcortical" %in% brainstructures) {
    xifti_mean$data$subcort <- template_mean[nleft+nright+(1:nsub),] #[flat_bs_mask == "subcortical",, drop=FALSE]
    xifti_var$data$subcort <- template_var[nleft+nright+(1:nsub),] #[flat_bs_mask == "subcortical",, drop=FALSE]
  }

  xifti_mean$meta$cifti$names <- paste0('IC ',inds)
  xifti_var$meta$cifti$names <- paste0('IC ',inds)

  if(!is.null(out_fname)){
    out_fname_mean <- paste0(out_fname, '_mean.dscalar.nii')
    out_fname_var <- paste0(out_fname, '_var.dscalar.nii')
    write_cifti(xifti_mean, out_fname_mean, verbose=verbose)
    write_cifti(xifti_var, out_fname_var, verbose=verbose)
  }

  if(FC){

    mean_FC <- (apply(FC1, c(2,3), mean, na.rm=TRUE) + apply(FC2, c(2,3), mean, na.rm=TRUE))/2
    var_FC_tot  <- (apply(FC1, c(2,3), var, na.rm=TRUE) + apply(FC2, c(2,3), var, na.rm=TRUE))/2
    var_FC_within  <- 1/2*(apply(FC1-FC2, c(2,3), var, na.rm=TRUE))
    var_FC_between <- var_FC_tot - var_FC_within
    var_FC_between[var_FC_between < 0] <- NA

    #function to minimize w.r.t. k
    fun <- function(nu, p, var_ij, xbar_ij, xbar_ii, xbar_jj){
      LHS <- var_ij
      phi_ij <- xbar_ij*(nu-p-1)
      phi_ii <- xbar_ii*(nu-p-1)
      phi_jj <- xbar_jj*(nu-p-1)
      RHS_numer <- (nu-p+1)*phi_ij^2 + (nu-p-1)*phi_ii*phi_jj
      RHS_denom <- (nu-p)*((nu-p-1)^2)*(nu-p-3)

      sq_diff <- (LHS - RHS_numer/RHS_denom)^2
      return(sq_diff)
    }

    nu_est <- matrix(NA, L, L)
    for(q1 in 1:(L-1)){
      for(q2 in (q1+1):L){

        #estimate k = nu - p - 1
        nu_opt <- optimize(f=fun, interval=c(L+1,100), p=L, var_ij=var_FC_between[q1,q2], xbar_ij=mean_FC[q1,q2], xbar_ii=mean_FC[q1,q1], xbar_jj=mean_FC[q2,q2])
        nu_est[q1,q2] <- nu_opt$minimum
      }
    }
    nu_est[lower.tri(nu_est, diag=TRUE)] <- NA
    nu_est_avg <- mean(nu_est[upper.tri(nu_est)], na.rm=TRUE)
    template_FC <- list(nu = nu_est_avg,
                        psi = mean_FC/(nu_est_avg - L - 1))
  } else {
    template_FC <- NULL
  }

  var2 <- matrix(NA, L, L)
  for(q1 in 1:(L-1)){
    for(q2 in (q1+1):L){
      nu <- nu_est_avg
      num_qq <- (nu-L+1)*mean_FC[q1,q2]*(nu-L-1)^2 + (nu - L - 1)*1*1*(nu - L - 1)^2
      denom_qq <- (nu-L)*((nu-L-1)^2)*(nu-L-3)
      var2[q1,q2] <- num_qq/denom_qq
    }
  }

  result <- list(template_mean=xifti_mean, template_var=xifti_var, scale=scale, inds=inds)
  class(result) <- 'template.cifti'
  result
}

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
  if((length(dim(mask)) > 3)){
    if(dim(mask)[4] > 1){
      stop('Mask should not contain more than three dimesions')
    }
  }
  dims <- dim(mask)[1:3]
  vals <- sort(unique(as.vector(mask)))
  if(any(!(vals %in% c(0,1)))) warning('Values other than 0/1 or TRUE/FALSE detected in the mask.  Coercing to binary.')
  mask@.Data <- array(1*as.logical(mask), dim=dims)
  V <- sum(mask)
  Q <- dim(GICA)[4]
  if(any(dim(GICA)[1:3] != dims)) stop('First three dimensions of GICA and mask do not match.')
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

  # ESTIMATE MEAN

  if(verbose) cat('\n Estimating Template Mean')
  mean1 <- apply(DR1, c(2,3), mean, na.rm=TRUE)
  mean2 <- apply(DR2, c(2,3), mean, na.rm=TRUE)
  template_mean <- t((mean1 + mean2)/2)

  # ESTIMATE SIGNAL (BETWEEN-SUBJECT) VARIANCE

  # total variance
  if(verbose) cat('\n Estimating Total Variance')
  var_tot1 <- apply(DR1, c(2,3), var, na.rm=TRUE)
  var_tot2 <- apply(DR2, c(2,3), var, na.rm=TRUE)
  var_tot <- t((var_tot1 + var_tot2)/2)

  # noise (within-subject) variance
  if(verbose) cat('\n Estimating Within-Subject Variance')
  DR_diff <- DR1 - DR2;
  var_noise <- t((1/2)*apply(DR_diff, c(2,3), var, na.rm=TRUE))

  # signal (between-subject) variance
  if(verbose) cat('\n Estimating Template (Between-Subject) Variance \n')
  template_var <- var_tot - var_noise
  template_var[template_var < 0] <- 0

  rm(DR1, DR2, mean1, mean2, var_tot1, var_tot2, var_tot, DR_diff)

  if(!is.null(out_fname)){
    out_fname_mean <- paste0(out_fname, '_mean')
    out_fname_var <- paste0(out_fname, '_var')
    out_fname_mask <- paste0(out_fname, '_mask')
    GICA@.Data <- GICA@.Data[,,,inds] #remove non-template ICs
    GICA@dim_[5] <- length(inds)
    template_mean_nifti <- template_var_nifti <- GICA #copy over header information from GICA
    img_tmp <- mask2
    for(l in 1:L){
      img_tmp[mask2==1] <- template_mean[,l]
      template_mean_nifti@.Data[,,,l] <- img_tmp
      img_tmp[mask2==1] <- template_var[,l]
      template_var_nifti@.Data[,,,l] <- img_tmp
    }
    writeNIfTI(template_mean_nifti, out_fname_mean)
    writeNIfTI(template_var_nifti, out_fname_var)
    writeNIfTI(mask2, out_fname_mask)
  }

  result <- list(template_mean=template_mean, template_var=template_var, scale=scale, mask=mask, mask2=mask2, inds=inds)
  class(result) <- 'template.nifti'
  return(result)
}



