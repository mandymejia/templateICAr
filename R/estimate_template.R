#' Estimate template from DR estimates
#'
#' @param DR, the test and retest dual regression estimates (\eqn{M \times N \times (L \times V)})
#' @param var_method \code{"unbiased"} (default) or \code{"non-negative"}
#' @param LV A length-two integer vector giving the dimensions L and V to reshape the result. 
#'  Default: \code{NULL} (do not reshape the result).
#'
#' @return List of two elements: the mean and variance templates
#' @export
estimate_template_from_DR <- function(
  DR, var_method=c("unbiased", "non-negative"), LV=NULL){

  # Check arguments.
  stopifnot(length(dim(DR)) == 4)
  nM <- dim(DR)[1]
  nN <- dim(DR)[2]
  nLV <- dim(DR)[3]
  var_method <- match.arg(var_method, c("unbiased", "non-negative"))
  if (!is.null(LV)) {
    stopifnot(is.numeric(nLV) && all(nLV > 0) && all(nLV == round(nLV)))
    stopifnot(prod(LV) == nLV)
  }

  vd <- var_decomp(DR)
  template <- list(
    mean = vd$grand_mean,
    var = switch(var_method,
      unbiased = (vd$SSB / (nN-1)) / 2,
      `non-negative` = ((vd$SSB / (nN-1) - vd$SSR / ((nM-1)*(nN-1)))) / 2
    )
  )

  if (!is.null(LV)) {
    template <- lapply(template, function(x){ matrix(x, nrow=LV[1], ncol=LV[2]) })
  }

  template
}

#' Estimate template from DR estimates (when there are two measurements)
#'
#' @param DR1,DR2 the test and retest dual regression estimates (\eqn{N \times L \times V})
#' @param var_method \code{"unbiased"} (default) or \code{"non-negative"}
#'
#' @return List of two elements: the mean and variance templates
#' @export
estimate_template_from_DR_two <- function(
  DR1, DR2, var_method=c("unbiased", "non-negative")){

  # Check arguments.
  stopifnot(length(dim(DR1)) == length(dim(DR2)))
  stopifnot(all(dim(DR1) == dim(DR2)))
  N <- dim(DR1)[1]
  var_method <- match.arg(var_method, c("unbiased", "non-negative"))

  template <- list(mean=NULL, var=NULL)

  # Mean.
  template$mean <- t(colMeans(DR1 + DR2, na.rm=TRUE) / 2)

  # Variance.
  SSB <- 2 * colSums(((DR1 + DR2)/2 - rep(t(template$mean), each=N))^2, na.rm=TRUE)
  MSB_div2 <- t(SSB / (N-1)) / 2
  if (var_method == "unbiased") {
    # 1. Fastest method.
    var_noise <- t( (1/2) * apply(DR1 - DR2, c(2,3), var, na.rm=TRUE) )
    template$var <- MSB_div2 - var_noise/2

    # # 2. Previous, equivalent calculation.
    # var_tot1 <- apply(DR1, c(2,3), var, na.rm=TRUE)
    # var_tot2 <- apply(DR2, c(2,3), var, na.rm=TRUE)
    # var_tot <- t((var_tot1 + var_tot2)/2)
    # # noise (within-subject) variance
    # DR_diff <- DR1 - DR2;
    # var_noise <- t((1/2)*apply(DR_diff, c(2,3), var, na.rm=TRUE))
    # # signal (between-subject) variance
    # template$var <- var_tot - var_noise
    # 
    # # 3. Another equivalent calculation.
    # template$var <- t(apply(
    #   abind::abind(DR1, DR2, along=1),
    #   seq(2, 3),
    #   function(q){ cov(q[seq(N)], q[seq(N+1, 2*N)], use="complete.obs") }
    # ))

    # Make negative estimates equal to zero.
    template$var[template$var < 0] <- 0

  } else {
    template$var <- MSB_div2
  }

  template
}

#' Estimate template
#'
#' Estimate template for Template or Diagnostic ICA based on fMRI data
#' 
#' All fMRI data (entries in \code{BOLD} and \code{BOLD2}, and \code{GICA}) must be in
#'  the same spatial resolution.
#'
#' @param BOLD,BOLD2 Vector of subject-level fMRI data in one of the following formats: 
#'  CIFTI file paths, \code{"xifti"} objects, NIFTI file paths, \code{"nifti"} objects, or
#'  \eqn{V \times T} numeric matrices, where \eqn{V} is the number of data locations and
#'  \eqn{T} is the number of timepoints.
#' 
#'  If \code{BOLD2} is provided it must be in the same format as \code{BOLD}; 
#'  \code{BOLD} will be the test data and \code{BOLD2} will be the retest data. 
#'  \code{BOLD2} should be the same length as \code{BOLD} and have the same subjects in the same order.
#'  If \code{BOLD2} is not provided, \code{BOLD} will be split in half; 
#'  the first half will be the test data and the second half will be the retest data.
#' @param GICA Group ICA maps in a format compatible with \code{BOLD}. Can also be a
#'  (vectorized) numeric matrix (\eqn{V \times Q}) no matter the format of \code{BOLD}.
#' @param inds Numeric indices of the group ICs to include in the template. If 
#'  \code{NULL}, use all group ICs (default).
#' 
#'  If \code{inds} is provided, the ICs not included will be removed after calculating
#'  dual regression, not before. This is because removing the ICs prior to dual 
#'  regression would leave unmodeled signals in the data, which could bias the 
#'  templates.
#' @param center_rows,center_cols Center BOLD data across rows (each data location's time series) or columns (each time point's image)? Default: \code{TRUE} for both.
#' @param center_Gcols Center GICA across columns (each ICA)? Default: \code{TRUE}.
#' @param scale A logical value indicating whether the fMRI timeseries should be scaled by the image standard deviation.
#' @param detrend_DCT Detrend the data? This is the number of DCT bases to use for detrending. If \code{0} (default), do not detrend.
#' @param normA Scale each IC timeseries (column of \eqn{A}) in the dual regression 
#'  estimates? Default: \code{FALSE}. (The opposite scaling will be applied to \eqn{S}
#'  such that the product \eqn{A \times S} remains the same).
#' @param brainstructures Only applies if the entries of \code{BOLD} are CIFTI file paths. 
#'  Character vector indicating which brain structure(s)
#'  to obtain: \code{"left"} (left cortical surface), \code{"right"} (right
#'  cortical surface) and/or \code{"subcortical"} (subcortical and cerebellar
#'  gray matter). Can also be \code{"all"} (obtain all three brain structures).
#'  Default: \code{c("left","right")} (cortical surface only).
#' @param mask Required if and only if the entries of \code{BOLD} are NIFTI file paths or
#'  \code{"nifti"} objects. This is a brain map formatted as a binary array of the same 
#'  size as the fMRI data, with \code{TRUE} corresponding to in-mask voxels.
#' @param var_method Method for estimating the template variance: \code{"unbiased"} 
#'  (default), \code{"non-negative"}, or \code{"both"}. The unbiased template variance is
#'  based on the assumed mixed effects/ANOVA model, whereas the non-negative template
#'  variance adds to it to account for greater potential between-subjects variation.
#'  (The template mean is the same for either choice of \code{var_method}.) 
#' @param keep_DR Keep the DR estimates? If \code{FALSE} (default), do not save the DR 
#'  estimates and only return the templates. If \code{TRUE}, the DR estimates are
#'  returned too. If a single file path, save the DR estimates as an RDS file at
#'  that location. If a list of two vectors of file paths with the same lengths as
#'  \code{BOLD}, save the DR estimates as individual files at these locations in
#'  the appropriate format (CIFTI, NIFTI, or RDS files, depending on \code{BOLD}). 
#' @param Q2,maxQ Obtain dual regression estimates after denoising? Denoising is based on modeling and
#'  removing nuisance ICs. It may result in a cleaner estimate for smaller datasets, but it may be unnecessary (and time-consuming) for larger datasets.
#'  If both arguments are \code{NULL}, denoising will be performed, with the number of nuisance 
#'  ICs estimated for \code{BOLD} and \code{BOLD2} separately. Otherwise, specify one or the other:
#'  use \code{Q2} to specify the number of nuisance ICs, or \code{maxQ} to specify the number of
#'  total ICs (group + nuisance, or \eqn{Q + Q2}). Set either to zero to skip denoising.
#'  Default: \code{Q2==0} (do not denoise).
#'  
#'  The valid inputs are \eqn{Q <= (Q+Q2) = maxQ <= T}, where \eqn{Q} is the number
#'  of group ICs and \eqn{T} is the number of timepoints in each fMRI scan. 
#' @param out_fname Character vector of file path(s) to write the mean and variance templates to.
#'  If one file name is provided, it will be appended with \code{"_mean.dscalar.nii"} for the
#'  template mean map and \code{"_var.dscalar.nii"} for the template variance map. If two
#'  file names are provided, the first will be used for the template mean and the second will
#'  be used for the template variance. If \code{var_method=="both"}, the non-negative template 
#'  variance will be appended with \code{"_var_nn.dscalar.nii"}.
#' @param FC Include the functional connectivity template?
#' @param verbose Display progress updates? Default: \code{TRUE}.
#'
#' @importFrom stats cov quantile
#' @importFrom ciftiTools read_xifti is.xifti write_cifti
#'
#' @return A list with two entries, \code{"template_mean"} and \code{"template_var"}. There
#'  may be more entries too, depending on the function arguments. 
#'
#' @export
#'
estimate_template <- function(
  BOLD, BOLD2=NULL, 
  GICA, inds=NULL,
  center_rows=TRUE, center_cols=TRUE, scale=TRUE, detrend_DCT=0, 
  center_Gcols=TRUE, normA=FALSE,
  Q2=0, maxQ=NULL,
  brainstructures=c("left","right"), mask=NULL,
  var_method=c("unbiased", "non-negative"),
  keep_DR=FALSE,
  out_fname=NULL,
  FC=FALSE, 
  verbose=TRUE) {

  # Check arguments ------------------------------------------------------------
  
  # Simple argument checks.
  stopifnot(is.logical(center_rows) && length(center_rows)==1)
  stopifnot(is.logical(center_cols) && length(center_cols)==1)
  stopifnot(is.logical(scale) && length(scale)==1)
  stopifnot(is.numeric(detrend_DCT) && length(detrend_DCT)==1)
  stopifnot(detrend_DCT >=0 && detrend_DCT==round(detrend_DCT))
  stopifnot(is.logical(normA) && length(normA)==1)
  var_method <- match.arg(var_method, c("unbiased", "non-negative")) 
  if (!is.null(Q2) && !is.null(maxQ)) { stop("Specify one of `Q2` or `maxQ`.") }
  stopifnot(is.logical(FC) && length(FC)==1)
  stopifnot(is.logical(verbose) && length(verbose)==1)
  retest <- !is.null(BOLD2)

  # `keep_DR`
  if (is.logical(keep_DR)) {
    stopifnot(length(keep_DR)==1)
  } else {
    if (is.character(keep_DR)) {
      stopifnot(length(keep_DR)==1)
      if (!endsWith(keep_DR, ".rds")) { keep_DR <- paste0(keep_DR, ".rds") }
      if (!dir.exists(dirname(keep_DR))) { stop('Directory part of `keep_DR` does not exist.') }
    } else if (is.list(keep_DR)) {
      if (length(keep_DR) != 2) {
        stop("If `keep_DR` is a list it must have two entries, each being a vector of file paths the same length as `BOLD`.")
      }
      if (length(keep_DR[[1]]) != nN || length(keep_DR[[2]]) != nN) {
        stop("If `keep_DR` is a list it must have two entries, each being a vector of file paths the same length as `BOLD`.")
      }
      if (!all(dir.exists(dirname(do.call(c, keep_DR))))) { stop('At least one directory part of `keep_DR` does not exist.') }
    }
  }

  # `out_fname`
  if (!is.null(out_fname)) {
    out_fname <- as.character(out_fname)
    if (!dir.exists(dirname(out_fname))) { stop('Directory part of `out_fname` does not exist.') }
    if (!endsWith(out_fname, FORMAT_extn)) { out_fname <- paste0(out_fname, FORMAT_extn) }
    if (length(out_fname) == 1) {
      out_fname <- c(
        gsub(FORMAT_extn, paste0("_mean", FORMAT_extn), out_fname),
        gsub(FORMAT_extn, paste0("_var", FORMAT_extn), out_fname),
        gsub(FORMAT_extn, paste0("_varNN", FORMAT_extn), out_fname)
      )
      if (var_method != "both") { out_fname <- out_fname[seq(2)] }
    } else if (var_method != "both") {
      if (length(out_fname) > 2) { warning("Using the first two entries of `out_fname`."); out_fname <- out_fname[seq(2)] }
    } else {
      if (length(out_fname) == 2) {
        out_fname <- c(out_fname, gsub(FORMAT_extn, paste0("_varNN", FORMAT_extn), out_fname))
      } else if (length(out_fname) > 3) { warning("Using the first three entries of `out_fname`."); out_fname <- out_fname[seq(3)] }
    }
  }

  # `xii1` will be used to format output
  xii1 <- NULL

  # `BOLD` and `BOLD2` ---------------------------------------------------------
  # Determine the format of `BOLD` and `BOLD2`. 
  format <- infer_BOLD_format(BOLD)
  if (retest) {
    format2 <- infer_BOLD_format(BOLD2)
    if (format2 != format) {
      stop("`BOLD` format is ", format, ", but `BOLD2` format is ", format2, ".")
    }
  }
  FORMAT <- switch(format,
    CIFTI = "CIFTI",
    xifti = "CIFTI",
    NIFTI = "NIFTI",
    nifti = "NIFTI",
    data = "DATA"
  )
  FORMAT_extn <- switch(FORMAT, CIFTI=".dscalar.nii", NIFTI=".nii", DATA=".rds")
  nN <- length(BOLD)

  # Ensure `BOLD2` is the same length.
  if (retest) {
    if (length(BOLD) != length(BOLD2)) {
      stop("`BOLD2` represents corresponding retest data for `BOLD`, so it must have the same length as `BOLD`.")
    }
  }

  # If BOLD (and BOLD2) is a CIFTI or NIFTI file, check that the file paths exist.
  if (format %in% c("CIFTI", "NIFTI")) {
    missing_BOLD <- (!file.exists(BOLD)) | (!file.exists(BOLD2))
    if (all(missing_BOLD)) stop('The files in `BOLD` do not exist.')
    if (retest) {
      missing_BOLD2 <- !file.exists(BOLD2)
      if (all(missing_BOLD2)) stop('The files in `BOLD2` do not exist.')
      # determine pairs with at least one missing scan
      missing_BOLD <- missing_BOLD | missing_BOLD2
      if (all(missing_BOLD)) stop('Files in `BOLD` and/or `BOLD2` are missing such that no complete pair of data exists.')
    }
    if (any(missing_BOLD)) {
      if (retest) {
        warning('There are ', missing_BOLD, ' pairs of `BOLD` and `BOLD2` with at least one non-existent scan. These pairs will be excluded from template estimation.')
        BOLD <- BOLD[!missing_BOLD]
        BOLD2 <- BOLD[!missing_BOLD]
      } else {
        warning('There are ', missing_BOLD, ' scans in `BOLD` that do not exist. These scans will be excluded from template estimation.')
        BOLD <- BOLD[!missing_BOLD]
      }
    }
  }

  # `GICA` ---------------------------------------------------------------------
  # Conver `GICA` to a numeric data matrix or array.
  if (FORMAT == "CIFTI") {
    if (is.character(GICA)) { GICA <- ciftiTools::read_xifti(GICA, brainstructures=brainstructures) }
    if (is.xifti(GICA)) { 
      xii1 <- select_xifti(GICA, 1) # for formatting output
      GICA <- as.matrix(GICA)
    }
    stopifnot(is.matrix(GICA))
  } else if (FORMAT == "NIFTI") {
    if (is.character(GICA)) { GICA <- oro.nifti::readNIfTI(GICA, reorient=FALSE) }
    stopifnot(length(dim(GICA)) > 1)
  } else {
    stopifnot(is.matrix(GICA))
  }
  nQ <- dim(GICA)[length(dim(GICA))]

  # `inds`.
  if (!is.null(inds)) {
    if (!all(inds %in% seq(nQ))) stop('Invalid entries in inds argument.')
    nL <- length(inds)
  } else {
    inds <- seq(nQ)
    nL <- nQ
  }

  # `mask` ---------------------------------------------------------------------
  # Get `mask` as a logical array.
  # Check `GICA` and `mask` dimensions match.
  # Vectorize `GICA`.
  if (FORMAT == "NIFTI") {
    if (is.null(mask)) { stop("`mask` is required.") }
    if (is.character(mask)) { mask <- oro.nifti::readNIfTI(mask, reorient=FALSE) }
    if (dim(mask)[length(dim(mask))] == 1) { mask <- array(mask, dim=dim(mask)[length(dim(mask))-1]) }
    if (is.numeric(mask)) { 
      cat("Coercing `mask` to a logical array with `as.logical`.\n")
      mask[] <- as.logical(mask)
    }
    nI <- dim(mask); nV <- sum(mask)
    stopifnot(length(dim(GICA)) %in% c(2, length(nI)+1))
    if (length(dim(GICA)) == length(nI)+1) {
      if (length(dim(GICA)) != 2) {
        stopifnot(all(dim(GICA)[length(dim(GICA))-1] == nI))
      }
      if (all(dim(GICA)[length(dim(GICA))-1] == nI)) {
        GICA <- matrix(GICA[rep(mask, nQ)], ncol=nQ)
        stopifnot(nrow(GICA) == nV)  
      }
    }
  } else {
    nI <- nV <- nrow(GICA)
  }

  # Center `GICA` columns.
  if (center_Gcols) { GICA - rep(colMeans(GICA), rep.int(nV, nQ)) }

  # Process each scan ----------------------------------------------------------
  if (verbose) {
    cat('Number of data locations:      ', nV, "\n")
    cat('Number of original group ICs:  ', nQ, "\n")
    cat('Number of template ICs:        ', nL, "\n")
    cat('Number of training subjects:   ', nN, "\n")
  }

  nM <- 2
  DR0 <- array(NA, dim=c(nM, nN, nL, nV)) # measurements by subjects by components by locations
  if(FC) FC0 <- array(NA, dim=c(nM, nN, nL, nL)) #for functional connectivity template

  for (ii in seq(nN)) {
    if(verbose) cat('\n Reading and analyzing data for subject ', ii,' of ', nN, ".\n")
    if (retest) { B2 <- BOLD2[[ii]] } else { B2 <- NULL }

    DR_ii <- dual_reg2(
      BOLD[ii], BOLD2=B2, 
      format=format,       
      GICA=GICA,
      center_rows=center_rows, center_cols=center_cols, 
      scale=scale, detrend_DCT=detrend_DCT,
      normA=normA,
      Q2=Q2, maxQ=maxQ,
      brainstructures=brainstructures, mask=mask,
      verbose=verbose
    )

    DR0[1,ii,,] <- DR_ii$test[inds,]
    DR0[2,ii,,] <- DR_ii$retest[inds,]
    if(FC) {
      FC0[1,ii,,] <- cov(DR_ii$test[,inds])
      FC0[2,ii,,] <- cov(DR_ii$retest[,inds])
    }
  }
  mask2 <- NULL # [TO DO]: fix

  # Vectorize components and locations
  DR0 <- array(DR0, dim=c(nM, nN, nL*nV))
  # FC0 <- array(FC1, dim=c(nM, nN, nL*nL))

  # Estimate mean and var template
  if (verbose) { cat("Estimating template.\n") }
  template <- estimate_template_from_DR(DR0, c(nL, nV))

  # Keep DR
  if (!isFALSE(keep_DR)) {
    if (is.character(keep_DR)) {
      if (length(keep_DR) > 1) {
        warning("Using first entry of `keep_DR`.")
        keep_DR <- keep_DR[1]
      }
      if (!endsWith(keep_DR, ".rds")) { keep_DR <- paste0(keep_DR, ".rds") }
      DR0 <- array(DR0, dim=c(nM, nN, nL, nV)) # Undo vectorize
      saveRDS(DR0, keep_DR)
      keep_DR <- FALSE # no longer need it.
    } else if (!isTRUE(keep_DR)) {
      warning("`keep_DR` should be `TRUE`, `FALSE`, or a file path. Using `FALSE`.")
      keep_DR <- FALSE
    }
  }
  if (!keep_DR) { rm(DR0) }

  # Estimate FC template
  if(FC){

    mean_FC <- var_FC_tot <- var_FC_within <- NULL
    # mean_FC <- (apply(FC1, c(2,3), mean, na.rm=TRUE) + apply(FC2, c(2,3), mean, na.rm=TRUE))/2
    # var_FC_tot  <- (apply(FC1, c(2,3), var, na.rm=TRUE) + apply(FC2, c(2,3), var, na.rm=TRUE))/2
    # var_FC_within  <- 1/2*(apply(FC1-FC2, c(2,3), var, na.rm=TRUE))
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

    nu_est <- matrix(NA, nL, nL)
    for(q1 in 1:nL){
      for(q2 in q1:nL){

        #estimate k = nu - p - 1
        nu_opt <- optimize(f=fun, interval=c(nL+1,nL*10), p=nL, var_ij=var_FC_between[q1,q2], xbar_ij=mean_FC[q1,q2], xbar_ii=mean_FC[q1,q1], xbar_jj=mean_FC[q2,q2])
        nu_est[q1,q2] <- nu_opt$minimum
      }
    }
    nu_est[lower.tri(nu_est, diag=FALSE)] <- NA
    nu_est1 <- quantile(nu_est[upper.tri(nu_est, diag=TRUE)], 0.1, na.rm = TRUE)

    template_FC <- list(nu = nu_est1,
                        psi = mean_FC*(nu_est1 - nL - 1))
  } else {
    template_FC <- NULL
  }

  # Format and save template
  if (FORMAT == "CIFTI" && !is.null(xii1)) {
    # Format template as "xifti"s
    GICA <- newdata_xifti(select_xifti(xii1, rep(1, nL)), GICA)
    GICA$meta$cifti$names <- paste0("IC ", inds)
    template$mean <- newdata_xifti(GICA, template$mean)
    template$mean$meta$cifti$misc <- list(template="mean")
    template$var <- newdata_xifti(GICA, template$var)
    template$var$meta$cifti$misc <- list(template="var")
    if (!is.null(out_fname)) {
      write_cifti(template$mean, out_fname[1], verbose=verbose)
      write_cifti(template$var, out_fname[2], verbose=verbose)
      # if (FORMAT == "NIFTI") { writeNIfTI(mask2, 'mask2') }
    }
  } else if (FORMAT == "NIFTI") {
    # GICA@.Data <- GICA@.Data[,,,inds]
    # GICA@dim_[5] <- length(inds) 
    # img_tmp <- mask2
    # for(l in 1:L){
    #   img_tmp[mask2==1] <- template_mean[,l]
    #   template_mean_nifti@.Data[,,,l] <- img_tmp
    #   img_tmp[mask2==1] <- template_var[,l]
    #   template_var_nifti@.Data[,,,l] <- img_tmp
    # }
    nii_temp <- array(mask, dim=c(dim(mask), nL))
    nii_temp[nii_temp[]] <- GICA
    GICA <- nii_temp
    nii_temp <- array(mask, dim=c(dim(mask), nL))
    nii_temp[nii_temp[]] <- template$mean
    template$mean <- nii_temp
    nii_temp <- array(mask, dim=c(dim(mask), nL))
    nii_temp[nii_temp[]] <- template$var
    template$var <- nii_temp
    if (!is.null(out_fname)) {
      writeNIfTI(template$mean, out_fname[1])
      writeNIfTI(template$var, out_fname[2])
      writeNIfTI(mask2, 'mask2') # [TO DO] fix
    }

  } else {
    if (!is.null(out_fname)) {
      saveRDS(template$mean, out_fname[1])
      saveRDS(template$var, out_fname[2])
      if (FORMAT == "NIFTI") { saveRDS(mask2, 'mask2') }
    }
  }

  # Results list.
  result <- list(
    template_mean=template$mean, 
    template_var=template$var, 
    template_FC=template$FC,
    opts=list(
      center_rows=center_rows, 
      center_cols=center_cols, 
      scale=scale, 
      detrend_DCT=detrend_DCT, 
      center_Gcols=center_Gcols,
      normA=normA,
      inds=inds,
      var_method=var_method
    )
  )

  # Return results.
  class(result) <- paste0(template, "_", tolower(FORMAT))
  result
}

#' @rdname estimate_template
#' 
estimate_template.cifti <- function(
  BOLD, BOLD2=NULL, 
  GICA, inds=NULL,
  center_rows=TRUE, center_cols=TRUE, scale=TRUE, detrend_DCT=0, 
  center_Gcols=TRUE, normA=FALSE,
  brainstructures=c("left","right"), 
  var_method=c("unbiased", "non-negative"),
  keep_DR=FALSE,
  Q2=0, maxQ=NULL,
  out_fname=NULL,
  FC=FALSE, 
  verbose=TRUE) {

  estimate_template(
    BOLD=BOLD, BOLD2=BOLD2,
    GICA=GICA, inds=inds,
    center_rows=center_rows, center_cols=center_cols, scale=scale, detrend_DCT=detrend_DCT, 
    center_Gcols=center_Gcols, normA=normA,
    brainstructures=brainstructures, 
    var_method=var_method,
    keep_DR=keep_DR,
    Q2=Q2, maxQ=maxQ,
    out_fname=out_fname,
    FC=FC, 
    verbose=verbose
  )
}

#' Summarize a \code{"template_cifti"} object
#'
#' Summary method for class \code{"template_cifti"}
#'
#' @param object Object of class \code{"template_cifti"}. 
#' @param ... further arguments passed to or from other methods.
#' @export
#' @method summary template_cifti
summary.template_cifti <- function(object, ...) {

  x <- c(
    summary(object$template_mean),
    list(has_DR="DR" %in% names(object)),
    object[c("scale", "inds", "var_method")]
  )

  class(x) <- "summary.template_cifti"
  return(x)
}

#' @rdname summary.template_cifti
#' @export
#' 
#' @param x The template from \code{estimate_template.cifti}
#' @param ... further arguments passed to or from other methods.
#' @method print summary.template_cifti
print.summary.template_cifti <- function(x, ...) {

  cat("====TEMPLATE INFO====================\n")
  cat("Spatially scaled:", x$scale, "\n")
  cat("Variance method: ", x$var_method, "\n")
  cat("\n")

  class(x) <- "summary.xifti"
  print(x) 
}

#' @rdname summary.template_cifti
#' @export
#' 
#' @method print template_cifti
print.template_cifti <- function(x, ...) {
  print.summary.template_cifti(summary(x))
}

#' Plot template
#' 
#' @param x The template from \code{estimate_template.cifti}
#' @param stat \code{"mean"} (default) or \code{"var"}
#' @param ... Additional arguments to \code{view_xifti}
#' @return The plot
#' @export
#' @importFrom ciftiTools view_xifti
#' @method plot template_cifti
plot.template_cifti <- function(x, stat=c("mean", "var"), ...) {
  stopifnot(inherits(x, "template_cifti"))
  stat <- match.arg(stat, c("mean", "var"))
  view_xifti(x[[paste0("template_", stat)]], ...)
}

#' @rdname estimate_template
#' 
estimate_template.nifti <- function(
  BOLD, BOLD2=NULL, 
  GICA, inds=NULL,
  center_rows=TRUE, center_cols=TRUE, scale=TRUE, detrend_DCT=0, 
  center_Gcols=TRUE, normA=FALSE,
  mask=mask, 
  var_method=c("unbiased", "non-negative"),
  keep_DR=FALSE,
  Q2=0, maxQ=NULL,
  out_fname=NULL,
  FC=FALSE, 
  verbose=TRUE) {

  estimate_template(
    BOLD=BOLD, BOLD2=BOLD2,
    GICA=GICA, inds=inds,
    center_rows=center_rows, center_cols=center_cols, scale=scale, detrend_DCT=detrend_DCT, 
    center_Gcols=center_Gcols, normA=normA,
    mask=mask, 
    var_method=var_method,
    keep_DR=keep_DR,
    Q2=Q2, maxQ=maxQ,
    out_fname=out_fname,
    FC=FC, 
    verbose=verbose
  )
}