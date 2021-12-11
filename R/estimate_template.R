#' Estimate template from DR estimates
#'
#' @param DR1,DR2 the test and retest lists of dual regression estimates
#' @param var_method \code{"unbiased"} (default) or \code{"non-negative"}
#'
#' @return List of two elements: the mean and variance templates
#' @export
#' @importFrom abind abind
estimate_template_from_DR <- function(
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
  # Unbiased: hat(sigmasq_btwn)
  # Non-negative: MSB_div2 === hat(sigmasq_btwn) + hat(sigmasq_noise) / k
  SSB <- 2 * colSums(((DR1 + DR2)/2 - rep(t(template$mean), each=N))^2, na.rm=TRUE)
  MSB_div2 <- t(SSB / (N-1)) / 2
  if (var_method == "unbiased") {
    # Fastest method.
    var_noise <- t( (1/2) * apply(DR1 - DR2, c(2,3), var, na.rm=TRUE) )
    template$var <- MSB_div2 - var_noise/2

    # # Previous, equivalent calculation.
    # var_tot1 <- apply(DR1, c(2,3), var, na.rm=TRUE)
    # var_tot2 <- apply(DR2, c(2,3), var, na.rm=TRUE)
    # var_tot <- t((var_tot1 + var_tot2)/2)
    # # noise (within-subject) variance
    # DR_diff <- DR1 - DR2;
    # var_noise <- t((1/2)*apply(DR_diff, c(2,3), var, na.rm=TRUE))
    # # signal (between-subject) variance
    # template$var <- var_tot - var_noise

    # # Another equivalent calculation.
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
#' @param BOLD Vector of fMRI timeseries data for template estimation. The data
#'  format can be one of the following: a character vector of CIFTI file paths
#'  (*.dtseries.nii), a list of \code{"xifti"} objects, a character vector of
#'  NIFTI file paths (*.nii), a list of \code{"nifti"} objects, or a list of 
#'  \eqn{V \times T} numeric matrices.
#' @param BOLD2 Vector of fMRI timeseries retest data for template estimation. Must
#'  be from the same subjects, in the same order, and in the same format as \code{BOLD}.
#'  If \code{NULL} (default), a pseudo test-retest dataset will be made by splitting
#'  each scan of \code{BOLD} in half. 
#' @param GICA Group ICA maps in a format compatible with \code{BOLD}.
#' @param inds Numeric indices of the group ICs to include in the template. If 
#'  \code{NULL}, use all group ICs (default).
#' 
#'  If \code{inds} is provided, the ICs not included will be removed after calculating
#'  dual regression, not before. This is because removing the ICs prior to dual 
#'  regression would leave unmodeled signals in the data, which could bias the 
#'  templates.
#' @param scale Scale each entry of \code{BOLD} and \code{BOLD2} by its mean spatial
#'  standard deviation before computing dual regression? Default: \code{TRUE}.
#' @param normA Normalize the A matrix (spatial maps)?
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
#'  \code{"nifti"} objects. A binary brain map of the same size as the fMRI data, with
#'  \code{1} corresponding to in-mask voxels.
#' @param var_method Method for estimating the template variance: \code{"unbiased"} 
#'  (default), \code{"non-negative"}, or \code{"both"}. The unbiased template variance is
#'  based on the assumed mixed effects/ANOVA model, whereas the non-negative template
#'  variance adds to it to account for greater potential between-subjects variation.
#'  (The template mean is the same for either choice of \code{var_method}.) 
#' @param keep_DR Keep the DR estimates? If \code{FALSE} (default), do not save the DR 
#'  estimates and only return the templates. If \code{TRUE}, the DR estimates are
#'  returned too. If a single file path, save the DR estimates as an .rds file at
#'  that location. If a list of two vectors of file paths with the same lengths as
#'  \code{BOLD}, save the DR estimates as individual files at these locations in
#'  the appropriate format (CIFTI, NIFTI, or .rds, depending on \code{BOLD}). 
#' @param Q2,maxQ Denoise the dual regression estimates? Denoising is based on modeling and
#'  removing nuisance ICs. It may result in cleaner estimates for smaller datasets, but may
#'  be unnecessary (and time-consuming) for larger datasets. If both are \code{NULL},
#'  denoising will be performed, with the number of nuisance ICs estimated for each fMRI scan
#'  separately. If \code{Q2==0}, do not denoise (default). Otherwise, specify one or the other:
#'  use \code{Q2} to specify the number of nuisance ICs, or \code{maxQ} to specify the number of
#'  total ICs (template/group + nuisance). \eqn{L <= (L+Q2) = maxQ <= T}, where \eqn{L} is the number
#'  of template ICs and \eqn{T} is the number of timepoints in each fMRI scan. 
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
#'
#' @return A list with two entries, \code{"template_mean"} and \code{"template_var"}. There
#'  may be more entries too, depending on the function arguments. 
#'
#' @export
#'
estimate_template <- function(
  BOLD, BOLD2=NULL, 
  GICA, inds=NULL,
  scale=TRUE, normA=FALSE,
  brainstructures=c("left","right"),
  mask=NULL,
  var_method=c("unbiased", "non-negative"),
  keep_DR=FALSE,
  Q2=0, maxQ=NULL,
  out_fname=NULL,
  FC = FALSE, 
  verbose=TRUE) {

  # Check arguments ------------------------------------------
  retest <- !is.null(BOLD2)
  # Determine format of `BOLD`: CIFTI, xifti, NIFTI, nifti, or data
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
  nS <- length(BOLD)

  # Ensure `BOLD2` is the same length.
  if (retest) {
    if (length(BOLD) != length(BOLD2)) {
      stop("`BOLD2` represents corresponding retest data for `BOLD`, so it must have the same length as `BOLD`.")
    }
  }

  # For CIFTI or NIFTI input, check that the file paths exist
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

  # `mask`
  if (FORMAT == "NIFTI") {
    if (is.null(mask)) { stop("`mask` is required.") }
    if (is.character(mask)) { mask <- oro.nifti::readNIfTI(mask) }
    if (length(dim(mask)) == 4 && dim(mask)[4] == 1) { mask <- mask[,,,1] }
    if (length(dim(mask)) != 3) { stop("`mask` should be a 3D binary image.") }
    if (!is.logical(mask)) {
      if (verbose) { cat("Coercing `mask` to a logical array with `as.logical`.\n)" }
      mask[] <- as.logical(mask)
    }
    if (!all(dim(mask) != nV)) {
      stop( 
        "The group ICA images have dimensions", paste(nV, collapse="x") 
        " but the mask has dimensions", paste(dim(mask, collapse="x")), "."
      )
    }
  }

  # `GICA` and `inds`
  if (FORMAT == "CIFTI") {
    if (is.character(GICA)) {
      if (length(GICA) > 1) { warning("Using the first entry of `GICA`."); GICA <- GICA[1] }
      GICA <- read_xifti(GICA, brainstructures=brainstructures)
    }
    stopifnot(is.xifti(GICA))
  } else if (FORMAT == "NIFTI") {
    if (is.character(GICA)) {
      if (length(GICA) > 1) { warning("Using the first entry of `GICA`."); GICA <- GICA[1] }
      GICA <- readNIfTI(GICA)
    }
    stopifnot(length(dim(GICA)) == 4)
    # stopifnot(inherits(GICA, "nifti")) # could be a 4D matrix
  } else {
    stopifnot(length(dim(GICA)) == 2)
  }
  nQ <- dim(GICA)[length(dim(GICA))]
  nV <- dim(GICA)[c(1, length(dim(GICA))-1)]
  GICA_mat <- switch(FORMAT,
    CIFTI = as.matrix(GICA),
    NIFTI = matrix(GICA[rep(mask, nQ)], ncol=nQ),
    DATA = GICA
  )
  GICA_mat <- scale(GICA_mat, scale=FALSE)

  # `keep_DR`
  if (is.logical(keep_DR)) {
    if (length(keep_DR) > 1) { warning("Using the first entry of `keep_DR`."); keep_DR <- keep_DR[1] }
  } else {
    if (is.character(keep_DR)) {
      if (length(keep_DR) > 1) { warning("Using the first entry of `keep_DR`."); keep_DR <- keep_DR[1] }
      if (!endsWith(keep_DR, ".rds")) { keep_DR <- paste0(keep_DR, ".rds") }
      if (!dir.exists(dirname(keep_DR))) { stop('Directory part of `keep_DR` does not exist.') }
    } else if (is.list(keep_DR)) {
      if (length(keep_DR) != 2) {
        stop("If `keep_DR` is a list it must have two entries, each being a vector of file paths the same length as `BOLD`.")
      }
      if (length(keep_DR[[1]]) != nS || length(keep_DR[[2]]) != nS) {
        stop("If `keep_DR` is a list it must have two entries, each being a vector of file paths the same length as `BOLD`.")
      }
      if (!all(dir.exists(dirname(do.call(c, keep_DR))))) { stop('At least one directory part of `keep_DR` does not exist.') }
    }
  }

  # `inds`.
  if (!is.null(inds)) {
    if (!all(inds %in% seq(nQ))) stop('Invalid entries in inds argument.')
    nL <- length(inds)
  } else {
    inds <- seq(nQ)
    nL <- nQ
  }

  # `out_fname`
  if (!is.null(out_fname)) {
    out_fname <- as.character(out_fname)
    if (!dir.exists(dirname(out_fname))) { stop('Directory part of `out_fname` does not exist.') }
    if (!endsWith(out_fname, FORMAT_extn)) { out_fname <- paste0(out_fname, FORMAT_extn) }
    if (length(out_fname) == 1) {
      out_fname <- c(
        gsub(FORMAT_extn, "_mean.dscalar.nii", out_fname),
        gsub(FORMAT_extn, "_var.dscalar.nii", out_fname)
        gsub(FORMAT_extn, "_var_nn.dscalar.nii", out_fname)
      )
      if (var_method != "both") { out_fname <- out_fname[seq(2)] }
    } else if (var_method != "both") {
      if (length(out_fname) > 2) { warning("Using the first two entries of `out_fname`."); out_fname <- out_fname[seq(2)] }
    } else {
      if (length(out_fname) == 2) {
        out_fname <- c(out_fname, gsub(FORMAT_extn, "_var_nn.dscalar.nii", out_fname))
      } else if (length(out_fname) > 3) { warning("Using the first three entries of `out_fname`."); out_fname <- out_fname[seq(3)] }
    }
  }

  # `brainstructures` is handled when needed
  # `Q2` and `maxQ` handled later

  # Rest of arguments
  scale <- as.logical(scale, "scale")
  var_method <- match.arg(var_method, c("unbiased", "non-negative")) 
  normA <- as.logical(normA, "normA")
  verbose <- as.logical(verbose, "verbose")
  FC <- as.logical(FC, "FC")

  retest <- !is.null(BOLD2)
  if(retest){
    if(length(BOLD) != length(BOLD2)) stop('If provided, BOLD2 must have same length as BOLD and be in the same subject order.')
  }

  # Process each scan ---------------------
  if (verbose) {
    nDV <- ifelse(FORMAT=="NIFTI", sum(mask), nV)
    cat('Number of data locations: ', nDV, "\n")
    cat('Number of original group ICs: ', nQ, "\n")
    cat('Number of template ICs: ', nL, "\n")
    cat('Number of training subjects: ', nS, "\n")
  }

  DR1 <- DR2 <- array(NA, dim=c(nS, nL, nV))
  if (FC) { FC1 <- FC2 <- array(NA, dim=c(nS, nL, nL)) }

  for (ii in seq(nS)) {
    if(verbose) cat('\n Reading and analyzing data for subject ', ii,' of ', nS, ".\n")
    if (retest) { B2 <- BOLD2[ii] } else { B2 <- NULL }

    DR_ii <- dual_reg2(
      BOLD[ii], BOLD2=B2, GICA=GICA_mat,
      scale=scale, normA=normA,
      format=format, dim_expect=nDV,
      brainstructures=brainstructures, mask=mask,
      Q2=Q2, maxQ=maxQ, verbose=verbose
    )

    # Shouldn't happen, since missing files are checked for at start of function,
    # unless they were deleted since then. 
    if (!isFALSE(DR_ii$missing)) { next }

    DR1[ii,,] <- DR_ii$test[inds,]
    DR2[ii,,] <- DR_ii$retest[inds,]

    if(FC) {
      FC1[ii,,] <- cov(DR_ii$test[inds,])
      FC2[ii,,] <- cov(DR_ii$retest[inds,])
    }
  }

  # FC template estimation

  # Format and save template
  if (FORMAT == "CIFTI") {
    # Format template as "xifti"s
    GICA <- select_xifti(GICA, inds)
    GICA$meta$cifti$names <- paste0("IC ", inds)
    template$mean <- newdata_xifti(GICA, template$mean)
    template$mean$meta$cifti$misc <- list(template="mean")
    template$var <- newdata_xifti(GICA, template$var)
    template$var$meta$cifti$misc <- list(template="var")

  } else if (FORMAT == "NIFTI") {
    # Format template as "nifti"s
    # GICA@.Data <- GICA@.Data[,,,inds]
    # GICA@dim_[5] <- length(inds) 
    # img_tmp <- mask2
    # for(l in 1:L){
    #   img_tmp[mask2==1] <- template_mean[,l]
    #   template_mean_nifti@.Data[,,,l] <- img_tmp
    #   img_tmp[mask2==1] <- template_var[,l]
    #   template_var_nifti@.Data[,,,l] <- img_tmp
    # }
    GICA <- GICA[,,,inds,drop=FALSE]
    nii_temp <- GICA
    nii_temp[] <- template$mean
    template$mean <- nii_temp
    nii_temp[] <- template_var
    template$var <- nii_temp

  }

  if(!is.null(out_fname)){
    write_cifti(template$mean, out_fname[1], verbose=verbose)
    write_cifti(template$var, out_fname[2], verbose=verbose)
    if (FORMAT == "NIFTI") { writeNIfTI(mask2, 'mask2') }
  }
  # [TO DO] var both

  # Results list.
  result <- list(
    template_mean=xifti_mean, template_var=xifti_var, template_FC=template_FC,
    scale=scale, inds=inds, var_method=var_method
  )
  
  # Keep DR?
  if (is.character(keep_DR) {
    saveRDS(list(DR1=DR1, DR2=DR2), keep_DR)
    keep_DR <- FALSE # no longer need it.
  }
  if (keep_DR) { 
    result <- c(result, list(DR=list(DR1=DR1, DR2=DR2)))
  }

  # Return results.
  class(result) <- paste0(template, "_", tolower(FORMAT))
  result
}