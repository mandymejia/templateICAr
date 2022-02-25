#' Estimate variance decomposition and templates from DR estimates
#'
#' @param DR the test/retest(s) dual regression estimates, as an array with
#'  dimensions \eqn{M \times N \times (L \times V)}, where \eqn{M} is the number
#'  of visits (2), \eqn{N} is the number of subjects, \eqn{L} is the number of
#'  IC networks, and \eqn{V} is the number of data locations.
#'
#'  (\eqn{L} and \eqn{V} are collapsed because they are treated equivalently
#'  in the context of the variance decomposition).
#' @param LV A length-two integer vector giving the dimensions \eqn{L} and
#'  \eqn{V} to reshape the result. Default: \code{NULL} (do not reshape the
#'  result).
#'
#' @return List of two elements: the templates and the variance decomposition.
#'
#'  There are two version of the variance template: \code{varUB} gives the
#'  unbiased variance estimate with values clamped to above zero, and
#'  \code{varNN} gives the upwardly-biased non-negative variance estimate.
#' @export
estimate_template_from_DR <- function(
  DR, LV=NULL){

  # Check arguments.
  stopifnot(length(dim(DR)) == 3)
  nM <- dim(DR)[1]  # visits
  nN <- dim(DR)[2]  # subjects
  nLV <- dim(DR)[3] # locations & networks
  if (!is.null(LV)) {
    stopifnot(is.numeric(nLV) && all(nLV > 0) && all(nLV == round(nLV)))
    stopifnot(prod(LV) == nLV)
  }

  # Variance decomposition
  vd <- var_decomp(DR)

  # Template calculation
  # Below true for M==2. Double check correct for M > 3? (Not used currently.)
  MSB_divM <- (vd$SSB / (nN-1)) / nM
  MSE_divM <- (vd$SSR / ((nM-1)*(nN-1))) / nM
  template <- list(
    mean = vd$grand_mean,
    varUB = MSB_divM - MSE_divM,
    varNN = MSB_divM
  )

  # Format `vd`
  vd$nM <- vd$grand_mean <- NULL # Get rid of redundant entries
  names(vd)[names(vd) == "nS"] <- "num_subjects"
  # Add the variance templates in matrix form to `vd`.
  vd$tmean <- template$mean
  vd$tvarUB <- template$varUB
  vd$tvarNN <- template$varNN

  # Format `template`: clamp var est above zero.
  template$varNN <- pmax(0, template$varNN)

  # Format as matrix if applicable.
  if (!is.null(LV)) {
    template <- lapply(template, function(x){ matrix(x, nrow=LV[1], ncol=LV[2]) })
    vd[names(vd)!="num_subjects"] <- lapply(vd[names(vd)!="num_subjects"],
      function(x){ matrix(x, nrow=LV[1], ncol=LV[2]) }
    )
  }

  # Return
  list(template=template, var_decomp=vd)
}

#' Estimate template from DR estimates (when there are two measurements)
#'
#' Legacy version of \code{\link{estimate_template_from_DR}}
#'
#' @param DR1,DR2 the test and retest dual regression estimates (\eqn{N \times L \times V})
#' @param var_method \code{"non-negative"} (default) or \code{"unbiased"}
#'
#' @return List of two elements: the mean and variance templates
#' @keywords internal
estimate_template_from_DR_two <- function(
  DR1, DR2, var_method=c("non-negative", "unbiased")){

  # Check arguments.
  stopifnot(length(dim(DR1)) == length(dim(DR2)))
  stopifnot(all(dim(DR1) == dim(DR2)))
  N <- dim(DR1)[1]
  var_method <- match.arg(var_method, c("non-negative", "unbiased"))

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
#'  Its columns will be centered.
#' @param inds Numeric indices of the group ICs to include in the template. If
#'  \code{NULL}, use all group ICs (default).
#'
#'  If \code{inds} is provided, the ICs not included will be removed after calculating
#'  dual regression, not before. This is because removing the ICs prior to dual
#'  regression would leave unmodelled signals in the data, which could bias the
#'  templates.
#' @param scale \code{"global"} (default), \code{"local"}, or \code{"none"}.
#'  Global scaling will divide the entire data matrix by the image standard 
#'  deviation (\code{sqrt(mean(rowVars(BOLD)))}). Local scaling will divide each
#'  data location's time series by its estimated standard deviation. 
#' @param scale_sm_FWHM Only applies if \code{scale=="local"}. To
#'  smooth the standard deviation estimates used for local scaling, provide the 
#'  smoothing FWHM (default: \code{2}). if \code{0}, do not smooth.
#' @param detrend_DCT Detrend the data? This is an integer number of DCT bases
#'  to use for detrending. If \code{0} (default), do not detrend.
#' @param center_Bcols Center BOLD across columns (each image)? Default: \code{FALSE}
#'  (not recommended).
#' @param normA Scale each IC timeseries (column of \eqn{A}) in the dual
#'  regression estimates? Default: \code{FALSE} (not recommended). Note that the
#'  product \eqn{A \times S} remains the same with either option.
#' @param brainstructures Only applies if the entries of \code{BOLD} are CIFTI file paths.
#'  Character vector indicating which brain structure(s)
#'  to obtain: \code{"left"} (left cortical surface), \code{"right"} (right
#'  cortical surface) and/or \code{"subcortical"} (subcortical and cerebellar
#'  gray matter). Can also be \code{"all"} (obtain all three brain structures).
#'  Default: \code{c("left","right")} (cortical surface only).
#' @param mask Required if and only if the entries of \code{BOLD} are NIFTI file paths or
#'  \code{"nifti"} objects. This is a brain map formatted as a binary array of the same
#'  size as the fMRI data, with \code{TRUE} corresponding to in-mask voxels.
#' @param var_method Method for estimating the template variance: \code{"non-negative"}
#'  (default) or \code{"unbiased"}. The unbiased template variance is
#'  based on the assumed mixed effects/ANOVA model, whereas the non-negative template
#'  variance adds to it to account for greater potential between-subjects variation.
#'  (The template mean is the same for either choice of \code{var_method}.)
#'
#'  In either case, both variance template version will be returned in matrix
#'  form in the \code{var_decomp} entry of the output.
#' @param keep_DR Keep the DR estimates? If \code{FALSE} (default), do not save the DR
#'  estimates and only return the templates. If \code{TRUE}, the DR estimates are
#'  returned too. If a single file path, save the DR estimates as an RDS file at
#'  that location rather than returning them.
#   [TO DO] If a list of two vectors of file paths with the same lengths as
#   \code{BOLD}, save the DR estimates as individual files at these locations in
#   the appropriate format (CIFTI, NIFTI, or RDS files, depending on \code{BOLD}).
#' @param Q2,Q2_max Obtain dual regression estimates after denoising? Denoising is
#'  based on modeling and removing nuisance ICs. It may result in a cleaner
#'  estimate for smaller datasets, but it may be unnecessary (and time-consuming)
#'  for larger datasets.
#'
#'  Set \code{Q2} to control denoising: use a positive integer to specify the
#'  number of nuisance ICs, \code{NULL} to have the number of nuisance ICs
#'  estimated by PESEL, or zero (default) to skip denoising.
#'
#'  If \code{is.null(Q2)}, use \code{Q2_max} to specify the maximum number of
#'  nuisance ICs that should be estimated by PESEL. \code{Q2_max} must be less
#'  than \eqn{T * .75 - Q} where \eqn{T} is the minimum number of timepoints in
#'  each fMRI scan and \eqn{Q} is the number of group ICs. If \code{NULL}
#'  (default), \code{Q2_max} will be set to \eqn{T * .50 - Q}, rounded.
#' @param out_fname Length-3 character vector of file path(s) to save the output to:
#'  the mean template, the variance template, and the variance decomposition,
#'  in that order. If one file name is provided, it will be appended with
#'  \code{"_mean.[file_ext]"} for the template mean map,
#'  \code{"_var.[file_ext]"} for the template variance map, and
#'  \code{"_varDecomp.rds"} for the variance decomposition, where \code{[file_ext]}
#'  will be \code{"dscalar.nii"} for CIFTI input, \code{"nii"} for NIFTI input,
#'  and \code{"rds"} for data input.
#' @param FC Include the functional connectivity template? Default: \code{FALSE}
#'  (work in progress, not available yet).
#' @param varTol Tolerance for variance of each data location. For each scan,
#'  locations which do not meet this threshold are masked out of the analysis.
#'  Default: \code{1e-6}. Variance is calculated on the original data, before
#'  any normalization.
#' @param maskTol For computing the dual regression results for each subject:
#'  tolerance for number of locations masked out due to low
#'  variance or missing values. If more than this many locations are masked out,
#'  a subject is skipped without calculating dual regression. \code{maskTol}
#'  can be specified either as a proportion of the number of locations (between
#'  zero and one), or as a number of locations (integers greater than one).
#'  Default: \code{.1}, i.e. up to 10 percent of locations can be masked out.
#'
#'  If \code{BOLD2} is provided, masks are calculated for both scans and then
#'  the intersection of the masks is used, for each subject.
#' @param missingTol For computing the variance decomposition across all subjects:
#'  tolerance for number of subjects masked out due to low variance or missing
#'  values at a given location. If more than this many subjects are masked out,
#'  the location's value will be \code{NA} in the templates. \code{missingTol}
#'  can be specified either as a proportion of the number of locations (between
#'  zero and one), or as a number of locations (integers greater than one).
#'  Default: \code{.1}, i.e. up to 10 percent of subjects can be masked out.
#' @param verbose Display progress updates? Default: \code{TRUE}.
#'
#' @importFrom stats cov quantile
#' @importFrom ciftiTools read_xifti is.xifti write_cifti
#'
#' @return A list with \code{"template_mean"} and \code{"template_var"}, as well
#'  as the \code{var_decomp}, \code{mask}, and \code{params}. The dual 
#'  regression results are included too if \code{keep_DR}.
#'
#' @export
#'
estimate_template <- function(
  BOLD, BOLD2=NULL,
  GICA, inds=NULL,
  scale=c("global", "local", "none"), scale_sm_FWHM=2, 
  detrend_DCT=0,
  center_Bcols=FALSE, normA=FALSE,
  Q2=0, Q2_max=NULL,
  brainstructures=c("left","right"), mask=NULL,
  var_method=c("non-negative", "unbiased"),
  keep_DR=FALSE,
  out_fname=NULL,
  FC=FALSE,
  varTol=1e-6, maskTol=.1, missingTol=.1,
  verbose=TRUE) {

  # Check arguments ------------------------------------------------------------

  # Simple argument checks.
  stopifnot(is.logical(center_Bcols) && length(center_Bcols)==1)
  if (is.null(scale) || isFALSE(scale)) { scale <- "none" }
  if (isTRUE(scale)) { 
    warning(
      "Setting `scale='global'`. Use `'global'` or `'local'` ",
      "instead of `TRUE`, which has been deprecated."
    )
    scale <- "global"
  }
  scale <- match.arg(scale, c("global", "local", "none"))
  stopifnot(is.numeric(scale_sm_FWHM) && length(scale_sm_FWHM)==1)
  if (isFALSE(detrend_DCT)) { detrend_DCT <- 0 }
  stopifnot(is.numeric(detrend_DCT) && length(detrend_DCT)==1)
  stopifnot(detrend_DCT >=0 && detrend_DCT==round(detrend_DCT))
  stopifnot(is.logical(normA) && length(normA)==1)
  var_method <- match.arg(var_method, c("non-negative", "unbiased"))
  # if (var_method=="both") { stop("`var_method=='both' not supported yet.") }
  var_name <- c(`non-negative`="varNN", unbiased="varUB")[var_method]
  var_name_alt <- c(`non-negative`="varUB", unbiased="varNN")[var_method]
  if (!is.null(Q2)) { stopifnot(Q2 >= 0) } # Q2_max checked later.
  stopifnot(is.logical(FC) && length(FC)==1)
  stopifnot(is.numeric(maskTol) && length(maskTol)==1 && maskTol >= 0)
  stopifnot(is.numeric(missingTol) && length(missingTol)==1 && missingTol >= 0)
  stopifnot(is.logical(verbose) && length(verbose)==1)
  real_retest <- !is.null(BOLD2)

  # `keep_DR`
  if (is.logical(keep_DR)) {
    stopifnot(length(keep_DR)==1)
  } else {
    if (is.character(keep_DR)) {
      stopifnot(length(keep_DR)==1)
      if (!dir.exists(dirname(keep_DR))) { stop('Directory part of `keep_DR` does not exist.') }
      if (!endsWith(keep_DR, ".rds")) { keep_DR <- paste0(keep_DR, ".rds") }
    } else if (is.list(keep_DR)) {
      stop("Not supported: `keep_DR` must be `TRUE`, `FALSE`, or a single file path.")
      # [TO DO]
      # if (length(keep_DR) != 2) {
      #   stop("If `keep_DR` is a list it must have two entries, each being a vector of file paths the same length as `BOLD`.")
      # }
      # if (length(keep_DR[[1]]) != nN || length(keep_DR[[2]]) != nN) {
      #   stop("If `keep_DR` is a list it must have two entries, each being a vector of file paths the same length as `BOLD`.")
      # }
      # if (!all(dir.exists(dirname(do.call(c, keep_DR))))) { stop('At least one directory part of `keep_DR` does not exist.') }
    }
  }

  # `xii1` will be used to format the output
  xii1 <- NULL

  # `BOLD` and `BOLD2` ---------------------------------------------------------
  # Determine the format of `BOLD` and `BOLD2`.
  format <- infer_BOLD_format(BOLD)
  if (real_retest) {
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

  if (FORMAT == "NIFTI") {
    if (!requireNamespace("RNifti", quietly = TRUE)) {
      stop("Package \"RNifti\" needed to read NIFTI data. Please install it.", call. = FALSE)
    }
  }

  # Ensure `BOLD2` is the same length.
  if (real_retest) {
    if (length(BOLD) != length(BOLD2)) {
      stop("`BOLD2` represents corresponding retest data for `BOLD`, so it must have the same length as `BOLD`.")
    }
  }

  # If BOLD (and BOLD2) is a CIFTI or NIFTI file, check that the file paths exist.
  if (format %in% c("CIFTI", "NIFTI")) {
    missing_BOLD <- !file.exists(BOLD)
    if (all(missing_BOLD)) stop('The files in `BOLD` do not exist.')
    if (real_retest) {
      missing_BOLD2 <- !file.exists(BOLD2)
      if (all(missing_BOLD2)) stop('The files in `BOLD2` do not exist.')
      # determine pairs with at least one missing scan
      missing_BOLD <- missing_BOLD | missing_BOLD2
      rm(missing_BOLD2)
      if (all(missing_BOLD)) stop('Files in `BOLD` and/or `BOLD2` are missing such that no complete pair of data exists.')
    }
    if (any(missing_BOLD)) {
      if (real_retest) {
        warning('There are ', missing_BOLD, ' pairs of `BOLD` and `BOLD2` with at least one non-existent scan. These pairs will be excluded from template estimation.')
        BOLD <- BOLD[!missing_BOLD]
        BOLD2 <- BOLD[!missing_BOLD]
      } else {
        warning('There are ', missing_BOLD, ' scans in `BOLD` that do not exist. These scans will be excluded from template estimation.')
        BOLD <- BOLD[!missing_BOLD]
      }
    }
  }

  # `out_fname` ----------------------------------------------------------------
  if (!is.null(out_fname)) {
    out_fname <- as.character(out_fname)
    if (!all(dir.exists(dirname(out_fname)))) { stop('Directory part of `out_fname` does not exist.') }
    if (length(out_fname) == 1) {
      if (!endsWith(out_fname, FORMAT_extn)) { out_fname <- paste0(out_fname, FORMAT_extn) }
      out_fname <- c(
        gsub(FORMAT_extn, paste0("_mean", FORMAT_extn), out_fname),
        gsub(FORMAT_extn, paste0("_var", FORMAT_extn), out_fname),
        gsub(FORMAT_extn, paste0("_varDecomp.rds"), out_fname)
      )
    } else if (length(out_fname) == 3) {
      if (!all(endsWith(out_fname[seq(2)], FORMAT_extn))) {
        out_fname[seq(2)] <- paste0(out_fname[seq(2)], FORMAT_extn)
      }
      if (!endsWith(out_fname[3], ".rds")) {
        out_fname[3] <- paste0(out_fname[3], ".rds")
      }
    } else {
      stop(
        "`out_fname` should be a length 1 or 3 character vector giving the ",
        "names for:\n\tThe mean template,\n\tThe variance template,",
        "\n\tand the variance decomposition.\n"
      )
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
    if (is.character(GICA)) { GICA <- RNifti::readNifti(GICA) }
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

  # [TO DO]: NA in GICA?

  # `mask` ---------------------------------------------------------------------
  # Get `mask` as a logical array.
  # Check `GICA` and `mask` dimensions match.
  # Vectorize `GICA`.
  if (FORMAT == "NIFTI") {
    if (is.null(mask)) { stop("`mask` is required.") }
    if (is.character(mask)) { mask <- RNifti::readNifti(mask) }
    if (dim(mask)[length(dim(mask))] == 1) { mask <- array(mask, dim=dim(mask)[length(dim(mask))-1]) }
    if (is.numeric(mask)) {
      cat("Coercing `mask` to a logical array with `as.logical`.\n")
      mask[] <- as.logical(mask)
    }
    nI <- dim(mask); nV <- sum(mask)
    stopifnot(length(dim(GICA)) %in% c(2, length(nI)+1))
    if (length(dim(GICA)) == length(nI)+1) {
      if (length(dim(GICA)) != 2) {
        stopifnot(all(dim(GICA)[seq(length(dim(GICA))-1)] == nI))
      }
      if (all(dim(GICA)[seq(length(dim(GICA))-1)] == nI)) {
        GICA <- matrix(GICA[rep(mask, nQ)], ncol=nQ)
        stopifnot(nrow(GICA) == nV)
      }
    }
  } else {
    nI <- nV <- nrow(GICA)
  }

  # Center `GICA` columns.
  center_Gcols <- TRUE
  if (center_Gcols) { GICA <- colCenter(GICA) }

  # Process each scan ----------------------------------------------------------
  if (verbose) {
    cat("Data input format:             ", format, "\n")
    cat('Number of data locations:      ', nV, "\n")
    if (FORMAT == "NIFTI") {
      cat("Unmasked dimensions:           ", paste(nI, collapse=" x "), "\n")
    }
    cat('Number of original group ICs:  ', nQ, "\n")
    cat('Number of template ICs:        ', nL, "\n")
    cat('Number of training subjects:   ', nN, "\n")
  }

  nM <- 2
  DR0 <- array(NA, dim=c(nM, nN, nL, nV)) # measurements by subjects by components by locations
  if(FC) FC0 <- array(NA, dim=c(nM, nN, nL, nL)) # for functional connectivity template

  for (ii in seq(nN)) {
    if(verbose) { cat(paste0(
      '\nReading and analyzing data for subject ', ii,' of ', nN, ".\n"
    )) }
    if (real_retest) { B2 <- BOLD2[[ii]] } else { B2 <- NULL }

    DR_ii <- dual_reg2(
      BOLD[[ii]], BOLD2=B2,
      format=format,
      GICA=GICA,
      center_Bcols=center_Bcols,
      scale=scale, scale_sm_FWHM=scale_sm_FWHM, 
      detrend_DCT=detrend_DCT,
      normA=normA,
      Q2=Q2, Q2_max=Q2_max,
      brainstructures=brainstructures, mask=mask,
      varTol=varTol, maskTol=maskTol,
      verbose=verbose
    )

    # Add results if this subject was not skipped.
    # (Subjects are skipped if too many locations are masked out.)
    if (!is.null(DR_ii)) {
      DR0[1,ii,,] <- DR_ii$test[inds,]
      DR0[2,ii,,] <- DR_ii$retest[inds,]
      if(FC) {
        FC0[1,ii,,] <- cov(DR_ii$test[,inds])
        FC0[2,ii,,] <- cov(DR_ii$retest[,inds])
      }
    }
  }
  rm(DR_ii)

  # Aggregate results, and compute templates. ----------------------------------

  # Mask out locations for which too many subjects do not have data
  nVm <- nV # Number of locations after data masking.
  if (missingTol < 1) { missingTol <- missingTol * nN }
  # Assume is.na(DR0[mm,,,]) is the same for all mm
  # Assume is.na(DR0[,,cc,]) is the same for all cc
  NA_counts <- colSums(is.na(DR0[1,,1,]))
  maskAll <- NA_counts < missingTol
  use_mask <- !all(maskAll)
  if (use_mask) {
    if (all(!maskAll)) { stop(
      "No locations meet the minimum number of subjects with data (",
      round(missingTol, digits=3), "). Check `maskTol` and `missingTol`."
    ) }
    DR0 <- DR0[,,,maskAll,drop=FALSE]
    nVm <- sum(maskAll)
  }
  # Note that `NA` values may still exist in `DR0`.

  # Vectorize components and locations
  DR0 <- array(DR0, dim=c(nM, nN, nL*nVm))
  # FC0 <- array(FC1, dim=c(nM, nN, nL*nL))

  if (verbose) { cat("\nEstimating template.\n") }
  # Estimate the mean and variance templates.
  # Also obtain the variance decomposition.
  x <- estimate_template_from_DR(DR0, c(nL, nVm))
  template <- x$template
  var_decomp <- x$var_decomp
  rm(x)

  # Unmask the templates.
  if (use_mask) {
    for (tname in c("mean", "varUB", "varNN")) {
      ttemp <- matrix(NA, nrow=nrow(template[[tname]]), ncol=length(maskAll))
      ttemp[,maskAll] <- template[[tname]]
      template[[tname]] <- ttemp
    }
  }

  # Estimate FC template
  if(FC){

    # [TO DO]: move to a separate function.
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

  # Format and save template ---------------------------------------------------
  # Keep DR
  if (!isFALSE(keep_DR)) {
    DR0 <- array(DR0, dim=c(nM, nN, nL, nVm)) # Undo vectorize
    if (use_mask) {
      DR0temp <- array(NA, dim=c(dim(DR0)[seq(3)], length(maskAll)))
      DR0temp[,,,maskAll] <- DR0
      DR0 <- DR0temp
    }
    if (is.character(keep_DR)) {
      saveRDS(DR0, keep_DR)
      keep_DR <- FALSE # no longer need it.
    } else if (!isTRUE(keep_DR)) {
      warning("`keep_DR` should be `TRUE`, `FALSE`, or a file path. Using `FALSE`.")
      keep_DR <- FALSE
    }
  }
  if (!keep_DR) { rm(DR0) }

  # Params, formatted as length-one character vectors to put in "xifti" metadata
  indsp <- if (all(seq(nQ) %in% inds)) { paste("all", nQ) } else { inds }
  tparams <- list(
    num_subjects=nN, num_visits=nM,
    inds=indsp, center_Bcols=center_Bcols,
    scale=scale, detrend_DCT=detrend_DCT, normA=normA,
    Q2=Q2, Q2_max=Q2_max,
    brainstructures=brainstructures,
    var_method=var_method,
    pseudo_retest=!real_retest
  )
  tparams <- lapply(
    tparams,
    function(x) {
      if (is.null(x)) { x <- "NULL"};
      paste0(as.character(x), collapse=" ")
    }
  )

  # Save
  if (FORMAT == "CIFTI" && !is.null(xii1)) {
    # Format template as "xifti"s
    GICA <- newdata_xifti(select_xifti(xii1, rep(1, nL)), GICA[,inds])
    GICA$meta$cifti$names <- paste0("IC ", inds)
    for (tname in c("mean", "varUB", "varNN")) {
      template[[tname]] <- newdata_xifti(GICA, t(template[[tname]]))
      template[[tname]]$meta$cifti$misc <- c(
        list(template=tname), tparams
      )
    }
    # Save
    if (!is.null(out_fname)) {
      if (verbose) { cat("\nWriting result to files.\n") }
      write_cifti(template$mean, out_fname[1], verbose=verbose)
      write_cifti(template[[var_name]], out_fname[2], verbose=verbose)
      saveRDS(var_decomp, out_fname[3])
    }
  } else if (FORMAT == "NIFTI") {
    for (tname in c("mean", "varUB", "varNN")) {
      template[[tname]] <- RNifti::asNifti(
        unmask_subcortex(t(template[[tname]]), mask, fill=NA)
      )
    }
    if (!is.null(out_fname)) {
      if (verbose) { cat("\nWriting result to files.\n") }
      writeNIfTI(template$mean, out_fname[1])
      writeNIfTI(template[[var_name]], out_fname[2])
      saveRDS(var_decomp, out_fname[3])
    }
  } else {
    if (!is.null(out_fname)) {
      if (verbose) { cat("\nWriting result to files.\n") }
      saveRDS(template$mean, out_fname[1])
      saveRDS(template$var, out_fname[2])
      saveRDS(var_decomp, out_fname[3])
    }
  }

  # If masking was performed, and if verbose, report the NA count in templates.
  # If no masking was performed, remove unnecessary metadata.
  if (use_mask) {
    if (verbose) {
      cat("Number of template locations with missing values:", sum(!maskAll), "\n")
    }
  } else {
    maskAll <- NULL
  }

  # Results list.
  result <- list(
    template_mean=template$mean,
    template_var=template[[var_name]],
    # template_FC=template$FC,
    var_decomp=var_decomp,
    mask=maskAll,
    params=tparams
  )

  # Add paths to written files if applicable.
  if (!is.null(out_fname)) { result$saved_files <- out_fname }
  # Add DR if applicable.
  if (keep_DR) { result$DR <- DR0 }
  # Add mask if applicable
  if (FORMAT == "NIFTI") { result$mask <- mask }

  # Return results.
  class(result) <- paste0("template.", tolower(FORMAT))
  result
}

#' @rdname estimate_template
#'
estimate_template.cifti <- function(
  BOLD, BOLD2=NULL,
  GICA, inds=NULL,
  scale=c("global", "local", "none"), scale_sm_FWHM=2, 
  detrend_DCT=0,
  center_Bcols=FALSE, normA=FALSE,
  Q2=0, Q2_max=NULL,
  brainstructures=c("left","right"),
  var_method=c("non-negative", "unbiased"),
  keep_DR=FALSE,
  out_fname=NULL,
  FC=FALSE,
  varTol=1e-6, maskTol=.1, missingTol=.1,
  verbose=TRUE) {

  estimate_template(
    BOLD=BOLD, BOLD2=BOLD2,
    GICA=GICA, inds=inds,
    scale=scale, scale_sm_FWHM=scale_sm_FWHM,
    detrend_DCT=detrend_DCT, 
    center_Bcols=center_Bcols, normA=normA,
    Q2=Q2, Q2_max=Q2_max, 
    brainstructures=brainstructures,
    var_method=var_method,
    keep_DR=keep_DR,
    out_fname=out_fname,
    FC=FC,
    varTol=varTol, maskTol=maskTol, missingTol=missingTol,
    verbose=verbose
  )
}

#' @rdname estimate_template
#'
estimate_template.nifti <- function(
  BOLD, BOLD2=NULL,
  GICA, inds=NULL,
  scale=c("global", "local", "none"), scale_sm_FWHM=2, 
  detrend_DCT=0,
  center_Bcols=FALSE, normA=FALSE,
  Q2=0, Q2_max=NULL,
  mask=NULL,
  var_method=c("non-negative", "unbiased"),
  keep_DR=FALSE,
  out_fname=NULL,
  FC=FALSE,
  varTol=1e-6, maskTol=.1, missingTol=.1,
  verbose=TRUE) {

  estimate_template(
    BOLD=BOLD, BOLD2=BOLD2,
    GICA=GICA, inds=inds,
    scale=scale, scale_sm_FWHM=scale_sm_FWHM,
    detrend_DCT=detrend_DCT, 
    center_Bcols=center_Bcols, normA=normA,
    Q2=Q2, Q2_max=Q2_max, 
    mask=mask,
    var_method=var_method,
    keep_DR=keep_DR,
    out_fname=out_fname,
    FC=FC,
    varTol=varTol, maskTol=maskTol, missingTol=missingTol,
    verbose=verbose
  )
}
