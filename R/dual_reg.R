#' Dual Regression
#'
#' @param dat Subject-level fMRI data (\eqn{VxT})
#' @param GICA Group-level independent components (\eqn{VxQ})
#' @param center,scale A logical value indicating whether the fMRI timeseries should
#'  be centered and/or scaled. See \code{\link{scale_BOLD}}.
#' @param normA Normalize the A matrix (spatial maps)?
#'
#' @importFrom matrixStats colVars
#' 
#' @return A list containing the subject-level independent components \strong{S} (\eqn{VxQ}), 
#'  and subject-level mixing matrix \strong{A} (\eqn{TxQ}).
#' 
#' @export
#'
dual_reg <- function(dat, GICA, center=TRUE, scale=FALSE, normA=FALSE){

  nT <- ncol(dat) #length of timeseries
  nvox <- nrow(dat) #number of data locations
  if(nT > nvox) warning('More time points than voxels. Are you sure?')
  if(nvox != nrow(GICA)) {
    stop('The number of voxels in dat (', nvox, ') and GICA (', nrow(GICA), ') must match')
  }

  Q <- ncol(GICA) #number of ICs
  if(Q > nvox) warning('More ICs than voxels. Are you sure?')
  if(Q > nT) warning('More ICs than time points. Are you sure?')

  # center timeseries data across space and time if `center`
  # standardize data scale if `scale`
  # transpose it
  dat <- t(scale_BOLD(dat, center=center, scale=scale))

  #center each group IC over voxels
  GICA - rep(colMeans(GICA), rep.int(nvox, Q))

	#estimate A (IC timeseries)
	A <- (dat %*% GICA) %*% chol2inv(chol(crossprod(GICA)))
	if (normA) { A <- scale(A) }

	#estimate S (IC maps)
	S <- solve(a=crossprod(A), b=crossprod(A, dat))

	#return result
	list(S = S, A = A)
}

#' Dual Regression wrapper
#' 
#' Wrapper to \code{dual_reg}. Handles NIFTI or CIFTI file input, halves the 
#' data in absense of retest data, and may remove nuisance regressors.
#' 
#' @param BOLD subject-level fMRI timeseries data in one of the following formats: 
#'  a CIFTI file path, a \code{"xifti"}, a NIFTI file path, a \code{"nifti"}, or
#'  a \eqn{V \times T} numeric matrix.
#' 
#'  If \code{BOLD2} is provided it must be in the same format as \code{BOLD}; 
#'  \code{BOLD} will be the test data and \code{BOLD2} will be the retest data. 
#'  If \code{BOLD2} is not provided, \code{BOLD} will be split in half; 
#'  the first half will be the test data and the second half will be the retest data.
#' @param GICA Group ICA maps in a format compatible with \code{BOLD}. Can also be
#'  vectorized (a \eqn{V \times T} numeric matrix) no matter the format of \code{BOLD}.
#' @param scale Scale \code{BOLD} (and \code{BOLD2}) by its mean spatial
#'  standard deviation before computing dual regression? Default: \code{TRUE}.
#' @param normA Normalize the A matrix (spatial maps)?
#' @param format Expected format of \code{BOLD} and \code{BOLD2}. Should be one of the following:
#'  \code{"infer"} (default), a \code{"CIFTI"} file, a \code{"xifti"} object, a
#'  \code{"NIFTI"} file, a \code{"nifti"} object, or a \code{"data"} matrix.
#' @param expectedV The number of rows to expect in the data. Or if the fMRI data
#'  is a NIFTI or "nifti", this is the number of in-mask locations.
#' @param brainstructures Only applies if \code{format} is \code{"CIFTI"}. 
#'  Character vector indicating which brain structure(s)
#'  to obtain: \code{"left"} (left cortical surface), \code{"right"} (right
#'  cortical surface) and/or \code{"subcortical"} (subcortical and cerebellar
#'  gray matter). Can also be \code{"all"} (obtain all three brain structures).
#'  Default: \code{c("left","right")} (cortical surface only).
#' @param mask Required if and only if the entries of \code{BOLD} are NIFTI file paths or
#'  \code{"nifti"} objects. A binary brain map of the same size as the fMRI data, with
#'  \code{1} corresponding to in-mask voxels.
#' @param Q2,maxQ Denoise the dual regression estimate? Denoising is based on modeling and
#'  removing nuisance ICs. It may result in a cleaner estimate, but will take longer to compute. 
#'  If both are \code{NULL},
#'  denoising will be performed, with the number of nuisance ICs estimated for \code{BOLD} and \code{BOLD2}
#'  separately. If \code{Q2==0}, do not denoise (default). Otherwise, specify one or the other:
#'  use \code{Q2} to specify the number of nuisance ICs, or \code{maxQ} to specify the number of
#'  total ICs (template/group + nuisance). \eqn{L <= (L+Q2) = maxQ <= T}, where \eqn{L} is the number
#'  of template ICs and \eqn{T} is the number of timepoints in each fMRI scan. 
#' @param verbose Display progress updates? Default: \code{TRUE}.
#' 
#' @keywords internal
dual_reg2 <- function(
  BOLD, BOLD2=NULL, 
  GICA, 
  scale=TRUE, normA=FALSE,
  format=c("infer", "CIFTI", "xifti", "NIFTI", "nifti", "data"), expectedV=NULL,
  brainstructures=c("left", "right"), mask=NULL, 
  Q2=0, maxQ=NULL, 
  verbose=TRUE){

  # TO DO: add removal of nuisance ICs.
  # If !retest, must be careful about the number of ICs since nT is halved.

  # Initialize the return value
  out <- list(
    missing = NULL,
    test = NULL,
    retest = NULL
  )

  # Retest or pseudo retest?
  retest <- !is.null(BOLD2)

  # Determine the format of `BOLD` and `BOLD2`
  format <- match.arg(format, c("infer", "CIFTI", "xifti", "NIFTI", "nifti", "data"))
  if (format == "infer") { format <- infer_BOLD_format(BOLD) }
  if (format %in% c("CIFTI", "NIFTI") && length(BOLD) > 1) { stop("`BOLD` should be length one.") }
  if (retest) {
    format2 <- infer_BOLD_format(BOLD2)
    if (format2 != format) {
      stop("`BOLD` format is ", format, ", but `BOLD2` format is ", format2, ".")
    }
    if (format %in% c("CIFTI", "NIFTI") && length(BOLD2) > 1) { stop("`BOLD2` should be length one.") }
  }
  FORMAT <- switch(format,
    CIFTI = "CIFTI",
    xifti = "CIFTI",
    NIFTI = "NIFTI",
    nifti = "NIFTI",
    data = "DATA"
  )

  # If BOLD (and BOLD2) is a CIFTI or NIFTI, check that it exists.
  if (format %in% c("CIFTI", "NIFTI")) {
    BOLD_exists <- file.exists(c(BOLD, BOLD2))
    if (!any(BOLD_exists)) {
      if (verbose) { warning("Missing BOLD file(s).") }
      out$missing <- c(BOLD, BOLD2)[BOLD_exists]
      return(out)
    }
  }

  # If BOLD (and BOLD2) is a CIFTI or NIFTI, read it in. Also read in the mask.
  # Get CIFTI data as matrix.
  #   CIFTI
  if (format == "CIFTI") {
    # [TO DO]: check for ciftiTools
    BOLD <- as.matrix(read_cifti(BOLD, brainstructures=brainstructures))
    if (retest) {
      BOLD2 <- as.matrix(read_cifti(BOLD2, brainstructures=brainstructures))
    }
  #   NIFTI
  } else if (format == "NIFTI") {
    BOLD <- oro.nifti::readNIfTI(BOLD, reorient=FALSE)
    if (retest) {
      BOLD2 <- oro.nifti::readNIfTI(BOLD2, reorient=FALSE)
    }
    if (is.null(mask)) { stop("`mask` is required.") }
    if (is.character(mask)) { mask <- oro.nifti::readNIfTI(mask) }
  }

  dBOLD <- BOLD # ???

  # Check dims.
  nV <- dim(GICA)[1]
  nQ <- dim(GICA)[2]
  nDV <- ifelse(FORMAT=="NIFTI", sum(mask), nV)
  if (!is.null(expectedV)) {
    expectedV <- as.numeric(expectedV)
    if (length(expectedV) > 1) {
      warning("`expectedV` should be the number of data locations. Using the first entry only.\n")
      expectedV <- expectedV[1]
    }
    if (format == "NIFTI") {
      if (length(dim(mask)) == 4 && dim(mask)[4] == 1) { mask <- mask[,,,1] }
      if (length(dim(mask)) != 3) { stop("`mask` should be a 3D binary image.") }
      if (!is.logical(mask)) {
        if (verbose) { cat("Coercing `mask` to a logical array with `as.logical`.\n") }
        mask[] <- as.logical(mask)
      }
      if (sum(mask) != expectedV) {
        stop("sum(mask) is ", sum(mask), ", but `expectedV` is ", expectedV, ".")
      }
      if (!all(dim(mask)[seq(3)] == dim(BOLD[seq(3)]))) {
        stop("`BOLD` and `mask` do not have the same spatial dimensions (first three dims).")
      }
      if (retest && !all(dim(mask)[seq(3)] == dim(BOLD2[seq(3)]))) {
        stop("`BOLD2` and `mask` do not have the same spatial dimensions (first three dims).")
      }
    } else {
      if (dBOLD[1] != expectedV) {
        stop("nrow(BOLD) is ", dBOLD[1], ", but `expectedV` is ", expectedV, ".")
      }
      if (retest && dim(BOLD2)[1] != expectedV) {
        stop("nrow(BOLD2) is ", dBOLD[1], ", but `nrow(BOLD)` is ", expectedV, ".")
      }
    }
  }

  # Get NIFTI data as matrix.
  if (format == "NIFTI") {
    BOLD <- matrix(BOLD[rep(mask, nQ)], ncol=nQ)
    if (retest) { BOLD2 <- matrix(BOLD2[rep(mask, nQ), ncol=nQ]) }
  }

  # If no retest data, halve the test data
  if (!retest) {
    nT <- dBOLD[length(dBOLD)]
    part1 <- seq(round(nT/2))
    part2 <- setdiff(seq(nT), part1)
    if (format=="CIFTI") {
      BOLD2 <- BOLD[, part2, drop=FALSE]
      BOLD <- BOLD[, part1, drop=FALSE]
    } else if (format=="NIFTI") {
      BOLD2 <- BOLD[,,, part2, drop=FALSE]
      BOLD <- BOLD[,,, part1, drop=FALSE]
    }
  }

  # [TO DO]: Remove NA values?

  # Perform dual regression on test and retest data
  out$test <- dual_reg(BOLD, as.matrix(GICA), scale=scale)
  out$retest <- dual_reg(BOLD2, as.matrix(GICA), scale=scale)
  out$missing <- FALSE

  nT <- min(dim(BOLD)[length(dim(BOLD))], dim(BOLD2)[length(dim(BOLD2))])
  skip_rm_nuisIC <- (!is.null(Q2) && Q2 <= 0) || (!is.null(maxQ) && maxQ <=0)
  if (!skip_rm_nuisIC) {

    L <- ncol(GICA)

    # Check `maxQ`
    if (!is.null(maxQ)) { 
      if(round(maxQ) != maxQ || maxQ <= 0) stop('maxQ must be NULL or a round positive number.')
    } else {
      maxQ <- round(nT/2)
    }
    if (maxQ < L) {
      warning('maxQ must be at least L.  Setting maxQ=L.')
      maxQ <- L
    }
    # This is to avoid the area of the pesel objective function that spikes close 
    #   to rank(X), which often leads to nPC close to rank(X)
    if (maxQ > nT*0.75) {
      warning('maxQ too high, setting to 75% of T.')
      maxQ <- round(nT*0.75)
    }

    # Scale
    BOLD <- scale_BOLD(BOLD, scale=scale)
    BOLD2 <- scale_BOLD(BOLD2, scale=scale)

    # Estimate and deal with nuisance ICs
    if (maxQ > L) {
      BOLD <- rm_nuisIC(BOLD, DR=out$test, Q2=Q2, Q2_max=maxQ-L, verbose=verbose)
      BOLD2 <- rm_nuisIC(BOLD2, DR=out$retest, Q2=Q2, Q2_max=maxQ-L, verbose=verbose)
    }

    # Do DR again.
    out$test_preclean <- out$test$S # [TO DO]: Rename?
    out$test <- dual_reg(BOLD, as.matrix(GICA), scale=scale)
    out$retest_preclean <- out$retest$S
    out$retest <- dual_reg(BOLD2, as.matrix(GICA), scale=scale)
  }

  out$test <- out$test$S
  out$retest <- out$retest$S
  out
}