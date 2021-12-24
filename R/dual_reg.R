#' Dual Regression
#'
#' @param BOLD Subject-level fMRI data matrix (\eqn{V \times T})
#' @param GICA Group-level independent components (\eqn{V \times Q})
#' @param center_rows,center_cols Center BOLD data across rows (each data location's time series) or columns (each time point's image)? Default: \code{TRUE} for both.
#' @param center_Gcols Center GICA across columns (each ICA)? Default: \code{TRUE}.
#' @param scale A logical value indicating whether the fMRI timeseries should be scaled by the image standard deviation.
#' @param detrend_DCT Detrend the data? This is the number of DCT bases to use for detrending. If \code{0} (default), do not detrend.
#' @param normA Scale each IC timeseries (column of \eqn{A}) in the dual regression 
#'  estimates? Default: \code{FALSE}. (The opposite scaling will be applied to \eqn{S}
#'  such that the product \eqn{A \times S} remains the same).
#'
#' @importFrom matrixStats colVars
#' 
#' @return A list containing the subject-level independent components \strong{S} (\eqn{V \times Q}), 
#'  and subject-level mixing matrix \strong{A} (\eqn{TxQ}).
#' 
#' @export
#'
dual_reg <- function(
  BOLD, GICA, 
  center_rows=TRUE, center_cols=TRUE, scale=FALSE, detrend_DCT=0, 
  center_Gcols=TRUE, normA=FALSE){

  stopifnot(is.matrix(BOLD))
  stopifnot(is.matrix(GICA))
  stopifnot(is.logical(center_rows) && length(center_rows)==1)
  stopifnot(is.logical(center_cols) && length(center_cols)==1)
  stopifnot(is.logical(scale) && length(scale)==1)
  stopifnot(is.logical(normA) && length(normA)==1)

  nT <- ncol(BOLD) #length of timeseries
  nV <- nrow(BOLD) #number of data locations
  if(nT > nV) warning('More time points than voxels. Are you sure?')
  if(nV != nrow(GICA)) {
    stop('The number of voxels in dat (', nV, ') and GICA (', nrow(GICA), ') must match')
  }

  nQ <- ncol(GICA) #number of ICs
  if(nQ > nV) warning('More ICs than voxels. Are you sure?')
  if(nQ > nT) warning('More ICs than time points. Are you sure?')

  # Center timeseries data across space and time if `center`
  # Standardize data scale if `scale`
  # Transpose it
  BOLD <- t(norm_BOLD(
    BOLD, center_rows=center_rows, center_cols=center_cols, 
    scale=scale, detrend_DCT=detrend_DCT
  ))

  # Center each group IC over voxels
  if (center_Gcols) { GICA - rep(colMeans(GICA), rep.int(nV, nQ)) }

	# Estimate A (IC timeseries)
	A <- (BOLD %*% GICA) %*% chol2inv(chol(crossprod(GICA)))
	if (normA) { A <- scale(A) }

	# Estimate S (IC maps)
	S <- solve(a=crossprod(A), b=crossprod(A, BOLD))

	#return result
	list(S = S, A = A)
}

#' Dual Regression wrapper
#' 
#' Wrapper to \code{dual_reg} used by `estimate_template`. The format of `BOLD`
#'  (and `BOLD2`) must be provided, and `GICA` must be vectorized if applicable.
#' 
#' @param BOLD,BOLD2 Subject-level fMRI data in one of the following formats: 
#'  a CIFTI file path, a \code{"xifti"} object, a NIFTI file path, a \code{"nifti"} object, or
#'  \eqn{V \times T} numeric matrices, where \eqn{V} is the number of data locations and
#'  \eqn{T} is the number of timepoints.
#' 
#'  If \code{BOLD2} is provided it must be in the same format as \code{BOLD}; 
#'  \code{BOLD} will be the test data and \code{BOLD2} will be the retest data. 
#'  \code{BOLD2} should be the same length as \code{BOLD} and have the same subjects in the same order.
#'  If \code{BOLD2} is not provided, \code{BOLD} will be split in half; 
#'  the first half will be the test data and the second half will be the retest data.
#' @param GICA Group ICA maps in as a (vectorized) numeric matrix 
#'  (\eqn{V \times Q}). Columns should be centered.
#' @param center_rows,center_cols Center BOLD data across rows (each data location's time series) or columns (each time point's image)? Default: \code{TRUE} for both.
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
#' @param format Expected format of \code{BOLD} and \code{BOLD2}. Should be one of the following:
#'  a \code{"CIFTI"} file path, a \code{"xifti"} object, a
#'  \code{"NIFTI"} file path, a \code{"nifti"} object, or a \code{"data"} matrix.
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
#'  than \eqn{T * .75 - Q} where \eqn{T} is the number of timepoints in each 
#'  fMRI scan and \eqn{Q} is the number of group ICs. If \code{NULL} (default),
#'  \code{Q2_max} will be set to \eqn{T * .50 - Q}, rounded.
#' @param verbose Display progress updates? Default: \code{TRUE}.
#' 
#' @keywords internal
dual_reg2 <- function(
  BOLD, BOLD2=NULL, 
  format=c("CIFTI", "xifti", "NIFTI", "nifti", "data"),
  GICA, 
  center_rows=TRUE, center_cols=TRUE, scale=TRUE, detrend_DCT=0, 
  normA=FALSE,
  Q2=0, Q2_max=NULL, 
  brainstructures=c("left", "right"), mask=NULL, 
  verbose=TRUE){

  # Prepare output.
  out <- list(test = NULL, retest = NULL)

  # Load helper variables.
  retest <- is.null(BOLD2)
  FORMAT <- switch(format,
    CIFTI = "CIFTI",
    xifti = "CIFTI",
    NIFTI = "NIFTI",
    nifti = "NIFTI",
    data = "DATA"
  )
  nQ <- ncol(GICA)
  nI <- nV <- nrow(GICA)

  # Get `BOLD` (and `BOLD2`) as a data matrix or array. 
  if (FORMAT == "CIFTI") {
    if (is.character(BOLD)) { BOLD <- ciftiTools::read_cifti(BOLD, brainstructures=brainstructures) }
    if (is.xifti(BOLD)) { BOLD <- as.matrix(BOLD) }
    stopifnot(is.matrix(BOLD))
    if (retest) {
      BOLD2 <- ciftiTools::read_cifti(BOLD2, brainstructures=brainstructures)
      if (is.xifti(BOLD2)) { BOLD <- as.matrix(BOLD2) }
      stopifnot(is.matrix(BOLD2))
    }
  } else if (format == "NIFTI") {
    if (is.character(BOLD)) { BOLD <- oro.nifti::readNIfTI(BOLD, reorient=FALSE) }
    stopifnot(length(dim(BOLD)) > 1)
    if (retest) {
      BOLD2 <- oro.nifti::readNIfTI(BOLD2, reorient=FALSE)
      stopifnot(length(dim(BOLD)) > 1)
    }
  } else {
    stopifnot(is.matrix(BOLD))
    if (retest) { stopifnot(is.matrix(BOLD2)) }
  }
  dBOLD <- dim(BOLD)
  ldB <- length(dim(BOLD))
  nT <- dim(BOLD)[ldB]

  # If `retest`, ensure that spatial dimensions of `BOLD2` match with `BOLD`.
  if (retest) {
    stopifnot(length(dim(BOLD)) == length(dim(BOLD2)))
    stopifnot(all(dBOLD[seq(ldB-1)] == dim(BOLD2)[seq(ldB-1)]))
  }

  # Check BOLD (and BOLD2) dimensions correspond with `GICA` and `mask`.
  stopifnot(ldB-1 == length(nI))
  stopifnot(all(dBOLD[seq(ldB-1)] == nI))

  # Vectorize `BOLD` (and `BOLD2`).
  if (FORMAT=="NIFTI") {
    BOLD <- matrix(BOLD[rep(mask, dBOLD[ldB])], ncol=nT)
    stopifnot(nrow(BOLD) == nV)
    if (retest) { 
      BOLD2 <- matrix(BOLD2[rep(mask, dBOLD[ldB])], ncol=nT)
      stopifnot(nrow(BOLD2) == nV)
    }
  }

  # If no retest data, halve the test data.
  # [TO DO] consider doing this after centering, scaling, detrending, and denoising?
  if (!retest) {
    part1 <- seq(round(nT/2))
    part2 <- setdiff(seq(nT), part1)
    BOLD2 <- BOLD[, part2, drop=FALSE]
    BOLD <- BOLD[, part1, drop=FALSE]
  }

  # [TO DO]: Check for `NA` values?

  # Normalize BOLD (and BOLD2) -------------------------------------------------
  # (Center, scale, and detrend)

  BOLD <- norm_BOLD(
    BOLD, center_rows=center_rows, center_cols=center_cols, 
    scale=scale, detrend_DCT=detrend_DCT
  )
  BOLD2 <- norm_BOLD(
    BOLD2, center_rows=center_rows, center_cols=center_cols, 
    scale=scale, detrend_DCT=detrend_DCT
  )

  # Perform dual regression on test and retest data. ---------------------------
  out$test <- dual_reg(
    BOLD, GICA, center_rows=FALSE, center_cols=FALSE, 
    scale=FALSE, center_Gcols=FALSE, detrend_DCT=0, normA=normA
  )$S
  out$retest <- dual_reg(
    BOLD2, GICA, center_rows=FALSE, center_cols=FALSE, 
    scale=FALSE, center_Gcols=FALSE, detrend_DCT=0, normA=normA
  )$S

  # Return now, if denoising is not needed. ------------------------------------
  nT2 <- min(ncol(BOLD), ncol(BOLD2))
  if (!is.null(Q2) && Q2==0) { return(out) }

  # Estimate and deal with nuisance ICs. ---------------------------------------
  BOLD <- rm_nuisIC(BOLD, DR=out$test, Q2=Q2, Q2_max=Q2_max, verbose=verbose)
  BOLD2 <- rm_nuisIC(BOLD2, DR=out$retest, Q2=Q2, Q2_max=Q2_max, verbose=verbose)

  # Center and scale `BOLD` (and `BOLD2`) again, but do not detrend again. -----
  BOLD <- norm_BOLD(
    BOLD, center_rows=center_rows, center_cols=center_cols, 
    scale=scale, detrend_DCT=0
  )
  BOLD2 <- norm_BOLD(
    BOLD2, center_rows=center_rows, center_cols=center_cols, 
    scale=scale, detrend_DCT=0
  )

  # Do DR again. ---------------------------------------------------------------
  out$test_preclean <- out$test
  out$test <- dual_reg(
    BOLD, GICA, center_rows=FALSE, center_cols=FALSE, 
    scale=FALSE, detrend_DCT=0, normA=normA
  )$S
  out$retest_preclean <- out$retest
  out$retest <- dual_reg(
    BOLD2, GICA, center_rows=FALSE, center_cols=FALSE, 
    scale=FALSE, detrend_DCT=0, normA=normA
  )$S

  out
}