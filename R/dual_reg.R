#' Dual Regression
#'
#' @param BOLD Subject-level fMRI data matrix (\eqn{V \times T}). Rows will be
#'  centered. 
#' @param GICA Group-level independent components (\eqn{V \times Q})
#' @param center_Bcols Center BOLD across columns (each image)? Default: \code{FALSE}
#'  (not recommended).
#' @param scale A logical value indicating whether the fMRI timeseries should be
#'  scaled by the image standard deviation. Default: \code{TRUE}.
#' @param detrend_DCT Detrend the data? This is an integer number of DCT bases 
#'  to use for detrending. If \code{0} (default), do not detrend.
#' @param normA Scale each IC timeseries (column of \eqn{A}) in the dual 
#'  regression estimates? Default: \code{FALSE} (not recommended). Note that the
#'  product \eqn{A \times S} remains the same with either option.
#'
#' @importFrom matrixStats colVars
#' 
#' @return A list containing 
#'  the subject-level independent components \strong{S} (\eqn{V \times Q}), 
#'  and subject-level mixing matrix \strong{A} (\eqn{TxQ}).
#' 
#' @export
#'
dual_reg <- function(
  BOLD, GICA, scale=TRUE, detrend_DCT=0, center_Bcols=FALSE, normA=FALSE){

  stopifnot(is.matrix(BOLD))
  stopifnot(is.matrix(GICA))
  stopifnot(is.logical(scale) && length(scale)==1)
  stopifnot(is.logical(normA) && length(normA)==1)

  nV <- nrow(BOLD) #number of data locations
  nT <- ncol(BOLD) #length of timeseries
  if(nV < nT) warning('More time points than voxels. Are you sure?')
  if(nV != nrow(GICA)) {
    stop('The number of voxels in dat (', nV, ') and GICA (', nrow(GICA), ') must match')
  }

  nQ <- ncol(GICA) #number of ICs
  if(nQ > nV) warning('More ICs than voxels. Are you sure?')
  if(nQ > nT) warning('More ICs than time points. Are you sure?')

  # Center each voxel timecourse. Do not center the image at each timepoint.
  # Standardize scale if `scale`, and detrend if `detrend_DCT`.
  # Transpose it: now `BOLD` is TxV.
  BOLD <- t(norm_BOLD(
    BOLD, center_rows=TRUE, center_cols=center_Bcols, 
    scale=scale, detrend_DCT=detrend_DCT
  ))

  # Center each group IC across space. (Used to be a function argument.)
  center_Gcols <- TRUE
  if (center_Gcols) { GICA <- colCenter(GICA) }

  # Estimate A (IC timeseries).
  # We need to center `BOLD` across space because the linear model has no intercept.
  A <- ((BOLD - rowMeans(BOLD)) %*% GICA) %*% chol2inv(chol(crossprod(GICA)))

  # Center each subject IC timecourse across time.
  # (Redundant. Since BOLD is column-centered, A is already column-centered.)
  # A <- colCenter(A)

  # Normalize each subject IC timecourse if `normA`.
  if (normA) { A <- scale(A) }

  # Estimate S (IC maps).
  # No worry about the intercept: `BOLD` and `A` are centered across time.
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
#' 
#'  If \code{BOLD2} is not provided, \code{BOLD} will be split in half; 
#'  the first half will be the test data and the second half will be the retest data.
#' @param GICA Group ICA maps as a (vectorized) numeric matrix 
#'  (\eqn{V \times Q}). Its columns will be centered.
#' @param scale A logical value indicating whether the fMRI timeseries should be
#'  scaled by the image standard deviation. Default: \code{TRUE}.
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
#' @param format Expected format of \code{BOLD} and \code{BOLD2}. Should be one 
#'  of the following: a \code{"CIFTI"} file path, a \code{"xifti"} object, a
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
#'  than \eqn{T * .75 - Q} where \eqn{T} is the minimum number of timepoints in
#'  each fMRI scan and \eqn{Q} is the number of group ICs. If \code{NULL} 
#'  (default), \code{Q2_max} will be set to \eqn{T * .50 - Q}, rounded.
#' @param verbose Display progress updates? Default: \code{TRUE}.
#' 
#' @return The dual regression S matrices.
#' 
#' @keywords internal
dual_reg2 <- function(
  BOLD, BOLD2=NULL, 
  format=c("CIFTI", "xifti", "NIFTI", "nifti", "data"),
  GICA, scale=TRUE, detrend_DCT=0, 
  center_Bcols=FALSE, normA=FALSE,
  Q2=0, Q2_max=NULL, 
  brainstructures=c("left", "right"), mask=NULL, 
  verbose=TRUE){

  # Prepare output.
  out <- list(test = NULL, retest = NULL)

  # Load helper variables.
  retest <- !is.null(BOLD2)
  FORMAT <- switch(format,
    CIFTI = "CIFTI",
    xifti = "CIFTI",
    NIFTI = "NIFTI",
    nifti = "NIFTI",
    data = "DATA"
  )
  nQ <- ncol(GICA)
  nI <- nV <- nrow(GICA)

  # Get `BOLD` (and `BOLD2`) as a data matrix or array.  -----------------------
  if (FORMAT == "CIFTI") {
    if (is.character(BOLD)) { BOLD <- ciftiTools::read_cifti(BOLD, brainstructures=brainstructures) }
    if (is.xifti(BOLD)) { BOLD <- as.matrix(BOLD) }
    stopifnot(is.matrix(BOLD))
    if (retest) {
      BOLD2 <- ciftiTools::read_cifti(BOLD2, brainstructures=brainstructures)
      if (is.xifti(BOLD2)) { BOLD2 <- as.matrix(BOLD2) }
      stopifnot(is.matrix(BOLD2))
    }
  } else if (format == "NIFTI") {
    if (is.character(BOLD)) { BOLD <- oro.nifti::readNIfTI(BOLD, reorient=FALSE) }
    stopifnot(length(dim(BOLD)) > 1)
    if (retest) {
      BOLD2 <- oro.nifti::readNIfTI(BOLD2, reorient=FALSE)
      stopifnot(length(dim(BOLD2)) > 1)
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

  # Vectorize `BOLD` (and `BOLD2`). --------------------------------------------
  if (FORMAT=="NIFTI") {
    BOLD <- matrix(BOLD[rep(mask, dBOLD[ldB])], ncol=nT)
    stopifnot(nrow(BOLD) == nV)
    if (retest) { 
      BOLD2 <- matrix(BOLD2[rep(mask, dBOLD[ldB])], ncol=nT)
      stopifnot(nrow(BOLD2) == nV)
    }
  }

  # [TO DO]: Check for `NA` values?

  # Normalize BOLD (and BOLD2) -------------------------------------------------
  # (Center, scale, and detrend)
  BOLD <- norm_BOLD(
    BOLD, center_rows=TRUE, center_cols=center_Bcols, 
    scale=scale, detrend_DCT=detrend_DCT
  )
  if (retest) {
    BOLD2 <- norm_BOLD(
      BOLD2, center_rows=TRUE, center_cols=center_Bcols, 
      scale=scale, detrend_DCT=detrend_DCT
    )
  }

  # Perform dual regression on test and retest data. ---------------------------
  # If no retest data, halve the test data temporarily.
  if (!retest) {
    part1 <- seq(round(nT/2))
    part2 <- setdiff(seq(nT), part1)
  }

  out$test <- dual_reg(
    if (retest) { BOLD } else { BOLD[, part1, drop=FALSE] }, 
    GICA, scale=FALSE, 
    center_Bcols=FALSE, detrend_DCT=0, normA=normA
  )
  out$retest <- dual_reg(
    if (retest) { BOLD2 } else { BOLD[, part2, drop=FALSE] }, 
    GICA, scale=FALSE, 
    center_Bcols=FALSE, detrend_DCT=0, normA=normA
  )

  # Estimate and deal with nuisance ICs. ---------------------------------------
  # (If !retest, we prefer to estimate nuisance ICs across the full scan.)
  BOLD <- rm_nuisIC(BOLD, DR=out$test, Q2=Q2, Q2_max=Q2_max, verbose=verbose)
  if (retest) {
    BOLD2 <- rm_nuisIC(BOLD2, DR=out$retest, Q2=Q2, Q2_max=Q2_max, verbose=verbose)
  }

  # If no retest data, halve the test data. ------------------------------------
  if (!retest) {
    BOLD2 <- BOLD[, part2, drop=FALSE]
    BOLD <- BOLD[, part1, drop=FALSE]
  }

  # Center and scale `BOLD` (and `BOLD2`) again, but do not detrend again. -----
  BOLD <- norm_BOLD(
    BOLD, center_rows=TRUE, center_cols=center_Bcols,
    scale=scale, detrend_DCT=0
  )
  BOLD2 <- norm_BOLD(
    BOLD2, center_rows=TRUE, center_cols=center_Bcols,
    scale=scale, detrend_DCT=0
  )

  # Do DR again. ---------------------------------------------------------------
  out$test_preclean <- out$test$S
  out$test <- dual_reg(
    BOLD, GICA, scale=FALSE, 
    center_Bcols=FALSE, detrend_DCT=0, normA=normA
  )$S
  out$retest_preclean <- out$retest$S
  out$retest <- dual_reg(
    BOLD2, GICA, scale=FALSE, 
    center_Bcols=FALSE, detrend_DCT=0, normA=normA
  )$S

  out
}